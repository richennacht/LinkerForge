use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::{Path, PathBuf},
    process::Stdio,
    sync::Mutex,
    time::{Duration, SystemTime, UNIX_EPOCH},
};
use tauri::Manager;
use tokio::process::Command;

const MAX_RESIDUES: usize = 2_100;
const MAX_PDB_BYTES: u64 = 15_000_000;

#[derive(Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct FoldStatus {
    phase: String,
    message: String,
    run_dir: Option<String>,
    sequence_sha256: Option<String>,
    pdb: Option<String>,
}

impl Default for FoldStatus {
    fn default() -> Self {
        Self {
            phase: "idle".into(),
            message: "No local prediction started.".into(),
            run_dir: None,
            sequence_sha256: None,
            pdb: None,
        }
    }
}

#[derive(Default)]
struct FoldState(Mutex<FoldStatus>);

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ComputeStatus {
    predictor_available: bool,
    predictor_backend: String,
    gpu: String,
    note: String,
}

fn sha256_bytes(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn validate_sequence(sequence: &str) -> Result<(), String> {
    if sequence.is_empty() || sequence.len() > MAX_RESIDUES {
        return Err(format!(
            "Full fusion must contain 1–{MAX_RESIDUES} residues."
        ));
    }
    if !sequence
        .bytes()
        .all(|aa| b"ACDEFGHIKLMNPQRSTVWY".contains(&aa))
    {
        return Err("Full fusion must contain only uppercase standard amino acids.".into());
    }
    Ok(())
}

fn read_ca_sequence(pdb: &str) -> Result<String, String> {
    let mut sequence = String::new();
    let mut previous = String::new();
    let mut chain_seen: Option<String> = None;
    for line in pdb.lines().filter(|line| line.starts_with("ATOM  ")) {
        if line.len() < 27 || line.get(12..16) != Some(" CA ") {
            continue;
        }
        let chain = line.get(21..22).ok_or("Malformed PDB chain field")?;
        if chain_seen.as_deref().is_some_and(|seen| seen != chain) {
            return Err("Predicted PDB contains multiple protein chains.".into());
        }
        chain_seen = Some(chain.to_string());
        let key = line.get(21..27).ok_or("Malformed PDB residue field")?;
        if key == previous {
            continue;
        }
        previous = key.to_string();
        let code = line.get(17..20).ok_or("Malformed PDB amino acid field")?;
        let aa = match code {
            "ALA" => 'A',
            "CYS" => 'C',
            "ASP" => 'D',
            "GLU" => 'E',
            "PHE" => 'F',
            "GLY" => 'G',
            "HIS" => 'H',
            "ILE" => 'I',
            "LYS" => 'K',
            "LEU" => 'L',
            "MET" => 'M',
            "ASN" => 'N',
            "PRO" => 'P',
            "GLN" => 'Q',
            "ARG" => 'R',
            "SER" => 'S',
            "THR" => 'T',
            "VAL" => 'V',
            "TRP" => 'W',
            "TYR" => 'Y',
            _ => return Err(format!("Unknown PDB residue {code}")),
        };
        sequence.push(aa);
    }
    if sequence.is_empty() {
        return Err("Predicted PDB contains no C-alpha trace.".into());
    }
    Ok(sequence)
}

fn find_matching_pdb(output: &Path, expected: &str) -> Result<(PathBuf, String), String> {
    let entries = fs::read_dir(output).map_err(|error| error.to_string())?;
    let mut files: Vec<PathBuf> = entries
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| {
            path.extension()
                .is_some_and(|ext| ext.to_string_lossy().eq_ignore_ascii_case("pdb"))
        })
        .collect();
    files.sort_by_key(|path| (!path.to_string_lossy().contains("rank_001"), path.clone()));
    for path in files {
        if fs::metadata(&path)
            .map_err(|error| error.to_string())?
            .len()
            > MAX_PDB_BYTES
        {
            continue;
        }
        let pdb = fs::read_to_string(&path).map_err(|error| error.to_string())?;
        if read_ca_sequence(&pdb).ok().as_deref() == Some(expected) {
            return Ok((path, pdb));
        }
    }
    Err(
        "No predicted single-chain PDB matches the entire submitted fusion. Inspect the run log."
            .into(),
    )
}

fn probe_compute() -> ComputeStatus {
    let native_available = std::process::Command::new("where.exe")
        .arg("colabfold_batch")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|status| status.success());
    // WSL's Windows stub can raise a blocking GUI error on machines without WSL.
    // Probe it only after the user explicitly enables this optional backend.
    let wsl_enabled = std::env::var("LINKERFORGE_ENABLE_WSL").as_deref() == Ok("1");
    let wsl_available = wsl_enabled
        && !native_available
        && std::process::Command::new("wsl.exe")
            .args(["--exec", "which", "colabfold_batch"])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|status| status.success());
    let gpu = std::process::Command::new("nvidia-smi")
        .args([
            "--query-gpu=name,memory.total,memory.free",
            "--format=csv,noheader",
        ])
        .output()
        .ok()
        .filter(|result| result.status.success())
        .map(|result| String::from_utf8_lossy(&result.stdout).trim().to_string())
        .filter(|value| !value.is_empty())
        .unwrap_or_else(|| "NVIDIA GPU not detected by nvidia-smi".into());
    ComputeStatus {
        predictor_available: native_available || wsl_available,
        predictor_backend: if native_available { "Windows PATH" } else if wsl_available { "WSL2 PATH" } else { "unavailable" }.into(),
        gpu,
        note: "LocalColabFold must be installed separately on Windows PATH or WSL2 PATH. WSL2 probing is opt-in via LINKERFORGE_ENABLE_WSL=1. GPU memory may be insufficient for a full fusion; no model is bundled.".into(),
    }
}

#[tauri::command]
async fn compute_status() -> Result<ComputeStatus, String> {
    tauri::async_runtime::spawn_blocking(probe_compute)
        .await
        .map_err(|error| error.to_string())
}

#[tauri::command]
fn fold_status(state: tauri::State<'_, FoldState>) -> Result<FoldStatus, String> {
    state
        .0
        .lock()
        .map(|status| status.clone())
        .map_err(|error| error.to_string())
}

#[tauri::command]
async fn start_fold(
    app: tauri::AppHandle,
    state: tauri::State<'_, FoldState>,
    sequence: String,
) -> Result<FoldStatus, String> {
    validate_sequence(&sequence)?;
    let capability = tauri::async_runtime::spawn_blocking(probe_compute)
        .await
        .map_err(|error| error.to_string())?;
    if !capability.predictor_available {
        return Err("ColabFold is unavailable. Install and validate LocalColabFold in Windows or WSL2 before starting a local prediction.".into());
    }
    let digest = sha256_bytes(sequence.as_bytes());
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|error| error.to_string())?
        .as_secs();
    let base = app
        .path()
        .app_local_data_dir()
        .map_err(|error| error.to_string())?;
    let run_dir = base
        .join("folds")
        .join(format!("{seconds}-{}", &digest[..12]));
    fs::create_dir_all(&run_dir).map_err(|error| error.to_string())?;
    let fasta = run_dir.join("fusion.fasta");
    fs::write(&fasta, format!(">LinkerForge_full_fusion\n{sequence}\n"))
        .map_err(|error| error.to_string())?;
    let output = run_dir.join("prediction");
    fs::create_dir(&output).map_err(|error| error.to_string())?;
    let running = FoldStatus {
        phase: "running".into(),
        message: "Predicting the complete fusion locally. This may take a long time.".into(),
        run_dir: Some(run_dir.display().to_string()),
        sequence_sha256: Some(digest.clone()),
        pdb: None,
    };
    {
        let mut slot = state.0.lock().map_err(|error| error.to_string())?;
        if slot.phase == "running" {
            return Err("A local fold is already running. Only one GPU job is allowed.".into());
        }
        *slot = running.clone();
    }
    let app_for_task = app.clone();
    tauri::async_runtime::spawn(async move {
        let log_path = run_dir.join("colabfold.log");
        let result = async {
            let log = fs::File::create(&log_path).map_err(|error| error.to_string())?;
            let errors = log.try_clone().map_err(|error| error.to_string())?;
            let (program, prefix, input, results) = if capability.predictor_backend == "WSL2 PATH" {
                let convert = |path: &Path| -> Result<String, String> {
                    let output = std::process::Command::new("wsl.exe")
                        .arg("--exec")
                        .arg("wslpath")
                        .arg("-a")
                        .arg(path)
                        .output()
                        .map_err(|error| error.to_string())?;
                    if !output.status.success() {
                        return Err("WSL2 could not translate a run path.".into());
                    }
                    Ok(String::from_utf8_lossy(&output.stdout).trim().to_string())
                };
                (
                    "wsl.exe",
                    vec![
                        "--exec",
                        "env",
                        "TF_FORCE_UNIFIED_MEMORY=1",
                        "TF_FORCE_GPU_ALLOW_GROWTH=true",
                        "XLA_PYTHON_CLIENT_MEM_FRACTION=4.0",
                        "XLA_PYTHON_CLIENT_PREALLOCATE=false",
                        "XLA_PYTHON_CLIENT_ALLOCATOR=platform",
                        "colabfold_batch",
                    ],
                    convert(&fasta)?,
                    convert(&output)?,
                )
            } else {
                (
                    "colabfold_batch",
                    Vec::new(),
                    fasta.display().to_string(),
                    output.display().to_string(),
                )
            };
            let mut command = Command::new(program);
            command
                .args(prefix)
                .arg("--use-fast-kernels")
                .arg(input)
                .arg(results)
                .env("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
                .env("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
                .stdout(Stdio::from(log))
                .stderr(Stdio::from(errors))
                .kill_on_drop(true);
            let mut child = command.spawn().map_err(|error| error.to_string())?;
            let exit = tokio::time::timeout(Duration::from_secs(120 * 60), child.wait())
                .await
                .map_err(|_| {
                    "ColabFold exceeded the 120-minute limit and was stopped.".to_string()
                })?
                .map_err(|error| error.to_string())?;
            if !exit.success() {
                return Err(format!(
                    "ColabFold exited with {exit}. Read {}",
                    log_path.display()
                ));
            }
            find_matching_pdb(&output, &sequence)
        }
        .await;
        let final_status = match result {
            Ok((pdb_path, pdb)) => FoldStatus {
                phase: "complete".into(),
                message: format!("Matching full-fusion PDB: {}", pdb_path.display()),
                run_dir: Some(run_dir.display().to_string()),
                sequence_sha256: Some(digest.clone()),
                pdb: Some(pdb),
            },
            Err(message) => FoldStatus {
                phase: "failed".into(),
                message,
                run_dir: Some(run_dir.display().to_string()),
                sequence_sha256: Some(digest.clone()),
                pdb: None,
            },
        };
        let manifest = serde_json::json!({
            "status": final_status.phase,
            "message": final_status.message,
            "predictor": "local ColabFold full-fusion sequence",
            "predictor_backend": capability.predictor_backend,
            "desktop_app_version": env!("CARGO_PKG_VERSION"),
            "sequence_sha256": digest,
            "sequence_length": sequence.len(),
            "selected_pdb_sha256": final_status.pdb.as_ref().map(|pdb| sha256_bytes(pdb.as_bytes())),
            "settings": {"fast_kernels": true, "jax_preallocate": false, "jax_allocator": "platform", "wsl_unified_memory": capability.predictor_backend == "WSL2 PATH", "timeout_minutes": 120},
            "evidence_quality": "predicted_unvalidated",
            "claim": "predicted structure, not experimental validation"
        });
        let _ = fs::write(run_dir.join("run_manifest.json"), format!("{}\n", manifest));
        if let Some(state) = app_for_task.try_state::<FoldState>() {
            if let Ok(mut slot) = state.0.lock() {
                *slot = final_status;
            }
        }
    });
    Ok(running)
}

pub fn run() {
    tauri::Builder::default()
        .manage(FoldState::default())
        .invoke_handler(tauri::generate_handler![
            compute_status,
            fold_status,
            start_fold
        ])
        .run(tauri::generate_context!())
        .expect("failed to run LinkerForge desktop");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequence_gate() {
        assert!(validate_sequence("MGS").is_ok());
        assert!(validate_sequence("MGX").is_err());
        assert!(validate_sequence("").is_err());
    }

    #[test]
    fn pdb_gate_requires_full_trace() {
        let line = "ATOM      1  CA  GLY A   1      0.000   0.000   0.000  1.00 50.00           C";
        assert_eq!(read_ca_sequence(line).unwrap(), "G");
    }
}
