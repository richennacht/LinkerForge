use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    ffi::OsStr,
    fs,
    fs::OpenOptions,
    io::Write,
    path::{Path, PathBuf},
    process::Stdio,
    sync::Mutex,
    time::{Duration, SystemTime, UNIX_EPOCH},
};
use tauri::Manager;
use tokio::process::Command;

#[cfg(windows)]
use std::os::windows::process::CommandExt;

const MAX_RESIDUES: usize = 2_100;
const MAX_PDB_BYTES: u64 = 15_000_000;
#[cfg(windows)]
const CREATE_NO_WINDOW: u32 = 0x08000000;

fn hidden_command<S: AsRef<OsStr>>(program: S) -> std::process::Command {
    let mut command = std::process::Command::new(program);
    #[cfg(windows)]
    command.creation_flags(CREATE_NO_WINDOW);
    command
}

fn hidden_async_command<S: AsRef<OsStr>>(program: S) -> Command {
    let mut command = Command::new(program);
    #[cfg(windows)]
    command.as_std_mut().creation_flags(CREATE_NO_WINDOW);
    command
}

#[derive(Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct FoldStatus {
    phase: String,
    message: String,
    run_dir: Option<String>,
    sequence_sha256: Option<String>,
    pdb: Option<String>,
    metrics: Option<FoldMetrics>,
    artifacts: Option<FoldArtifacts>,
}

#[derive(Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct FoldMetrics {
    predictor: String,
    model_name: Option<String>,
    mean_plddt: Option<f64>,
    ptm: Option<f64>,
    max_pae: Option<f64>,
    mean_interdomain_pae: Option<f64>,
}

#[derive(Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct FoldArtifacts {
    pdb_path: Option<String>,
    score_json_path: Option<String>,
    pae_plot_path: Option<String>,
    plddt_plot_path: Option<String>,
    log_path: Option<String>,
}

impl Default for FoldStatus {
    fn default() -> Self {
        Self {
            phase: "idle".into(),
            message: "No local prediction started.".into(),
            run_dir: None,
            sequence_sha256: None,
            pdb: None,
            metrics: None,
            artifacts: None,
        }
    }
}

#[derive(Default)]
struct FoldState(Mutex<FoldStatus>);

struct GpuJobLock {
    path: PathBuf,
    file: Option<fs::File>,
}

fn read_lock_pid(path: &Path) -> Result<u32, String> {
    let contents = fs::read_to_string(path)
        .map_err(|error| format!("Could not inspect GPU job lock {}: {error}", path.display()))?;
    contents
        .lines()
        .find_map(|line| line.strip_prefix("pid="))
        .ok_or_else(|| {
            format!(
                "GPU job lock {} has no PID and was left in place for safety.",
                path.display()
            )
        })?
        .parse::<u32>()
        .ok()
        .filter(|pid| *pid > 0)
        .ok_or_else(|| {
            format!(
                "GPU job lock {} has an invalid PID and was left in place for safety.",
                path.display()
            )
        })
}

#[cfg(windows)]
fn process_is_active(pid: u32) -> Result<bool, String> {
    let filter = format!("PID eq {pid}");
    let output = hidden_command("tasklist.exe")
        .args(["/FI", &filter, "/FO", "CSV", "/NH"])
        .output()
        .map_err(|error| format!("Could not verify GPU lock PID {pid}: {error}"))?;
    if !output.status.success() {
        return Err(format!(
            "Could not verify GPU lock PID {pid}; the lock was left in place."
        ));
    }
    Ok(String::from_utf8_lossy(&output.stdout).contains(&format!("\"{pid}\"")))
}

#[cfg(not(windows))]
fn process_is_active(pid: u32) -> Result<bool, String> {
    Ok(Path::new("/proc").join(pid.to_string()).exists())
}

#[cfg(windows)]
fn prediction_or_simulation_active() -> Result<bool, String> {
    let script = r#"$active = Get-CimInstance Win32_Process | Where-Object { $_.ProcessId -ne $PID -and (($_.Name -match '^(colabfold|alphafold|gmx|gromacs|rosetta)') -or ($_.CommandLine -match '(?i)(colabfold_batch|run_alphafold|openmm|gromacs|(?:^|[\\/ ])gmx(?:\.exe)?(?: |$)|rosetta)')) } | Select-Object -First 1; if ($null -ne $active) { 'active' } else { 'inactive' }"#;
    let output = hidden_command("powershell.exe")
        .args(["-NoProfile", "-NonInteractive", "-Command", script])
        .output()
        .map_err(|error| {
            format!("Could not inspect prediction and simulation processes: {error}")
        })?;
    if !output.status.success() {
        return Err("Could not verify whether a prediction or simulation is active; the GPU lock was left in place.".into());
    }
    Ok(String::from_utf8_lossy(&output.stdout).trim() == "active")
}

#[cfg(not(windows))]
fn prediction_or_simulation_active() -> Result<bool, String> {
    let needles = [
        "colabfold",
        "alphafold",
        "openmm",
        "gromacs",
        "gmx",
        "rosetta",
    ];
    for entry in fs::read_dir("/proc").map_err(|error| error.to_string())? {
        let entry = entry.map_err(|error| error.to_string())?;
        if entry.file_name().to_string_lossy().parse::<u32>().is_err() {
            continue;
        }
        let command_line = fs::read(entry.path().join("cmdline")).unwrap_or_default();
        let command_line = String::from_utf8_lossy(&command_line).to_ascii_lowercase();
        if needles.iter().any(|needle| command_line.contains(needle)) {
            return Ok(true);
        }
    }
    Ok(false)
}

fn create_gpu_lock(path: &Path) -> Result<fs::File, std::io::Error> {
    OpenOptions::new().write(true).create_new(true).open(path)
}

impl GpuJobLock {
    fn acquire(runs_root: &Path) -> Result<Self, String> {
        fs::create_dir_all(runs_root).map_err(|error| error.to_string())?;
        let path = runs_root.join(".linkerforge-gpu.lock");
        let mut file = match create_gpu_lock(&path) {
            Ok(file) => file,
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {
                let pid = read_lock_pid(&path)?;
                if process_is_active(pid)? || prediction_or_simulation_active()? {
                    return Err(format!(
                        "Another LinkerForge GPU job may be running (lock PID {pid}). Inspect {} and remove it only after confirming no prediction or simulation process is active.",
                        path.display()
                    ));
                }
                fs::remove_file(&path).map_err(|error| {
                    format!(
                        "Confirmed stale GPU job lock {} could not be removed: {error}",
                        path.display()
                    )
                })?;
                create_gpu_lock(&path).map_err(|error| {
                    format!(
                        "Could not claim GPU job lock {} after stale-lock recovery: {error}",
                        path.display()
                    )
                })?
            }
            Err(error) => {
                return Err(format!(
                    "Could not create GPU job lock {}: {error}",
                    path.display()
                ));
            }
        };
        writeln!(
            file,
            "pid={}\ncreated_unix_seconds={}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(|error| error.to_string())?
                .as_secs()
        )
        .map_err(|error| error.to_string())?;
        Ok(Self {
            path,
            file: Some(file),
        })
    }
}

impl Drop for GpuJobLock {
    fn drop(&mut self) {
        self.file.take();
        let _ = fs::remove_file(&self.path);
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ComputeStatus {
    predictor_available: bool,
    predictor_backend: String,
    gpu: String,
    pymol_available: bool,
    pymol_path: Option<String>,
    note: String,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct PymolLaunchRequest {
    pdb: String,
    sequence: String,
    domain_a_length: usize,
    linker_length: usize,
}

#[derive(Clone)]
enum PymolBackend {
    Native(PathBuf),
    WslOpenSource { executable: PathBuf },
}

impl PymolBackend {
    fn description(&self) -> String {
        match self {
            Self::Native(path) => path.display().to_string(),
            Self::WslOpenSource { executable } => {
                format!("WSL2 open-source PyMOL ({})", executable.display())
            }
        }
    }
}

fn sha256_bytes(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn fold_runs_root(app: &tauri::AppHandle) -> Result<PathBuf, String> {
    if let Ok(configured) = std::env::var("LINKERFORGE_RUNS_DIR") {
        let configured = configured.trim();
        if configured.is_empty() {
            return Err("LINKERFORGE_RUNS_DIR is set but empty.".into());
        }
        let path = PathBuf::from(configured);
        if !path.is_absolute() {
            return Err("LINKERFORGE_RUNS_DIR must be an absolute path.".into());
        }
        return Ok(path);
    }
    app.path()
        .app_local_data_dir()
        .map(|path| path.join("folds"))
        .map_err(|error| error.to_string())
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

fn first_matching_file(output: &Path, required_parts: &[&str], extension: &str) -> Option<PathBuf> {
    let mut files: Vec<PathBuf> = fs::read_dir(output)
        .ok()?
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| {
            let name = path.file_name().unwrap_or_default().to_string_lossy();
            path.extension()
                .is_some_and(|value| value.to_string_lossy().eq_ignore_ascii_case(extension))
                && required_parts.iter().all(|part| name.contains(part))
        })
        .collect();
    files.sort();
    files.into_iter().next()
}

fn colabfold_metrics(
    output: &Path,
    pdb_path: &Path,
    domain_a_length: usize,
    linker_length: usize,
) -> (Option<FoldMetrics>, FoldArtifacts) {
    let score_path = first_matching_file(output, &["scores_rank_001"], "json");
    let pae_plot = first_matching_file(output, &["_pae"], "png");
    let plddt_plot = first_matching_file(output, &["_plddt"], "png");
    let metrics = score_path.as_ref().and_then(|path| {
        let value: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(path).ok()?).ok()?;
        let plddt = value.get("plddt")?.as_array()?;
        let plddt_values: Vec<f64> = plddt.iter().filter_map(serde_json::Value::as_f64).collect();
        let mean_plddt = (!plddt_values.is_empty())
            .then(|| plddt_values.iter().sum::<f64>() / plddt_values.len() as f64);
        let pae = value.get("pae").and_then(serde_json::Value::as_array);
        let domain_b_start = domain_a_length + linker_length;
        let mut cross_domain = Vec::new();
        if let Some(rows) = pae {
            for row_value in rows.iter().take(domain_a_length) {
                if let Some(row) = row_value.as_array() {
                    for cell in row.iter().skip(domain_b_start) {
                        if let Some(number) = cell.as_f64() {
                            cross_domain.push(number);
                        }
                    }
                }
            }
            for row_value in rows.iter().skip(domain_b_start) {
                if let Some(row) = row_value.as_array() {
                    for cell in row.iter().take(domain_a_length) {
                        if let Some(number) = cell.as_f64() {
                            cross_domain.push(number);
                        }
                    }
                }
            }
        }
        let mean_interdomain_pae = (!cross_domain.is_empty())
            .then(|| cross_domain.iter().sum::<f64>() / cross_domain.len() as f64);
        let model_name = pdb_path
            .file_stem()
            .map(|name| name.to_string_lossy().to_string())
            .and_then(|name| name.split("rank_001_").nth(1).map(str::to_string));
        Some(FoldMetrics {
            predictor: "AlphaFold2 model executed by ColabFold".into(),
            model_name,
            mean_plddt,
            ptm: value.get("ptm").and_then(serde_json::Value::as_f64),
            max_pae: value.get("max_pae").and_then(serde_json::Value::as_f64),
            mean_interdomain_pae,
        })
    });
    let artifacts = FoldArtifacts {
        pdb_path: Some(pdb_path.display().to_string()),
        score_json_path: score_path.map(|path| path.display().to_string()),
        pae_plot_path: pae_plot.map(|path| path.display().to_string()),
        plddt_plot_path: plddt_plot.map(|path| path.display().to_string()),
        log_path: None,
    };
    (metrics, artifacts)
}

fn executable_from_path(name: &str) -> Option<PathBuf> {
    let output = hidden_command("where.exe")
        .arg(name)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    String::from_utf8_lossy(&output.stdout)
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(PathBuf::from)
        .find(|path| path.is_file())
}

fn windows_path_to_wsl(path: &Path) -> Result<String, String> {
    let output = hidden_command("wsl.exe")
        .arg("--exec")
        .arg("wslpath")
        .arg("-a")
        .arg(path)
        .output()
        .map_err(|error| error.to_string())?;
    if !output.status.success() {
        return Err(format!("WSL2 could not translate {}", path.display()));
    }
    let translated = String::from_utf8_lossy(&output.stdout).trim().to_string();
    if translated.is_empty() {
        return Err(format!(
            "WSL2 returned an empty path for {}",
            path.display()
        ));
    }
    Ok(translated)
}

fn find_pymol() -> Option<PymolBackend> {
    if let Some(project) = std::env::var_os("LINKERFORGE_PYMOL_WSL_PROJECT").map(PathBuf::from) {
        let executable = project
            .join(".pixi")
            .join("envs")
            .join("default")
            .join("bin")
            .join("pymol");
        if executable.is_file() {
            return Some(PymolBackend::WslOpenSource { executable });
        }
    }
    if let Some(path) = std::env::var_os("LINKERFORGE_PYMOL")
        .map(PathBuf::from)
        .filter(|path| path.is_file())
    {
        return Some(PymolBackend::Native(path));
    }
    for name in ["PyMOLWin.exe", "pymol.exe", "pymol"] {
        if let Some(path) = executable_from_path(name) {
            return Some(PymolBackend::Native(path));
        }
    }
    let mut candidates = Vec::new();
    if let Some(base) = std::env::var_os("LOCALAPPDATA") {
        candidates.push(
            PathBuf::from(base)
                .join("Schrodinger")
                .join("PyMOL2")
                .join("PyMOLWin.exe"),
        );
    }
    for variable in ["ProgramFiles", "ProgramFiles(x86)"] {
        if let Some(base) = std::env::var_os(variable) {
            let base = PathBuf::from(base);
            candidates.push(base.join("Schrodinger").join("PyMOL2").join("PyMOLWin.exe"));
            candidates.push(base.join("PyMOL").join("PyMOLWin.exe"));
        }
    }
    candidates
        .into_iter()
        .find(|path| path.is_file())
        .map(PymolBackend::Native)
}

fn pymol_script(domain_a_length: usize, linker_length: usize, sequence_length: usize) -> String {
    let linker_start = domain_a_length + 1;
    let linker_end = domain_a_length + linker_length;
    let domain_b_start = linker_end + 1;
    let linker_commands = if linker_length == 0 {
        String::new()
    } else {
        format!("select linker, resi {linker_start}-{linker_end}\ncolor magenta, linker\n")
    };
    format!(
        "hide everything, all\nshow cartoon, all\nselect domain_a, resi 1-{domain_a_length}\n{linker_commands}select domain_b, resi {domain_b_start}-{sequence_length}\ncolor green, domain_a\ncolor red, domain_b\nbg_color white\nset cartoon_transparency, 0\nzoom all\n"
    )
}

fn probe_compute() -> ComputeStatus {
    let native_available = hidden_command("where.exe")
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
        && hidden_command("wsl.exe")
            .args(["--exec", "which", "colabfold_batch"])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|status| status.success());
    let gpu = hidden_command("nvidia-smi")
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
    let pymol = find_pymol();
    ComputeStatus {
        predictor_available: native_available || wsl_available,
        predictor_backend: if native_available { "Windows PATH" } else if wsl_available { "WSL2 PATH" } else { "unavailable" }.into(),
        gpu,
        pymol_available: pymol.is_some(),
        pymol_path: pymol.map(|backend| backend.description()),
        note: "LocalColabFold must be installed separately on Windows PATH or WSL2 PATH. WSL2 probing is opt-in via LINKERFORGE_ENABLE_WSL=1. PyMOL is optional and used only for inspection; an explicit LINKERFORGE_PYMOL_WSL_PROJECT can select an isolated open-source pixi environment. GPU memory may be insufficient for a full fusion; no model is bundled.".into(),
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
async fn open_pdb_in_pymol(
    app: tauri::AppHandle,
    request: PymolLaunchRequest,
) -> Result<String, String> {
    validate_sequence(&request.sequence)?;
    if request.pdb.len() as u64 > MAX_PDB_BYTES {
        return Err("The local PDB limit is 15 MB.".into());
    }
    if request.domain_a_length == 0
        || request.domain_a_length + request.linker_length >= request.sequence.len()
    {
        return Err("Domain and linker boundaries do not describe the submitted fusion.".into());
    }
    let observed = read_ca_sequence(&request.pdb)?;
    if observed != request.sequence {
        return Err(
            "The PDB C-alpha sequence does not exactly match the selected full fusion.".into(),
        );
    }
    let pymol = find_pymol().ok_or(
        "PyMOL was not found. Set LINKERFORGE_PYMOL to a licensed Windows executable or LINKERFORGE_PYMOL_WSL_PROJECT to an open-source pixi project.",
    )?;
    let base = app
        .path()
        .app_local_data_dir()
        .map_err(|error| error.to_string())?;
    let digest = sha256_bytes(request.sequence.as_bytes());
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|error| error.to_string())?
        .as_millis();
    let view_dir = base
        .join("pymol-views")
        .join(format!("{timestamp}-{}", &digest[..12]));
    let script = pymol_script(
        request.domain_a_length,
        request.linker_length,
        request.sequence.len(),
    );
    let pdb_path = view_dir.join("fusion.pdb");
    let pml_path = view_dir.join("colour-fusion.pml");
    tauri::async_runtime::spawn_blocking(move || {
        fs::create_dir_all(&view_dir).map_err(|error| error.to_string())?;
        fs::write(&pdb_path, request.pdb).map_err(|error| error.to_string())?;
        fs::write(&pml_path, script).map_err(|error| error.to_string())?;
        match pymol {
            PymolBackend::Native(executable) => {
                hidden_command(executable)
                    .arg(&pdb_path)
                    .arg(&pml_path)
                    .spawn()
                    .map_err(|error| format!("Could not start PyMOL: {error}"))?;
            }
            PymolBackend::WslOpenSource { executable } => {
                let executable = windows_path_to_wsl(&executable)?;
                let pdb = windows_path_to_wsl(&pdb_path)?;
                let script = windows_path_to_wsl(&pml_path)?;
                hidden_command("wsl.exe")
                    .arg("--exec")
                    .arg(executable)
                    .arg(pdb)
                    .arg(script)
                    .spawn()
                    .map_err(|error| {
                        format!("Could not start open-source PyMOL through WSL2: {error}")
                    })?;
            }
        }
        Ok(format!(
            "Opened the sequence-matched fusion in PyMOL from {}",
            pdb_path.display()
        ))
    })
    .await
    .map_err(|error| error.to_string())?
}

#[tauri::command]
async fn start_fold(
    app: tauri::AppHandle,
    state: tauri::State<'_, FoldState>,
    sequence: String,
    domain_a_length: usize,
    linker_length: usize,
) -> Result<FoldStatus, String> {
    validate_sequence(&sequence)?;
    if domain_a_length == 0 || domain_a_length + linker_length >= sequence.len() {
        return Err("Domain and linker boundaries do not describe the submitted fusion.".into());
    }
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
    let runs_root = fold_runs_root(&app)?;
    let gpu_lock = GpuJobLock::acquire(&runs_root)?;
    let run_dir = runs_root.join(format!("{seconds}-{}", &digest[..12]));
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
        metrics: None,
        artifacts: None,
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
                    windows_path_to_wsl(&fasta)?,
                    windows_path_to_wsl(&output)?,
                )
            } else {
                (
                    "colabfold_batch",
                    Vec::new(),
                    fasta.display().to_string(),
                    output.display().to_string(),
                )
            };
            let mut command = hidden_async_command(program);
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
            Ok((pdb_path, pdb)) => {
                let (metrics, mut artifacts) =
                    colabfold_metrics(&output, &pdb_path, domain_a_length, linker_length);
                artifacts.log_path = Some(log_path.display().to_string());
                FoldStatus {
                    phase: "complete".into(),
                    message: format!("Matching full-fusion PDB: {}", pdb_path.display()),
                    run_dir: Some(run_dir.display().to_string()),
                    sequence_sha256: Some(digest.clone()),
                    pdb: Some(pdb),
                    metrics,
                    artifacts: Some(artifacts),
                }
            }
            Err(message) => FoldStatus {
                phase: "failed".into(),
                message,
                run_dir: Some(run_dir.display().to_string()),
                sequence_sha256: Some(digest.clone()),
                pdb: None,
                metrics: None,
                artifacts: Some(FoldArtifacts {
                    pdb_path: None,
                    score_json_path: None,
                    pae_plot_path: None,
                    plddt_plot_path: None,
                    log_path: Some(log_path.display().to_string()),
                }),
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
            "metrics": final_status.metrics.clone(),
            "artifacts": final_status.artifacts.clone(),
            "boundaries": {"domain_a_length": domain_a_length, "linker_length": linker_length, "domain_b_length": sequence.len() - domain_a_length - linker_length},
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
        drop(gpu_lock);
    });
    Ok(running)
}

pub fn run() {
    tauri::Builder::default()
        .manage(FoldState::default())
        .invoke_handler(tauri::generate_handler![
            compute_status,
            fold_status,
            start_fold,
            open_pdb_in_pymol
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
    fn configured_runs_root_must_be_absolute() {
        assert!(Path::new(r"M:\LinkerForge\folds").is_absolute());
        assert!(!Path::new("relative/folds").is_absolute());
    }

    #[test]
    fn gpu_job_lock_is_cross_process_and_recoverable_after_drop() {
        let directory = std::env::temp_dir().join(format!(
            "linkerforge-lock-test-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("system clock")
                .as_nanos()
        ));
        let first = GpuJobLock::acquire(&directory).expect("first lock");
        assert!(GpuJobLock::acquire(&directory).is_err());
        drop(first);
        let second = GpuJobLock::acquire(&directory).expect("lock after drop");
        drop(second);
        fs::remove_dir(&directory).expect("remove lock test directory");
    }

    #[test]
    fn gpu_job_lock_recovers_only_after_dead_pid_and_process_scan() {
        let directory = std::env::temp_dir().join(format!(
            "linkerforge-stale-lock-test-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("system clock")
                .as_nanos()
        ));
        fs::create_dir(&directory).expect("create stale lock test directory");
        let path = directory.join(".linkerforge-gpu.lock");
        fs::write(&path, "pid=2147483647\ncreated_unix_seconds=1\n").expect("write stale lock");
        let recovered = GpuJobLock::acquire(&directory).expect("recover stale lock");
        assert!(path.exists());
        drop(recovered);
        fs::remove_dir(&directory).expect("remove stale lock test directory");
    }

    #[test]
    fn gpu_job_lock_leaves_malformed_lock_in_place() {
        let directory = std::env::temp_dir().join(format!(
            "linkerforge-malformed-lock-test-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("system clock")
                .as_nanos()
        ));
        fs::create_dir(&directory).expect("create malformed lock test directory");
        let path = directory.join(".linkerforge-gpu.lock");
        fs::write(&path, "not-a-valid-lock\n").expect("write malformed lock");
        assert!(GpuJobLock::acquire(&directory).is_err());
        assert!(path.exists());
        fs::remove_file(&path).expect("remove malformed lock");
        fs::remove_dir(&directory).expect("remove malformed lock test directory");
    }

    #[test]
    fn pdb_gate_requires_full_trace() {
        let line = "ATOM      1  CA  GLY A   1      0.000   0.000   0.000  1.00 50.00           C";
        assert_eq!(read_ca_sequence(line).unwrap(), "G");
    }

    #[test]
    fn pymol_script_marks_exact_boundaries() {
        let script = pymol_script(238, 15, 489);
        assert!(script.contains("select domain_a, resi 1-238"));
        assert!(script.contains("select linker, resi 239-253"));
        assert!(script.contains("select domain_b, resi 254-489"));
    }

    #[test]
    fn pymol_script_supports_direct_fusion_control() {
        let script = pymol_script(20, 0, 62);
        assert!(script.contains("select domain_a, resi 1-20"));
        assert!(!script.contains("select linker"));
        assert!(script.contains("select domain_b, resi 21-62"));
    }
}
