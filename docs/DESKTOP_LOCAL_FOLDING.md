# LinkerForge desktop: local full-fusion folding

The Tauri desktop shell reuses the React interface. Enter both domain sequences,
generate a template, select it, and choose **Fold full fusion on this GPU**.
The app sends `Domain A + linker + Domain B` to a local `colabfold_batch`
process, waits for a complete single-chain PDB whose C-alpha sequence exactly
matches the submitted fusion, and colours the three segments in its 3D viewer.
When PyMOL is installed, **Open coloured fusion in PyMOL** writes an isolated
copy of that same sequence-checked PDB plus a colouring script and launches
PyMOL with Domain A, linker, and Domain B highlighted separately.
Unranked templates and a fold are *predictions*, not biological validation.

## What is and is not bundled

- The desktop application and local job controller are bundled.
- ColabFold, model weights, MSA databases/service access, CUDA, and GPU drivers
  are **not** bundled. `colabfold_batch` must be separately installed and on
  Windows `PATH` or WSL2 `PATH`; verify its version and GPU setup using its own
  instructions. The official LocalColabFold Windows route uses WSL2. To enable
  the optional WSL2 backend, set `LINKERFORGE_ENABLE_WSL=1` before launching
  the app. The app does not probe an uninstalled WSL stub by default.
- The audited workstation uses the smallest supported data configuration: no
  local sequence databases. It sends sequence searches to ColabFold's public
  MMseqs2 MSA service. Researchers must be told that submitted sequences leave
  the workstation before they start a run.
- The app does not fall back to CPU or cloud without saying so. If ColabFold is
  absent, the fold button is disabled; FASTA export and matching-PDB import
  remain usable.
- PyMOL is optional and is used only for interactive inspection, measurements,
  and figures; it is not a folding or molecular-dynamics engine. The desktop
  app detects `PyMOLWin.exe` on `PATH` and common Windows install locations.
  Set `LINKERFORGE_PYMOL` to the full executable path for another installation.
  Alternatively, set `LINKERFORGE_PYMOL_WSL_PROJECT` to a Windows path holding
  a pixi project whose `.pixi/envs/default/bin/pymol` is an open-source Linux
  PyMOL executable. This explicit setting takes priority over an automatically
  detected Windows installation and launches through WSL2/WSLg.
  A PDB is rechecked against the complete selected fusion before PyMOL starts.
- The full-fusion prediction uses one local GPU job at a time, fast kernels,
  JAX non-preallocation and the platform allocator, and a 120-minute limit.
  The WSL2 route also applies LocalColabFold's documented unified-memory
  settings, allowing some GPU-memory oversubscription into system RAM.
  These settings reduce peak memory but do not guarantee that this workstation's
  approximately 4 GB RTX 3050 Laptop GPU can
  fold a roughly 500-residue fusion. Out-of-memory exits remain visible as a
  failed run, with the log retained.
- Each run is isolated under the app's local data directory in `folds/` and
  writes `fusion.fasta`, `colabfold.log`, `run_manifest.json`, and ColabFold
  outputs. The UI shows the run path. No predicted structure is called an
  experimental result.
- Set `LINKERFORGE_RUNS_DIR` to an absolute directory to keep run inputs,
  outputs, logs and manifests on a dedicated data drive. The audited workstation
  uses `M:\ChatGPT\LinkerForge 2\runs\folds`; an empty or relative override is
  rejected rather than silently falling back to C:.

## Build and run on Windows

Install Node.js, Rust's MSVC toolchain, the Microsoft C++ Build Tools and
WebView2 as described by [Tauri's Windows prerequisites](https://v2.tauri.app/start/prerequisites/).
From `frontend/`:

```powershell
npm ci
npm run desktop
```

To build a local executable without an installer:

```powershell
npm run desktop:build
```

The executable is expected at `frontend/src-tauri/target/release/linkerforge-desktop.exe`.
Run it from a Windows environment in which `nvidia-smi` works and either
`colabfold_batch` is on Windows `PATH`, or WSL2 has `colabfold_batch` on its
default `PATH`. LocalColabFold's platform and CUDA instructions must be followed
separately. On the audited workstation, Ubuntu 24.04 and LocalColabFold 1.6.3
are installed on M:, `colabfold_batch --help` passes, and JAX 0.10.1 detects the
RTX 3050 as `CudaDevice(id=0)`. The repository's
`scripts/colabfold_batch_wsl.sh` wrapper is installed at
`/usr/local/bin/colabfold_batch`, and the Windows user setting
`LINKERFORGE_ENABLE_WSL=1` is enabled. The monomer/PTM parameter set is
installed, and a 62-residue, one-model/one-recycle, single-sequence GPU smoke
prediction completed in 22.1 seconds. A later 142-residue public
1UBQ–`GGGGSGGGGS`–2GB1 acceptance run completed in about 697 seconds, and its
ranked PDB rendered in PyMOL. This proves the local command/GPU/visualization
handoff for those inputs, not that the roughly 489-residue demonstration will
fit or be biologically valid. A separate conda-forge `pymol-open-source` 3.1.0
pixi environment on M: rendered the accepted structure headlessly without the
unlicensed Incentive executable. Its package metadata uses the
`LicenseRef-PyMOL` identifier, so redistribution/licence review remains a
release-owner responsibility. Visual Studio 2022 Build Tools 17.14 and Windows
SDK 10.0.26100 are now available and the native Rust test suite passes.

For the audited M:-only workstation, launch a development session with:

```powershell
$env:LINKERFORGE_ENABLE_WSL = "1"
$env:LINKERFORGE_RUNS_DIR = "M:\ChatGPT\LinkerForge 2\runs\folds"
$env:LINKERFORGE_PYMOL_WSL_PROJECT = "M:\ChatGPT\LinkerForgeTools\pymol-open-source"
$env:CARGO_INCREMENTAL = "0"
$env:CARGO_PROFILE_DEV_DEBUG = "0"
npm run desktop
```

The same preflight and environment setup is available as a one-command entry
point from the repository root:

```powershell
.\scripts\start_linkerforge_desktop.ps1
```

After a release build, launch it with
`.\scripts\start_linkerforge_desktop.ps1 -Mode built`.

## Memory strategy

Moving weights between CPU and GPU or freeing layer memory can help some
inference workloads, but AlphaFold/ColabFold also keeps sequence-pair and MSA
activations whose size grows with the *whole* sequence. Folding disconnected
chunks independently would not predict the interdomain orientation or the
linker's complete-chain geometry. Therefore this desktop path uses a complete
fusion, memory-saving inference settings, and a fail-closed OOM outcome. It
does not claim that streaming makes arbitrary fusion length fit in 4 GB.

See [ColabFold's GPU guidance](https://github.com/sokrypton/ColabFold) and
[Tauri's existing-frontend setup](https://v2.tauri.app/start/create-project/).
