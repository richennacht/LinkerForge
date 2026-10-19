# LinkerForge desktop: local full-fusion folding

The Tauri desktop shell reuses the React interface. Enter both domain sequences,
generate a template, select it, and choose **Fold full fusion on this GPU**.
The app sends `Domain A + linker + Domain B` to a local `colabfold_batch`
process, waits for a complete single-chain PDB whose C-alpha sequence exactly
matches the submitted fusion, and colours the three segments in its 3D viewer.
Unranked templates and a fold are *predictions*, not biological validation.

## What is and is not bundled

- The desktop application and local job controller are bundled.
- ColabFold, model weights, MSA databases/service access, CUDA, and GPU drivers
  are **not** bundled. `colabfold_batch` must be separately installed and on
  Windows `PATH` or WSL2 `PATH`; verify its version and GPU setup using its own
  instructions. The official LocalColabFold Windows route uses WSL2. To enable
  the optional WSL2 backend, set `LINKERFORGE_ENABLE_WSL=1` before launching
  the app. The app does not probe an uninstalled WSL stub by default.
- The app does not fall back to CPU or cloud without saying so. If ColabFold is
  absent, the fold button is disabled; FASTA export and matching-PDB import
  remain usable.
- The full-fusion prediction uses one local GPU job at a time, fast kernels,
  JAX non-preallocation and the platform allocator, and a 120-minute limit.
  The WSL2 route also applies LocalColabFold's documented unified-memory
  settings, allowing some GPU-memory oversubscription into system RAM.
  These settings reduce peak memory but do not guarantee that a 6 GB GPU can
  fold a roughly 500-residue fusion. Out-of-memory exits remain visible as a
  failed run, with the log retained.
- Each run is isolated under the app's local data directory in `folds/` and
  writes `fusion.fasta`, `colabfold.log`, `run_manifest.json`, and ColabFold
  outputs. The UI shows the run path. No predicted structure is called an
  experimental result.

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
separately. WSL2 is not installed on the current workstation, so that route
could not be exercised here. The browser app remains available through `npm run dev`.

## Memory strategy

Moving weights between CPU and GPU or freeing layer memory can help some
inference workloads, but AlphaFold/ColabFold also keeps sequence-pair and MSA
activations whose size grows with the *whole* sequence. Folding disconnected
chunks independently would not predict the interdomain orientation or the
linker's complete-chain geometry. Therefore this desktop path uses a complete
fusion, memory-saving inference settings, and a fail-closed OOM outcome. It
does not claim that streaming makes arbitrary fusion length fit in 6 GB.

See [ColabFold's GPU guidance](https://github.com/sokrypton/ColabFold) and
[Tauri's existing-frontend setup](https://v2.tauri.app/start/create-project/).
