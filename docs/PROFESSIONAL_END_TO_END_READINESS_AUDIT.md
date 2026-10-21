# LinkerForge professional end-to-end readiness audit

> Historical installation audit. For the updated 20 October implementation,
> evidence checklist, 28-construct detector comparison and publication decision,
> use `PROFESSIONAL_PUBLICATION_READINESS_AUDIT_2026-10-20.md`. Resource and build
> values below describe the earlier checkpoint and must not be treated as current.

Audit date: 19–20 October 2026 (Asia/Calcutta)
Audited checkout: `M:\ChatGPT\LinkerForge 2`, commit base `8dfefb0` plus the
uncommitted Tauri/PyMOL integration in this checkout.

## Decision

**Not yet ready for professional scientific use.** The implemented application
is a useful, tested software prototype. It can validate two sequences, construct
clearly labelled template-based fusion candidates, export a full-fusion FASTA,
accept an exactly sequence-matched PDB, display it with 3Dmol.js, and launch
PyMOL with domain/linker colouring. The code-only generator, attachment-site,
geometry, checkpoint and orchestration interfaces pass their synthetic software
tests.

It cannot yet make a validated biological recommendation. There is no reviewed
fusion-training dataset, scientifically trained checkpoint, calibrated ranker,
or completed structural funnel. A synthetic smoke checkpoint must never be used
as a research model. LocalColabFold's software environment and monomer/PTM
parameters are installed, and a short GPU smoke prediction completed. The
489-residue demonstration has not been folded and the smoke output is not
evidence that the full fusion will fit in 4 GB VRAM.

## Workstation audit

| Resource | Observed | Readiness consequence |
| --- | --- | --- |
| OS | Windows 11, build 26200 | WSL2 is enabled and operational. |
| GPU | NVIDIA GeForce RTX 3050 A Laptop GPU, 4094 MiB, driver 616.56 | CUDA-capable, but 4 GB VRAM is a serious risk for a roughly 500-residue complete fusion. A successful fold is not assumed. |
| RAM | 23.7 GB | Adequate for the UI and software tests; not enough for a complete local ColabFold sequence database. |
| C: free | 3.5 GB | Critically low. Do not install a default WSL distribution or model cache here. |
| M: free | 4.1 GB after runtime, parameters, smoke output and a lean 1.29 GB Tauri rebuild | Too little for a full-fusion attempt with a safe failure/output margin. The previous 4.5 GiB debug tree was replaced using no incremental compilation and stripped test/dev debug symbols. |
| WSL | WSL 3.0.3, kernel 6.18.54.1, Ubuntu 24.04 on WSL2 | Distribution is stored at `M:\WSL\LinkerForge-Ubuntu`; default user is the dedicated unprivileged `linkerforge` account. |
| PyMOL | Incentive PyMOL 3.1.8 starts and renders via OpenGL 4.6 | Visualization works, but no licence file is installed and rendered images carry an evaluation watermark. Acquire an appropriate licence or approve an open-source PyMOL build for professional output. |
| Python model environment | Repository `.venv`, Python 3.14.5, PyTorch 2.14.1 CPU | Synthetic generator and architecture verification run successfully. This is not the approved college-cloud GPU image. |

Official ColabFold supports Windows through WSL2 using LocalColabFold. Its local
database route requires about 940 GB of storage and around 128 GB RAM for batch
searches, so this workstation must use the serial public MSA server or an
institutional service. The MSA-server choice needs a privacy/data-governance
notice because sequences leave the workstation.

Primary references:

- <https://github.com/sokrypton/ColabFold/blob/main/README.md>
- <https://github.com/YoshitakaMo/localcolabfold/blob/main/README.md>
- <https://learn.microsoft.com/windows/wsl/install>
- <https://docs.nvidia.com/cuda/wsl-user-guide/index.html>

## Verified in this audit

| Area | Evidence | Result |
| --- | --- | --- |
| Desktop UI | Tauri window starts from `frontend/`; local page responds HTTP 200 | Pass |
| Frontend/native checks | Lint, 10 frontend tests, production build, 3 Rust tests | Pass |
| Design architecture | `verify_design_architecture.py --scope all` | 165/165 synthetic/software tests pass |
| Generator wiring | Three bounded CPU optimizer steps, checkpoint save/reload, one generated candidate | Pass as synthetic software only |
| Template design | Example pipeline returns `(GGGGS)3` for the synthetic tandem request | Pass as unranked sequence proposal |
| PDB safety gate | Complete single-chain C-alpha sequence must exactly match the selected fusion | Pass |
| PyMOL | Illustrative 489-residue PDB loads, colours, and renders | Pass for visualization only |
| Local ColabFold | LocalColabFold 1.6.1 installer at commit `e2e8c4a834a3d50ebaede3bb42786707fdc8f6ec`; runtime reports ColabFold 1.6.3 at `297da92fdff1348b8c99c2a946d7806a4eb4c25e`; JAX 0.10.1 reports `CudaDevice(id=0)` | One 62-residue, single-sequence, one-model/one-recycle GPU smoke run passed in 22.1 s; this is operational evidence only |
| ColabFold → PyMOL | Ranked smoke PDB rendered through PyMOL's Python API | Pass for visualization handoff; evaluation watermark remains |
| Real learned linker model | No approved data or trained scientific checkpoint | Blocked |
| Simulation | No OpenMM/GROMACS execution is connected to the desktop | Not implemented |

ColabFold predicts a structure; PyMOL inspects or renders it. Neither is a
molecular-dynamics simulation. If simulation is required, it must be a separate
OpenMM/GROMACS stage with preparation, force-field, solvent, equilibration,
replicates, convergence and trajectory-quality records.

## Professional user journey and required evidence

1. **Create a project.** Record researcher/project IDs, purpose, consent/privacy
   classification, organism/source IDs, exact domain sequences, boundary source,
   protected sites and input hashes.
2. **Run input preflight.** Show canonical-residue checks, lengths, duplicate or
   family exposure checks, estimated compute risk, available GPU/RAM/disk, tool
   versions and whether a sequence will be sent to an external MSA service.
3. **Choose a design method.** Distinguish fixed control, cited template,
   retrieval result, trained-model prediction and manually supplied linker.
   Display model/checkpoint/config/data identities and evidence quality. Disable
   trained-model claims when no approved checkpoint is loaded.
4. **Review candidates.** Show full sequence, junctions, linker properties,
   provenance, all rejection reasons, uncertainty/abstention, controls and why a
   candidate advanced. Never label a sequence proposal as successful.
5. **Predict the complete fusion.** Write an immutable FASTA and run manifest;
   record ColabFold version, backend, model type, MSA source/hash, templates,
   seeds, recycles, relaxation, environment, timings, log and every output hash.
6. **Assess structure.** Report pLDDT by domain/linker, PAE/interdomain
   uncertainty, clashes, chain completeness, domain preservation, attachment
   geometry, rejected/accepted status and matched baseline scores. Preserve all
   failed runs.
7. **Visualize.** Provide inline 3D viewing, confidence colouring and exact
   domain/linker selections; open the same hash-checked PDB in PyMOL. Record the
   visualization script and tool/licence version.
8. **Optionally simulate finalists.** Treat minimization/MD as a new, explicit
   job with its own validated protocol and outputs, not as a PyMOL action.
9. **Export a research bundle.** Include machine-readable JSON, FASTA, PDB/mmCIF,
   PAE/confidence files, images, logs, manifests, citations, comparison table,
   software versions, commit SHA/config hashes and a plain-language limitations
   statement.
10. **Reopen and reproduce.** A researcher must be able to reload a project,
    inspect run history, verify hashes and rerun without overwriting the original.

## Product gaps before professional use

### P0 — blocks a trustworthy end-to-end run

- Preserve the installed WSL2/Ubuntu and pinned LocalColabFold commit; add a
  reproducible lockfile/archive policy before treating it as a deployment.
- Free substantially more disk space; define cache/output retention and cleanup.
- Free enough space to retain a safe full-fusion result/failure margin, then
  exercise the short smoke job through the Tauri controller before any full
  fusion.
- Add ColabFold version/model/weights/MSA identities to the run manifest.
- Add cancellation, stale-job recovery and first-run dependency diagnostics.
- Add explicit public-MSA data disclosure/consent or use an approved private
  institutional MSA service.
- Provide reviewed training data and a frozen checkpoint, or keep the product in
  template-only mode with no learned recommendation language.
- Add formal JSON Schema for all scientific run/result records before emitting
  them as scientific evidence.
- Resolve the PyMOL licence for professional figures.

### P1 — required researcher results

- Project save/open and immutable run history.
- Candidate comparison against direct fusion and `(GGGGS)3` under identical
  prediction/scoring rules.
- pLDDT, PAE, clash, preservation and linker/attachment metrics in the UI.
- Downloadable provenance bundle and human-readable report.
- Clear failure diagnosis for MSA, compilation, CUDA OOM, timeout, invalid PDB,
  licensing and network errors.
- Resource estimator and a cloud handoff when a job is unsuitable for local GPU.

### P2 — later scientific workflow

- Rosetta linker conformational sampling and hard-rejection funnel.
- OpenMM minimization and replicated MD for only reviewed finalists.
- Optional GROMACS for 1–3 finalists.
- Validated learned generator/site scorer/ranker, calibration and prospective
  family-disjoint evaluation under the repository's benchmark firewall.

## Installation gate and safe order

Do not start a full fusion yet. The 3,722,752,000-byte monomer/PTM parameter set
is now installed, but M: has only about 4.1 GB free after a lean Tauri rebuild.
More space is required before a long run is attempted.

1. Preserve the current repository changes and free disk space. Target at least
   30–50 GB free on M: for WSL, the LocalColabFold environment, weights, caches,
   logs and several runs. Keep C: above a safe Windows operating margin.
2. Enable WSL without a distribution, restart Windows, then install Ubuntu with
   `--location` pointing to a dedicated directory on M:.
3. Update WSL and verify the Windows NVIDIA driver is exposed inside WSL. Do not
   install a Linux display driver in WSL.
4. Install LocalColabFold with its current `pixi` workflow inside the Linux file
   system, not a case-insensitive Windows-mounted project directory.
5. Verify `colabfold_batch --help`, versions and a short official sample.
6. Enable `LINKERFORGE_ENABLE_WSL=1`, restart LinkerForge, and confirm its
   compute-status card detects `WSL2 PATH`.
7. Run a short sequence before attempting the 489-residue demonstration. Treat
   a 4 GB OOM as an expected capacity finding and use college-cloud GPU compute.

### Current installation checkpoint

- Completed: WSL2 and Ubuntu 24.04 installed at
  `M:\WSL\LinkerForge-Ubuntu`; CPU virtualization, the Windows hypervisor and
  NVIDIA WSL passthrough are operational.
- Completed: LocalColabFold 1.6.1 installed in
  `/home/linkerforge/localcolabfold` from commit
  `e2e8c4a834a3d50ebaede3bb42786707fdc8f6ec` using Pixi 0.82.0.
- Completed: package cache removed, VHD trimmed/compacted and sparse mode
  enabled. The installed environment remains approximately 8.9 GB.
- Completed: repository wrapper `scripts/colabfold_batch_wsl.sh` installed as
  `/usr/local/bin/colabfold_batch`; default-user, command, help and JAX CUDA
  probes pass. Windows user setting `LINKERFORGE_ENABLE_WSL=1` is enabled.
- Intentionally absent: local sequence databases. MSA searches use the public
  MMseqs2 service, so sequences leave the workstation and require user consent.
- Completed: official AlphaFold2-PTM parameters installed under
  `/home/linkerforge/.cache/colabfold/params`; the upstream completion marker is
  present. No multimer parameter set was downloaded.
- Completed: the bundled 62-residue `1BJP_1.fasta` ran with `single_sequence`,
  one model, one recycle, no relaxation and seed 0. ColabFold reported GPU use,
  22.1 seconds, pLDDT 57.2 and pTM 0.375. These scores describe a smoke test,
  not a validated research result.
- Completed: the ranked smoke PDB rendered in PyMOL to
  `tmp/colabfold-smoke-20261019/pymol_confidence.png`; the evaluation watermark
  confirms that PyMOL licensing remains unresolved.
- Completed: frontend lint, 10 tests and production build pass; 3 Rust tests
  pass. A lean Tauri development build is running with the WSL opt-in enabled;
  the native target occupies approximately 1.29 GB.
- Pending: free more space, run the same short job through the Tauri controller,
  add full provenance fields, and only then assess whether the complete fusion
  should run locally or on the college cloud. No molecular-dynamics simulation
  result exists.

## Current evidence artifacts

- `tmp/generator-smoke-restart-20261019/smoke.json`
- `tmp/generator-smoke-restart-20261019/development-software-only.json`
- `tmp/design-software-example-restart-20261019.json`
- `tmp/design-checks-restart-20261019.json`
- `tmp/pymol-linkerforge-restart-20261019.png`
- `tmp/colabfold-smoke-20261019/`

These files are software/audit evidence only. The first PyMOL image uses an
explicitly illustrative geometric PDB. The ColabFold smoke directory contains a
real short-sequence prediction and its PyMOL render, but it is neither fusion
validation nor a molecular-dynamics simulation.
