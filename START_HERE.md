# LinkerForge: start here

`README.md` is the canonical scientific specification. This page is a map of
the executable parts of the repository and their evidence boundaries.

| Workstream | Start with | Current boundary |
| --- | --- | --- |
| Sequence-based DFL detector | `RESEARCH_HANDOFF.md`, `CLASSIFIER_DATA_GATE.md`, `src/linkerforge/detector/`, `scripts/check_detector_autoresearch.py` | The historical checkpoint is exploratory. The published B1 V2 split reproduces its preregistered sequence/family preflight, but that is not authorization for a new blind comparison: CAID3 is spent retrospective and the legacy/current prospective constraints remain binding. Do not mine or claim precision from it. |
| Fusion design components | `ARCHITECTURE_HANDOFF.md`, `DESIGN_PIPELINE.md`, `src/linkerforge/{geometry,properties,design,attachment,pipeline}/` | Geometry and software scaffolding have CPU tests. There is no biologically trained end-to-end fusion ranker or validated fusion-success predictor. |
| Web demonstration | `frontend/`, `docs/FULL_FUSION_LOCAL_WORKFLOW.md` | Enter two sequences, assemble labelled template candidates, export FASTA and inspect a sequence-matched predicted PDB. The browser does **not** fold proteins. |
| Windows desktop | `frontend/src-tauri/`, `docs/DESKTOP_LOCAL_FOLDING.md`, `scripts/start_linkerforge_desktop.ps1` | LocalColabFold, monomer/PTM parameters and an isolated open-source PyMOL environment are installed on M:. A documented 142-residue 1UBQ–linker–2GB1 run completed and rendered. Visual Studio Build Tools/Windows SDK are installed, native tests pass, and a release executable builds. The approximately 4 GB GPU remains unsuitable for assuming that a roughly 500-residue fusion will fit. |
| Professional/publication audit | `docs/PROFESSIONAL_PUBLICATION_READINESS_AUDIT_2026-10-20.md` | Current implementation, 28-construct published-detector panel, completed structural checklist, model comparison, publication boundary and unresolved scientific/product gates. |
| Data, evaluation and claims | `datasets/`, `schemas/`, `configs/`, `reports/`, `literature/DFL_UNTOUCHED_COMPARISON_CONTRACT_2026-10-17.md` | The GitHub release contains the historical 60/28/9 DFL split and a 27,766-record mined/natural source inventory. The latter fails closed for fusion-ranker training because it has no ranker-specific family splits or distinct compatibility negatives. Read the applicable provenance and benchmark contract before training or publishing a result. |

## Layout

`frontend/` contains the React/Vite UI and its optional `src-tauri/` desktop
host. `src/linkerforge/` is the Python package. `scripts/` contains thin entry
points; `tests/` contains CPU/synthetic checks. `schemas/` defines records,
`configs/` pins experiments, `data/curated/` holds small cited source records,
and `reports/` records results and failed gates. `literature/` holds research
contracts and surveys. Large `data/raw/`, `data/processed/`, `artifacts/`, and
desktop `target/` outputs are ignored; a clean clone does not contain them.

Historical proposals, papers and handoffs remain in the root or `paper/` for
traceability. They are not interchangeable with current validation evidence.

## Fast checks

From `frontend/`, run `npm ci`, `npm run lint`, `npm run test`, and
`npm run build`. With Rust and the Windows Tauri prerequisites installed,
`npm run desktop:build` builds the desktop executable. `cargo test` in
`frontend/src-tauri/` checks its native validation gates. Desktop installation
does not install ColabFold or model weights. On the audited workstation,
`.\scripts\start_linkerforge_desktop.ps1` performs the local WSL/GPU preflight
and starts the development shell; add `-Mode built` for the compiled release.

For Python work, use the component-specific commands in its handoff and report.
Do not infer that passing software tests unlocks the detector's B1 gate or a
biological performance claim.
