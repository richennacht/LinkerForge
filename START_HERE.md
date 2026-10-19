# LinkerForge: start here

`README.md` is the canonical scientific specification. This page is a map of
the executable parts of the repository and their evidence boundaries.

| Workstream | Start with | Current boundary |
| --- | --- | --- |
| Sequence-based DFL detector | `RESEARCH_HANDOFF.md`, `CLASSIFIER_DATA_GATE.md`, `src/linkerforge/detector/`, `scripts/check_detector_autoresearch.py` | The historical checkpoint is exploratory. The B1 data/split gate is **not passed**; CAID3 is spent retrospective and the old split fails the prospective firewall. Do not mine or claim precision from it. |
| Fusion design components | `ARCHITECTURE_HANDOFF.md`, `DESIGN_PIPELINE.md`, `src/linkerforge/{geometry,properties,design,attachment,pipeline}/` | Geometry and software scaffolding have CPU tests. There is no biologically trained end-to-end fusion ranker or validated fusion-success predictor. |
| Web demonstration | `frontend/`, `docs/FULL_FUSION_LOCAL_WORKFLOW.md` | Enter two sequences, assemble labelled template candidates, export FASTA and inspect a sequence-matched predicted PDB. The browser does **not** fold proteins. |
| Windows desktop | `frontend/src-tauri/`, `docs/DESKTOP_LOCAL_FOLDING.md` | Builds and launches locally. Full-fusion folding requires a separate working ColabFold installation; the 6 GB GPU has not completed a verified fold. |
| Data, evaluation and claims | `schemas/`, `configs/`, `reports/`, `literature/DFL_UNTOUCHED_COMPARISON_CONTRACT_2026-10-17.md` | Read the applicable provenance and benchmark contract before training or publishing a result. |

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
does not install ColabFold or model weights.

For Python work, use the component-specific commands in its handoff and report.
Do not infer that passing software tests unlocks the detector's B1 gate or a
biological performance claim.
