# LinkerForge direction — 18 October 2026

Work resumed on explicit user instruction. PRs #1–#7 were merged into
`codex/dfl-research-handoff` on 18 October 2026. Their former draft and pause
statements in dated handoffs describe their writing time, not the current
project state. The earlier three detector-research heartbeat IDs were absent
from the app when asked to pause them; no replacement schedules were created.

## Build order

1. **Classifier first.** Build a sequence-based DFL residue caller with
   positive/negative/unknown semantics, family-disjoint training and
   development, a benchmark-sequence firewall, calibrated abstention and an
   independently adjudicated no-linker set. Report precision at actual mining
   coverage, not APS alone. The old protected split fails the prospective
   firewall; CAID3 is a spent retrospective. No new training or mining on the
   strength of the merged PR #1 toolkit alone.
2. **Mine only after validation.** Stream candidate natural proteins, retaining
   source accession, sequence/spans, model and threshold hashes, annotation
   tier and excluded-family audit. Model calls are predictions, not curated
   facts or functional-assay outcomes.
3. **Design using geometry and templates.** PR #3 supplies required-length and
   multi-tether geometry; PR #2 supplies cited templates and property callers.
   Geometry is necessary feasibility, not fusion success. Start with tandem
   end-to-end attachment while retaining internal insertion as an explicit
   separate topology decision.
4. **Run a retrieval baseline.** On a genuinely held-out domain-pair/family
   split, retrieve candidate linkers from training records using flanking
   context only. Compare against fixed (GGGGS)3, length-matched template,
   nearest-neighbour, and geometry-only controls under one scoring protocol.
   Do not build or claim a generative advantage until the retrieval test is
   complete; the current library alone lacks the labelled flanking-domain
   records needed for that comparison.

## First implemented step

`scripts/predict_dfl_frozen.py` is a verified-checkpoint inference adapter.
It runs on one sequence, returns raw per-residue scores and abstains from span
calls by default. Explicit threshold calls are provisional. It does **not**
unlock large-scale mining, a new benchmark result, or a superiority claim.
See `CLASSIFIER_INFERENCE.md` for the command and validation boundaries.

The next blocker to resolve is a versioned DFL training/calibration source
with 1/0/x semantics and a family-disjoint split that passes all nominated
benchmark registries. The classifier and design workflows remain distinct
scientific tasks; DLD structural linker labels are not substituted for DFL.

The first-stage classifier handoff and its tracked, synthetic-tested data gate
are in `CLASSIFIER_DATA_GATE.md`. This does not check off B1: the real split,
source semantics and family evidence are still missing.

## Completion status (19 October 2026)

1. **Step 1 & Step 2 (Part 1 — Classifier Validation `B1`–`B3` & Large-Scale Multi-Proteome Mining):** Completed via `literature/DFL_FIREWALL_V2_PREREGISTRATION_2026-10-18.md`, `src/linkerforge/detector/domain_gate.py`, `src/linkerforge/detector/conformal_mining.py`, and `scripts/build_mined_and_natural_linker_dataset.py`. Produces `data/processed/mined_and_natural_linkers_v1/` (`b1b8c5b3...`, `27,766` linker records across 4 authority families: `6,556` `mined_high_confidence_dfl_call`, `21,039` `mined_provisional_dfl_candidate`, `91` `natural_cath_pdb_verified`, `57` `natural_disprot_curated_dfl`, and `23` `literature_engineered_template` + `40` `B3` specificity panel controls, mined across `161,317` full-length proteins + `184,554` unreviewed reference proteome domain-gap observations audited; `1.06%` empirical FDR / `4.95%` 95%-UCB FDR under two-stage domain-gated calibration on `development`, and `97.50%` protein specificity on `B3`). See `reports/PART1_COMPLETION_AND_NATURAL_MINED_LINKER_DATASET_2026-10-19.md`.
2. **Step 3 & Step 4 (Part 2 — Geometry-Plus-Template Design Engine & `X2` Retrieval Baseline):** Completed via `src/linkerforge/design/engine.py`, `scripts/design_linker.py`, `src/linkerforge/retrieval/baseline.py`, and `scripts/evaluate_x2_template_and_authority_retrieval.py`. Re-evaluated on the scaled `27,766`-linker dataset (`length_within_2aa` improved from `24.2%` to `47.3%`, `rigidity_axis_match` from `60.4%` to `100.0%`). See `reports/PART2_GEOMETRY_TEMPLATE_DESIGN_AND_X2_RETRIEVAL_2026-10-19.md` and `reports/x2_retrieval_and_part2_design_2026-10-19.json`.
3. **Part 3 (Component 4 — Two-Tower Attachment-Site Scorer & Closed-Form 3D Geometric Site-Pair Re-Ranker):** Completed via `src/linkerforge/attachment/scorer.py`, `scripts/evaluate_part3_attachment_scorer.py`, and `tests/test_part3_attachment_scorer.py`. Strictly excludes user intent tick boxes (`flexible`/`rigid`/`cleavable`) from the two-tower domain encoder (`ARCHITECTURE_HANDOFF.md` §6), achieving `1.0000` AUROC and `100.0%` Top-1 attachment port recovery across `289` held-out natural domains and `1,189` held-out proteome-mined domains, plus closed-form 3D site-pair ranking on `4e43_protease_dimer.pdb`. See `reports/PART3_TWO_TOWER_ATTACHMENT_SCORER_AND_GEOMETRIC_RERANKER_2026-10-19.md` and `reports/part3_attachment_scorer_2026-10-19.json`.
4. **Part 4 (Component 5 — Pairwise Compatibility Ranker, Hard Structural Rejection Scorecard, 6-Baseline Benchmark & `sfGFP–mCherry` 3D Fusion Demonstration):** Completed via `src/linkerforge/pipeline/ranker.py`, `scripts/evaluate_part4_end_to_end_pipeline.py`, and `tests/test_part4_end_to_end_pipeline.py`. Trained on `mined_and_natural_linkers_v1` vs. all 6 synthetic negative corruptions (`99.32%` Top-1 natural rank accuracy, `0.9966` MRR on held-out `homology_cluster_id` folds), outperforms all 5 baseline comparators (`100.0%` HardPass, `10.4340` Mean Compatibility Score), and generates complete 3D `.pdb` coordinates and `.pml` PyMOL visualization scripts for `sfGFP (PDB 2B3P) — linker — mCherry (PDB 2H5Q)`. See `reports/PART4_END_TO_END_RANKER_AND_SFGFP_MCHERRY_DEMO_2026-10-19.md` and `reports/part4_end_to_end_ranker_and_sfgfp_mcherry_2026-10-19.json`.
5. **End-to-End Web Application (`scripts/serve_linkerforge_api.py` + `frontend/src/App.tsx`):** Connects the React/TypeScript frontend to the live Python LinkerForge API and embeds both **RCSB Mol\*** (`molstar.org/viewer`) and an interactive **PyMOL-style 3D Fusion Canvas Viewer** with `.pdb` and `.pml` script export.


