# LinkerForge research handoff — 17 October 2026

Read this after `README.md` and `AGENTS.md`. The active scientific task is **sequence-based segmentation of disordered flexible linker (DFL) residues in natural proteins**. The older fusion-linker *design/ranking* specification in README Sections 1–16 is retained as historical/future context; it is not the current training target. DLD's structural independent/dependent-domain linker labels are a separate secondary task and must not silently become DFL ground truth.

**18 October update:** the user resumed work with a classifier-first sequence:
validated DFL caller → audited natural-linker mining → deterministic geometry
plus cited templates → retrieval baseline before any generative design claim.
PRs #1–#7 are merged into this branch. See
`CURRENT_DIRECTION_2026-10-18.md` and `CLASSIFIER_INFERENCE.md`. The old
training firewall and spent-CAID3 rule below are unchanged.

**25 October retrospective update:** at the user's explicit request, the fixed
B1 V2 three-layer CNN was reproduced on its hash-pinned, CAID3-excluded split
and evaluated once more against seven released predictors on the identical
CAID3 common-27 residues. Official APS/AUC-ROC/Fmax were 0.353/0.845/0.423 for
the CNN versus 0.482/0.898/0.491 for IPA-AF2-Linker. See
`reports/DFL_B1_V2_CAID3_RETROSPECTIVE_2026-10-25.md` and its frozen selection
manifest. This is an **additional post-hoc public retrospective**, not a new
untouched test, a mining-precision certificate, or grounds for further CAID3
selection/tuning. The checkpoint and downloaded evaluator inputs remain ignored
local artifacts; the compact official comparison and hashes are committed.

## What exists, and what does not

- There are local detector trainers, preprocessing scripts, experiment configs, tests, and retrospective reports under `scripts/`, `src/`, `configs/`, `tests/`, and `reports/`. The React frontend remains a labelled demonstration, not a deployed detector.
- The protected SAM-like baseline, variants and three-seed repeats were evaluated on filtered development. See `reports/DFL_SAM_REPEAT_AND_LOCKED_CAID3_2026-10-17.md`. Neither the variants nor the baseline demonstrated published-method superiority. The original CAID3 common-27 comparison and the later user-requested B1 V2 retrospective are **spent and retrospective**; never tune or select against either result.
- Current `data/processed/caid3_protected_dfl_v1` is **not eligible for new training**: a prospective sequence-only DP94 screen found overlap with the proposed holdout. `scripts/preflight_caid3_protected_training.py` and `scripts/check_detector_autoresearch.py` fail closed. DP94 is only a candidate, not an independently validated DFL benchmark; verify its label semantics and published-method exposure before considering it.
- The newer self-contained `data/processed/dfl_b1_v2_family_split` passes its
  preregistered sequence/family/hash checks without opening benchmark labels.
  This does **not** authorize a new active round on the present checkout: the
  root protected input `artifacts/caid3-linker.fasta` is ignored and currently
  absent, so `check_detector_autoresearch.py` correctly returns
  `data_integrity_audit_only`. Restore the exact hash-pinned protected input
  through its documented provenance before launching anything; never rebuild
  it from viewed labels or relax the root-registry cross-check.
- The DisoFLAG authors document a DFL annotation track for DP94 (line 9 of
  their released records), so it is inaccurate to call DP94 an IDR-only file.
  Its public availability, negative/unknown semantics, source overlap, and
  competitor-training lineage still prevent a fresh blind-test claim. The
  18 October B1 audit leaves only 24 positive sequences outside both the
  unchanged direct local firewall and the exploratory transitive
  Pfam/Gene3D/UniRef50 graph; see
  `reports/DFL_B1_PAPER_SURVEY_AND_SPLIT_DESIGN_2026-10-18.md`.
- No untouched, homology/domain-family-certified common-target benchmark has yet established a breakthrough. Do not call a development score, public CAID3 comparison, or local-alignment-only exclusion a new blind test.

## Binding scientific rules

1. Follow `literature/DFL_UNTOUCHED_COMPARISON_CONTRACT_2026-10-17.md` for benchmark isolation, family clustering, label sealing, matched published baselines, model freeze and one-time evaluation. Also read `literature/LINKER_DETECTION_FRONTIER_CHARTER_2026-10-11.md` and `literature/DFL_THREE_TRACK_FRONTIER_PROGRAM_2026-10-17.md`.
2. Benchmark **sequences and IDs** may be used for exclusion. Benchmark labels, per-protein errors, scores and thresholds may not be used for training, filtering, feature selection or model selection. Unknown residues are not negatives.
3. Before any new GPU job, create a new versioned family-disjoint train/development/auxiliary split and make its launcher pass a hash-checked preflight against *all* nominated benchmark sequence registries. Do not bypass the present failing preflight. One local GPU job at a time; each round has at most three architecture variants and six GPU-hours.
4. Retrain baseline and proposal on the same data, masks, encoder/input modality and scoring coverage. Select only on development; repeat promising variants across at least three seeds. Freeze model, weights, thresholds, evaluator and competitor coverage before opening a new benchmark's labels.
5. Report negative results, raw counts, hashes, provenance, coverage and failure reasons. Do not auto-promote weights or claim state of the art from a public/previously viewed set.

## Where to start

| Question | Source |
| --- | --- |
| Active task and benchmark policy | `literature/LINKER_DETECTION_FRONTIER_CHARTER_2026-10-11.md`; `literature/DFL_UNTOUCHED_COMPARISON_CONTRACT_2026-10-17.md` |
| Latest frozen result and caveats | `reports/DFL_SAM_REPEAT_AND_LOCKED_CAID3_2026-10-17.md` |
| Three parallel research responsibilities | `literature/DFL_THREE_TRACK_FRONTIER_PROGRAM_2026-10-17.md` and `reports/track-transfer/`, `reports/track-competitors/`, `reports/track-ml/` |
| Current machine-readable status | `configs/detector_autoresearch_active.json`, its named ledger, and `scripts/check_detector_autoresearch.py` |
| Data gate | `scripts/preflight_caid3_protected_training.py` and `scripts/prepare_prospective_dfl_holdout.py` |
| Older training history | `TRAINING.md` and dated reports; do not confuse DLD, DFL, or fusion ranking |

On a machine with the locally configured GPU environment, the read-only state check is:

```powershell
.venv-gpu\Scripts\python.exe scripts\check_detector_autoresearch.py
```

This environment is **not committed** and the command is not a clean-clone installation guarantee. `data/raw/`, `data/processed/`, `artifacts/`, `.venv-gpu/`, checkpoints, caches and downloaded third-party archives are ignored. Source scripts, configs, compact reports and hashes are published; reproduce or obtain source data according to each script's provenance and licence before running. Never substitute placeholder records or silently treat a missing source as a successful run.

## Next executable research work

- Transfer/data track: verify candidate DFL benchmark label semantics, accession lineage and a domain/family-level exclusion and cluster split. Preserve the benchmark label seal.
- Competitor track: document published predictor training corpora, versions, available frozen outputs and whether a truly common-target, equal-coverage comparison is possible.
- ML track: implement and test the new split/preflight and matched development baselines; launch training only after the firewall passes. Preserve the old negative results and one-time CAID3 seal.

The three former 15-minute detector heartbeats were not found as active
automations when the user requested pausing them on 18 October. No replacement
schedule was created. Local schedules and the GPU environment do **not**
transfer through Git; collaborators must not assume a running background job
after cloning.
