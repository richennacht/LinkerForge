# B1 V2 CNN versus CAID3 common-27: fixed retrospective

Date: 2026-10-25 UTC. Status: **post-hoc public retrospective, not an untouched test or an independent superiority claim**. CAID3 scores had been viewed in this project before this run. The run was requested to test whether excluding CAID3-overlapping proteins and retraining the fixed B1 CNN changed the retrospective result; these scores must not drive a subsequent model or threshold choice.

## Frozen procedure and training

- Clean managed checkout of `origin/main` at `0f9ebbbefeb609ea89e9e0421f441a81aee5a709`; ignored input data copied from the primary checkout without edits.
- B1 V2 family/benchmark preflight passed before GPU work. The versioned split has 60 training, 28 development, and 9 auxiliary proteins. The prohibited CAID3 and DP94 registries were checked by the V2 firewall. Split manifest SHA-256: `d62136ce699294a0c07f0fa01a5350d3571643796cedadaf37d8399824893452`.
- Fixed configuration: frozen offline ESM-2 35M residue encoder, three-layer 1/5/1 CNN (50 channels), 510-residue nonoverlapping inference windows, seed 31, 12 epochs. Configuration SHA-256: `c683afbff11088d902603c1f43bf8874ef1b062b7e1444db7856472da46c44ad`.
- Checkpoint selected by maximum **masked development APS** at epoch 4: APS `0.8177586569`, AUROC `0.7813274230`, grid Fmax `0.7258064516` on 547 known positive and 518 known negative development residues. All-window development diagnostic: APS `0.8177616154`, AUROC `0.7813309523`; 13,713 additional residues have unknown labels. The prior 99%-precision diagnostic is **not** a mining-precision certificate; this comparison does not use its threshold.
- Training completed all 12 epochs in 49.032 seconds on NVIDIA RTX 3050 6GB Laptop GPU; peak PyTorch allocated VRAM `186616832` bytes. Best-checkpoint SHA-256: `bc79e37c36cecf4d5b95527ec65d4324372a3a73c64d58cb5b836e1c1f1ceeb3`.
- The [selection lock](DFL_B1_V2_CAID3_RETROSPECTIVE_SELECTION_2026-10-25.json) records the checkpoint, config, split, development audit, scorer, CAID3 reference, seven released competitor files, and official evaluator hashes. It passed hash/lineage validation **before** the CAID3 reference was parsed. The CAID output's binary field uses fixed 0.5 only as an unoptimized display value; APS/AUROC/Fmax use continuous scores.

## Official CAID evaluator, identical common-27 residues

All methods cover the same 27 proteins and 15,570 evaluated residues, including 1,214 positive residues. Values below are from the official evaluator's `all.dataset.f1s.metrics.csv` (APS, AUC-ROC and best-F1 row); local unrounded cross-check is in the machine-readable report. The common-27 intersection is not the full CAID3 challenge.

| Predictor | APS | AUC-ROC | Fmax |
| --- | ---: | ---: | ---: |
| IPA-AF2-Linker | 0.482 | 0.898 | 0.491 |
| LINKER-Pred-Lite | 0.429 | 0.872 | 0.473 |
| DisoFLAG-DFL | 0.410 | 0.846 | 0.464 |
| LINKER-Pred | 0.398 | 0.884 | 0.467 |
| LINKER-Pred2 | 0.386 | 0.889 | 0.464 |
| APOD | 0.369 | 0.852 | 0.438 |
| **LinkerForge B1 V2 CNN** | **0.353** | **0.845** | **0.423** |
| DFLpred | 0.097 | 0.491 | 0.148 |

The B1 V2 CNN did **not** beat the strongest released comparator on any of the three requested metrics. The development-to-CAID3 APS drop is substantial, but these two measurements differ in source, labels, mask prevalence and sampling; this comparison alone does not isolate a causal explanation. No CAID3-based retraining, calibration or architecture choice is authorized by this result.

The committed [compact official metrics snapshot](DFL_B1_V2_CAID3_RETROSPECTIVE_OFFICIAL_2026-10-25.csv) preserves this common-target table. The full evaluator output, prediction files, checkpoint and downloaded reference remain ignored local artifacts under `artifacts/`, in accordance with repository policy; a clean clone must recover those exact hash-pinned inputs separately before it can rerun the scorer. The CSV is a result record, not a substitute for raw predictions or an independent reproduction.

## Independent-claim gate still open

Excluding CAID3 from training is necessary but insufficient for an *untouched* claim: its labels and earlier results were already public and viewed. A new test must be nominated using sequence and source metadata **before** any test labels or outcomes are inspected. Build a versioned graph over every train/development/auxiliary record and prospective test protein using accession aliases, exact/contained sequences, declared local alignment or clustering thresholds, and domain-family evidence. Assign entire connected components to one side and verify that no test family component occurs in training, development, auxiliary data, or retrievable candidates. Independently adjudicate DFL spans and verified no-linker mechanisms with unknown residues masked; seal the label table with a separate custodian. Freeze the model, decision rule, comparator versions, evaluator, target list, and hashes before the custodian scores once. Report protein- and family-level uncertainty and common-target coverage. Until that exists, restrict claims to this post-hoc retrospective.

Artifacts: `artifacts/dfl-b1-v2-cnn-seed31-retrospective-2026-10-25/manifest.json`, `best.pt`, `epochs.jsonl`, `development_all_windows.json`; `artifacts/dfl-b1-v2-caid3-retrospective-2026-10-25/report.json`, `predictions/LinkerForge_B1_V2_frozen.caid`, and `official_evaluator/`.

Methodological context: [CAID challenge results and released Linker reference](https://caid.idpcentral.org/challenge/results); [independent-set sequence splitting study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8929697/); [protein language-model pretraining leakage study](https://proceedings.mlr.press/v261/hermann24a.html). These motivate, but do not themselves certify, the project-specific untouched-test gate above.
