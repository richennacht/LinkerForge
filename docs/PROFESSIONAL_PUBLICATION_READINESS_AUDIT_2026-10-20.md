# LinkerForge professional and publication-readiness audit

Audit initiated: 20 October 2026; verification completed: 21 October 2026 (Asia/Calcutta)
Checkout: `codex/tauri-colabfold-shell`, reviewed-data base commit `43be3ab294c9551bf7b6b32494b71a5442ef419c`, modified working tree

## Executive decision

LinkerForge is suitable for release as an explicitly labelled **computational
prototype** after ordinary code review and packaging work. It is **not ready for
a scientific performance claim, clinical use, or a claim that it designs
experimentally successful linkers**. A software/protocol or work-in-progress
publication is feasible if it reports the limitations and reproducibility
artifacts below. A biological-method paper needs the blocked dataset, ranker,
protected prospective evaluation, matched structural controls and experimental
validation.

There are two distinct research tracks in the repository and they must not be
conflated:

1. The active DFL-detector research is governed by `RESEARCH_HANDOFF.md` and the
   protected benchmark firewall. Its current split fails prospective preflight;
   public CAID3 scores must not be used for model selection.
2. The desktop fusion-design workflow is currently a deterministic
   candidate-enumeration and structural-inspection prototype. It does not contain
   the validated fusion-compatibility ranker specified in the README.

## Changes verified in this audit

Follow-up: [evidence alternatives and future integrations](./EVIDENCE_ALTERNATIVES_2026-10-21.md)
documents the CPU surface-burial extension, MolProbity feasibility, existing
OpenMM/DFLpred alternatives and the remaining scientific and deployment scope.
Run completion is recorded in the referenced durable job artifacts.

- Individual domain structures are no longer fixed sfGFP/mCherry examples. The
  Structure page accepts a sequence-matched PDB independently for Domain A and
  Domain B, plus its source ID and experimental/predicted classification.
- Confidence provenance is explicit: local results are AlphaFold2 predictions
  executed by ColabFold. LinkerForge does not generate pLDDT, PAE or pTM and its
  future compatibility score is a different quantity.
- ColabFold rank-1 score JSON is ingested and checked for exact array dimensions;
  the UI reports region pLDDT, pTM and mean A↔B PAE.
- The UI performs a deterministic C-alpha continuity/overlap pre-screen,
  sequence-matched domain-reference RMSD and domain-relative protected-site
  numbering. These checks do not replace Rosetta or solvent/interface analysis.
- Project v1 save/open and result v1 JSON export are available. Public-MSA
  consent is intentionally not persisted. Desktop prediction directories are
  timestamp-and-sequence-hash addressed, preserving prior runs rather than
  overwriting them.
- Formal schemas now exist for the project, result snapshot, structure-evidence
  checklist and public detector audit. The result export embeds the Git commit
  and whether the build was made from a modified worktree.
- The reviewed GitHub dataset release is present and byte-validated. The
  historical DFL B1 V2 split has 60 training, 28 development and 9 auxiliary
  records and reproduces its preregistered sequence/family preflight. This does
  not authorize a new blind comparison or reopen spent CAID3 labels. The separate 27,766-record
  mined-and-natural-linker release is a source inventory, not a validated
  fusion-compatibility training set.
- Local history is SHA-256 hash chained and the application produces a
  deterministic checksummed research ZIP. This detects accidental or later
  modification but is not an institutional digital signature.
- A native Windows release executable builds with Visual Studio Build Tools and
  the Windows SDK installed on M:. Rust tests exercise the bounded command and
  cross-process GPU-lock behavior.
- Long training, prediction, simulation and comparison chains now use a global
  detached pipeline runner with persistent request, status, stdout and stderr
  records. Required downstream comparisons and schema checks remain in the same
  durable chain when the chat is inactive.

## Software and data verification

- Public-release validator: exact hashes and counts passed for the historical
  DFL B1 V2 split and the 27,766-record composite.
- Historical DFL V2 preflight: reproduced
  `pass_v2_preregistered_sequence_and_family_firewall` without opening benchmark
  labels; CAID3 remains spent and DP94 remains unopened.
- Fusion-ranker preflight: failed closed as `source_inventory_only` with
  `training_authorized: false`, as intended.
- Python: 417 tests passed, 6 dependency-aware tests skipped because their
  non-versioned protected CAID3/DisProt/DLD source artifacts are absent, and 35
  subtests passed. The current self-contained V2 firewall tests were executed.
- Frontend after a clean `npm ci`: lint passed, all 23 tests passed and the
  production build completed. The only build notices are the upstream 3Dmol.js
  direct-`eval` warning and its isolated >500 kB viewer chunk.
- Native Rust: formatting passed and all 6 tests passed, including the
  cross-process GPU lock.
- JSON contracts: all 69 schemas passed Draft 2020-12 meta-validation; the
  generated detector, protected-site, structure-evidence, ranker-preflight and
  exploratory-ranker records also validate against their respective schemas.

## Completed structure evidence: public 1UBQ–linker–2GB1 run

Artifact: `runs/folds/acceptance-20261020-ubiquitin-gb1-target10/evidence_checklist-linked-20261021.json`

The complete 142-residue construct contains human ubiquitin (PDB 1UBQ),
`GGGGSGGGGS`, and the protein G B1 domain (PDB 2GB1). The references were
downloaded from RCSB and exact sequence matching was required before RMSD.

| Check | Result | Interpretation |
| --- | ---: | --- |
| Full-fusion sequence | Exact 142-residue match | Passes identity gate |
| Whole-chain pLDDT | 81.74 | AlphaFold2/ColabFold local confidence, not LinkerForge confidence |
| Domain A pLDDT | 91.10 | High local-coordinate confidence |
| 10-aa linker pLDDT | 37.73 | Very low local confidence; compatible with flexibility but not proof of it |
| Domain B pLDDT | 76.90 | Moderate local-coordinate confidence |
| pTM | 0.53 | Global arrangement is not strongly supported |
| Mean cross-domain PAE | 22.77 Å | Relative domain orientation is highly uncertain |
| C-alpha chain breaks / severe overlaps | 0 / 0 | Passes only the deterministic pre-screen |
| Maximum adjacent C-alpha distance | 4.04 Å | No >4.5 Å pre-screen break |
| Domain A preservation | 0.94 Å C-alpha RMSD | Close to the experimental reference under this metric |
| Domain B preservation | 3.02 Å C-alpha RMSD | More deformation/difference; requires review |

The fold is a successful local compute and visualization acceptance test. The
low linker confidence and high cross-domain PAE mean it must not be presented as
a uniquely determined interdomain arrangement or as evidence of function.

## Published-linker detector sensitivity comparison

Artifacts: `runs/detector-audits/public-fusion-panel-20261020/`

The exact public DFLpred implementation at revision
`b62a4b9543d6fe25b4230da8bfd1ed66f0071d8f` was run at its documented 0.18
threshold. The panel has four public-sequence contexts—1UBQ/2GB1, 2B3P/2H5Q,
1CRN/4PTI and 1PGA/1LMB—and seven junction choices, for 28 constructs.

| Junction | Length | Mean inserted-span call fraction | Range across contexts | Mean domain call fraction |
| --- | ---: | ---: | ---: | ---: |
| Direct fusion | 0 | N/A | N/A | 0.073 |
| `GGGGS` | 5 | 0.550 | 0.000–1.000 | 0.140 |
| `GGGGSGGGGS` | 10 | 0.950 | 0.800–1.000 | 0.151 |
| `(GGGGS)3` fixed control | 15 | 0.883 | 0.533–1.000 | 0.141 |
| `GSAGSAAGSGEF` | 12 | 0.667 | 0.250–1.000 | 0.107 |
| `KESGSVSSEQLAQFRSLD` | 18 | 0.361 | 0.000–0.833 | 0.094 |
| `EAAAKEAAAKEAAAK` rigid control | 15 | 0.017 | 0.000–0.067 | 0.088 |

This answers only whether a published DFL detector recognizes the deliberately
inserted span. It is not an independent accuracy benchmark: the inserted spans
are known by construction, not experimentally annotated DFLs in those contexts.
The wide context ranges and 7.3–15.1% mean domain calls show why a positive call
cannot be treated as compatibility, expression, activity or structural success.
The rigid EAAAK negative result is expected for a detector specifically intended
for disordered flexible linkers.

## Exploratory historical ranker reproduction

Artifact: `runs/ranker-audits/part4-reproduction-20261020.json`

The historical 12-feature pairwise ranker was reproduced against six
rule-generated corruptions per record. Its record is schema-validated and now
states `exploratory_reproduction_not_validation` and
`training_authorized: false`. It is not the README-frozen ESM-2 plus
attachment-point geometric-GNN model.

On the 148 natural-authority records, the natural sequence ranked first in
99.32% of these synthetic seven-way comparisons (MRR 0.9966). This high number
is not an external biological result: the negatives were generated by rules
that overlap the hand features being scored. The mean score comparison was:

| Historical computational baseline | Mean score | Hard-pass fraction | Win rate versus generated negatives |
| --- | ---: | ---: | ---: |
| Direct fusion | -2.1164 | 0.0000 | 0.0000 |
| Fixed `(GGGGS)3` | 8.8044 | 0.8514 | 0.9403 |
| Random length-matched GS | 10.1562 | 1.0000 | 0.9977 |
| Geometry only | 9.3986 | 1.0000 | 1.0000 |
| Nearest-neighbour retrieval | 9.7191 | 1.0000 | 0.9966 |
| Historical LinkerForge hand-feature ranker | 10.4340 | 1.0000 | 1.0000 |

The small score separation from random length-matched GS and the saturated
synthetic-negative metrics reinforce the need for a curator-reviewed,
family-disjoint compatibility dataset and an untouched prospective evaluation.
These numbers must not appear in the user workflow as a calibrated probability
or a claim that the generated linker is best.

## Matched ColabFold structural sensitivity panel

Artifact: `runs/folds/matched-baselines-20261020-1ubq-2gb1/comparison_summary.json`

Seven complete 1UBQ-linker-2GB1 constructs, including direct fusion and the
existing 10-residue acceptance run, were compared with the same
AlphaFold2-ptm model set, seed, MSA mode, template setting and pLDDT ranking.
Each newly submitted candidate completed all five models. The values below are
from each candidate's rank-1 model; they are matched computational sensitivity
measurements, not biological success labels.

| Candidate | Linker length | Overall pLDDT | Linker pLDDT | pTM | Cross-domain PAE (Å) | Domain A RMSD (Å) | Domain B RMSD (Å) | DFLpred linker recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Direct fusion | 0 | 88.25 | n/a | 0.61 | 16.74 | 1.00 | 2.96 | n/a |
| `GGGGS` | 5 | 87.03 | 48.57 | 0.56 | 21.14 | 0.80 | 2.84 | 0.20 |
| `GGGGS` x2 | 10 | 81.74 | 37.73 | 0.53 | 22.77 | 0.94 | 3.02 | 0.80 |
| Fixed control `(GGGGS)3` | 15 | 82.75 | 38.99 | 0.53 | 23.40 | 0.87 | 2.79 | 1.00 |
| WALDO-style | 12 | 84.09 | 44.42 | 0.54 | 21.85 | 0.94 | 1.88 | 0.58 |
| BIRD-style | 18 | 80.21 | 44.87 | 0.52 | 24.34 | 1.17 | 2.96 | 0.56 |
| Rigid `EAAAK` x3 | 15 | 87.88 | 50.88 | 0.51 | 24.90 | 0.79 | 1.07 | 0.00 |

All seven rank-1 coordinates passed the same limited coordinate pre-screen:
zero C-alpha chain breaks and zero severe C-alpha overlaps. Direct fusion had
the highest overall pLDDT and pTM and the lowest cross-domain PAE, while the
rigid control preserved Domain B best by reference RMSD. Neither observation
means either construct is the desired flexible biological linker. Conversely,
DFLpred called every residue of `(GGGGS)3` and none of rigid `EAAAK`, showing
that detector sensitivity and fold confidence measure different properties.
No winner was selected because acceptance thresholds were not preregistered.
This single public domain pair does not establish family-disjoint
generalization, expression, activity, stability or an optimal linker.

## OpenMM execution smoke test

Artifact: `runs/openmm/acceptance-1ubq-gs10-2gb1-smoke/run_manifest.json`

The accepted 1UBQ-`GGGGSGGGGS`-2GB1 rank-1 prediction completed three explicit-
solvent OpenMM smoke replicas on CUDA using AMBER14/TIP3P-FB. Each replica used
10 ps equilibration followed by 50 ps production, with independent seeds
20261020–20261022. PDBFixer added only the missing terminal `OXT` atom and
hydrogens at pH 7; it added no residues or non-terminal heavy atoms. The
solvated system contained 36,751 atoms.

| Replica | Final potential energy (kJ/mol) | Final fusion C-alpha RMSD (Å) | Final Domain A RMSD (Å) | Final Domain B RMSD (Å) | Minimum linker/domain C-alpha distance observed (Å) | Maximum linker/domain C-alpha contacts |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | -610784.79 | 1.69 | 0.72 | 0.86 | 3.83 | 7 |
| 2 | -610590.53 | 1.96 | 0.82 | 1.25 | 5.81 | 3 |
| 3 | -608768.17 | 1.96 | 0.99 | 0.89 | 6.06 | 4 |

This establishes that preparation, minimization, replicated CUDA execution,
trajectory/state capture and structural-metric reporting work end to end. The
50 ps trajectories are far too short for convergence, equilibrium sampling or
a stability claim. The versioned 3 x 50 ns configuration remains the intended
reviewed production protocol and has not been executed here.

Operational preflight artifact:
`reports/openmm-production-preflight-2026-10-21.json`. Using the observed
same-construct smoke trajectories and production timestamps, it projects
15,000 stored frames, 6.162 GiB of DCD output and 31.76 hours of GPU execution.
The frozen launch policy requires 1.5 times projected trajectory storage plus
a 5 GiB recovery reserve. M: currently has 11.887 GiB free versus 14.243 GiB
required, so the detached production job is fail-closed with a 2.357 GiB
shortfall. This is an operational capacity decision, not scientific evidence.

A separate storage-aware v2 protocol has now passed a fresh preflight without
altering the simulated system or duration. It preserves 10 ps structural
metrics/state/checkpoints and stores full-system DCD coordinates every 50 ps.
The projection is 3,000 DCD frames, 1.232 GiB of trajectory storage and 31.76
GPU hours; its unchanged safety gate requires 6.849 GiB versus 10.545 GiB free.
Artifact: `reports/openmm-production-preflight-storage-aware-2026-10-21.json`.
This removes only the storage blocker. A durable production pipeline is queued
behind the active ColabFold comparison and will repeat the preflight, acquire
the shared GPU lock, execute three replicas, summarize descriptive five-block
and between-replicate diagnostics, link the evidence and validate schemas.
Until all those stages finish successfully, production evidence is incomplete.

## Method and product comparison

| Method | What it measures | Local status | What LinkerForge should use it for |
| --- | --- | --- | --- |
| DFLpred | Per-residue disordered flexible-linker propensity | Exact public implementation run on 28 constructs | Orthogonal sensitivity flag only |
| CAID3 released submissions | Per-residue DFL detection on an identical 27-protein common set | One-time retrospective completed: best historical LinkerForge APS/Fmax 0.089/0.165, below IPA-AF2-Linker 0.482/0.491, LINKER-Pred-Lite 0.429/0.473 and DisoFLAG-DFL 0.409/0.464 | Negative historical evidence only. The public set is spent, has no no-linker proteins, cannot guide tuning/selection and does not measure fusion compatibility |
| LINKER-Pred / Lite | Per-residue linker-region prediction using protein embeddings | Repository inspected; not reproduced because the Lite workflow requires external 1024-dimensional ProtTrans embeddings/model assets and a sealed common-target evaluation | Future detector comparator on a sealed evaluation set |
| TransDFL, DisoFLAG-DFL | Modern DFL/linker detectors | Not reproduced locally | Future frozen comparator; never use spent CAID3 results for selection |
| ColabFold / AlphaFold2 | Complete-sequence structure plus pLDDT, PAE and pTM | Operational; seven matched public constructs and 30 newly computed models plus the five-model existing run | Structural hypothesis and uncertainty, not linker success label |
| OpenFold | Reproduction/training implementation of AlphaFold2 | Not installed | Reproducibility/engineering alternative; not independent biological evidence |
| ESMFold | MSA-free protein structure prediction | Not installed; original repository archived | Orthogonal single-sequence sensitivity if run on adequate hardware |
| AlphaFold 3 | Broader biomolecular complex prediction | Not a local 4 GB workflow | Approved service/cloud sensitivity study where data governance permits |
| Boltz / Chai-1 | Open biomolecular structure models | Not installed; hardware envelope not established | Cloud/HPC comparator for selected finalists |
| Rosetta KIC/loop modelling | Linker conformational sampling, closure and all-atom scoring | Not installed | Required structural-funnel stage, under a frozen protocol |
| Reduce / Probe | Hydrogen-aware fixed-coordinate severe-overlap contacts | Pinned CPU environment and seven-prediction run completed; GGGGS10 has 48 unique severe pairs, 22.22 per 1000 explicit atoms, worst overlap 1.473 Å | Complementary geometry cross-check only; not canonical MolProbity clashscore, Rosetta sampling or a candidate acceptance test |
| OpenMM / GROMACS | Energy minimization and molecular dynamics | Three-replica explicit-solvent OpenMM CUDA smoke run completed; storage-aware 3 × 50 ns v2 preflight passes and a durable pipeline is queued behind the active ColabFold GPU panel | Optional reviewed-finalist validation, with force field, solvent, replicates and convergence records |
| PyMOL | Inspection, measurement and figure generation | Open-source PyMOL 3.1 environment on M: rendered the public acceptance structure headlessly; LinkerForge's distribution policy does not bundle PyMOL | Visualization only; never call it simulation; any future bundling requires a fresh package-specific licence review |

## Professional evidence checklist

| Requirement | Status after this audit | Release consequence |
| --- | --- | --- |
| Mined linker corpus, independent review and trained compatibility ranker | **Partial / ranker blocked**: the 27,766-record composite is a valid provenance-stratified discovery, retrieval and weak-supervision corpus (27,595 classifier-mined candidates, 148 natural authorities and 23 literature templates). Whole-corpus manual validation is neither expected nor required. A deterministic 971-record dual-curator review packet now targets linker identity and source evidence, but no completed adjudication or separate engineered-fusion outcome dataset exists and the frozen ESM-2/geometric-GNN compatibility ranker is not trained | The mined corpus may support retrieval and weak supervision; no personalized or “best linker” claim |
| PAE ingestion/interdomain analysis | **Implemented** | Report source and uncertainty; no frozen acceptance threshold yet |
| Rosetta closure/clash | **Rosetta future; complementary screen complete**: licensed Rosetta is unavailable, while pinned Reduce/Probe completed hydrogen-aware fixed-coordinate overlap measurements on seven public predictions | Cannot claim Rosetta linker reconstruction, conformational sampling or closure acceptance; report Reduce/Probe only under its narrower claim boundary |
| Reference-domain preservation | **Implemented and exercised** | C-alpha RMSD is useful but not sufficient |
| Protected-site mapping | **Implemented as a coordinate pre-screen**: domain-relative numbering, minimum linker distance and contact count; a public K48 example also reports SASA | Do not equate this with ligand/interface annotation or all-atom site preservation |
| Remaining baselines | **Partial**: sequence-detector panel and a seven-construct matched ColabFold panel on one public domain pair are complete; other structure-model families and family-disjoint targets remain absent | Cannot claim superiority over published/control linkers |
| Project save/open | **Implemented v1** | Reference PDBs must be reloaded by design |
| Immutable run history | **Partial**: timestamp/hash desktop directories and an in-app SHA-256 hash chain are implemented; no digital signature or institutional append-only store | Preserve exported history and verify its chain; do not call it cryptographically signed |
| Versioned schemas/result export | **Implemented**: v1 project/result/history/evidence schemas, deterministic checksummed research ZIP and schema meta-validation in CI | Retain schema versions and migrate explicitly |
| Institutional/cloud handoff | **Not performed**: no credentials or approved endpoint supplied | Local 4 GB failures need manual handoff |
| OpenMM/GROMACS workflow | **Smoke complete; production v2 queued**: explicit solvent, AMBER14/TIP3P-FB, three 50 ps smoke replicas and independent seeds are verified. The 3 × 50 ns v2 protocol preserves 10 ps metrics while storing DCD every 50 ps; it projects 1.232 GiB, 31.76 GPU hours and passes the 6.849 GiB safety gate with 10.545 GiB free. A durable pipeline is waiting for ColabFold, then repeats preflight and runs all downstream diagnostics | The smoke proves execution only. The v2 evidence remains incomplete until every queued stage succeeds; even completion is not experimental stability or function evidence |
| Experimental assays | **Outside this software audit by user direction** | Biological efficacy remains an external claim boundary |

### External P-Linker archive disposition

The separately supplied `Linker.zip` does not close the compatibility-data
blocker. Its 7,892 silver records are rule-mined Pfam gaps without fusion
outcomes; 7,890 join the same Pfam family, 5,057 involve seven-transmembrane
families, and 2,843 fall outside the frozen 5–30 aa linker range. All ten
literature/PDB records fail their stored boundary-consistency invariant. The
archive is therefore quarantined and was not merged. The reproducible findings
are in `reports/external-plinker-dataset-audit-2026-10-21.json` and
`docs/EXTERNAL_PLINKER_DATASET_AUDIT_2026-10-21.md`.

## Discrepancies and remaining publication blockers

1. The README describes the future fusion-compatibility ranker while the active
   research handoff concerns a DFL detector. Publication must separate these
   objectives, datasets, metrics and claims.
   The self-contained B1 V2 split currently passes its preregistered
   sequence/family firewall, but the active launcher remains fail-closed because
   the ignored root protected CAID3 FASTA is absent. This is a reproducibility
   input failure, not permission to bypass the cross-check or call DP94 an
   independently certified benchmark.
2. The labelled preset UI contains historical illustrative ranking numbers.
   They remain examples, not scientific records; a release build should either
   remove them or place them in a separate demo fixture with conspicuous
   watermarks.
3. The matched ColabFold panel now compares direct fusion, three GS lengths,
   rigid `EAAAK` and two literature-derived sequences under common settings,
   but it covers only one public domain pair and had no preregistered winner
   threshold. Family-disjoint target panels and independent structure-model
   sensitivity runs remain necessary for comparative performance claims.
4. Drive M: now has approximately 12.8 GB free. This is sufficient for the
   bounded checks in this audit but not a reason to install another large
   predictor/database without a storage plan; the GPU still has only 4 GB VRAM.
5. Public MMseqs2 use sends sequences off the workstation. Institutional data,
   confidential constructs and sequences with contractual restrictions require
   an approved private search/compute route.
6. PyMOL redistribution is no longer an installer blocker under the frozen
   non-bundling policy: LinkerForge ships neither Incentive nor Open-Source
   PyMOL and treats a separately installed executable as optional. Upstream
   Open-Source PyMOL uses a permissive BSD-like licence with notice and
   trademark conditions; any future bundled runtime still needs a fresh review
   of the exact package and dependency set.
7. The default Tauri configuration deliberately produces a local release
   executable with bundling disabled. A separate configuration now builds an
   unsigned, current-user NSIS preview for packaging verification. A public
   installer still needs the real publisher identity, code signing, update
   policy, final licence notices and a supported matrix; those cannot be
   fabricated from the repository.

## Minimum credible publication plan

1. Freeze the product claim as “template enumeration, full-fusion prediction
   orchestration, provenance capture and evidence review,” not validated linker
   design.
2. Remove or isolate illustrative scores; add CI schema validation and a
   checksummed research-bundle archive.
3. Complete the prospective DFL firewall remediation before training or model
   selection. Publish the sealed split protocol and all exclusions.
4. Complete and adjudicate the independent linker-identity review packet without
   relabelling mined predictions as compatibility outcomes. Separately build the
   fusion-compatibility dataset with source IDs, evidence classes,
   homology/domain-family clustering and synthetic-negative labels kept distinct.
5. Train/calibrate the specified pairwise ranker and evaluate against every
   baseline under identical rules on an untouched family-disjoint test.
6. Run the matched structural funnel for all shortlisted candidates and controls
   on an approved cloud/HPC image; retain failures, versions, hashes and seeds.
7. Treat expression, stability and activity assays as external evidence that is
   outside this software delivery. Do not support biological-performance claims
   until an independent experimental programme supplies that evidence.

## Primary external references

- DFLpred: <https://github.com/fanchi/DFLpred>
- LINKER-Pred: <https://github.com/deemeng/punch_linker>
- LINKER-Pred Lite: <https://github.com/deemeng/punch_linker_light>
- CAID methods: <https://caid.idpcentral.org/methods>
- DisoFLAG: <https://github.com/YihePang/DisoFLAG>
- TransDFL: <https://bliulab.net/TransDFL/>
- ColabFold: <https://github.com/sokrypton/ColabFold>
- AlphaFold 3 input specification: <https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md>
- Boltz: <https://github.com/jwohlwend/boltz>
- Chai-1: <https://github.com/chaidiscovery/chai-lab>
- OpenFold: <https://github.com/aqlaboratory/openfold>
- ESM/ESMFold: <https://github.com/facebookresearch/esm>
- Rosetta KIC loop modelling: <https://docs.rosettacommons.org/docs/latest/application_documentation/structure_prediction/loop_modeling/loopmodel-kinematic>
- OpenMM user guide: <https://docs.openmm.org/latest/userguide/>
