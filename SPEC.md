# LinkerForge: project specification (state as of 2026-10-17)

Read this first. It is written for an agent with no memory of the project. It states what is verified, what is not, what must not be claimed, and what has to happen before a publishable result exists. The step-by-step work order is in `CHECKLIST.md`, next to this file.

Status update (18 October 2026): work resumed on explicit user instruction. PRs #1–#7 have been merged into `codex/dfl-research-handoff`. This document remains a 17 October evidence snapshot; read `CURRENT_DIRECTION_2026-10-18.md` for the current classifier-first, geometry-plus-template and retrieval-baseline programme. Scientific gates below remain binding.

Sources of truth, in order: `RESEARCH_HANDOFF.md` on `codex/dfl-research-handoff`, the dated reports under `reports/`, then this file. If this file disagrees with a dated report, the report wins and this file needs fixing.

---

## 1. What the project is

LinkerForge is a repository (`richennacht/LinkerForge`, active branch `codex/dfl-research-handoff`) with two strands of work.

**Strand A (active): sequence-based detection of disordered flexible linker (DFL) residues in natural proteins.** A residue-level classifier that labels which residues are flexible linkers, using sequence input and a frozen protein language model encoder (ESM-2). Evaluated on DisProt-derived labels and the CAID3 common-27 benchmark. `RESEARCH_HANDOFF.md` names this as the active task.

**Strand B (downstream goal, not yet buildable): a fusion linker design model.** The classifier exists to mine linkers. Its output segments unannotated proteins (UniProt-wide in the long run). The mined linkers then train a design model with three questions: Q1 is feasibility given a linker specification, Q2 is attachment sites, Q3 is the linker sequence. Linker type is set by two tick-boxes, cleavable or non-cleavable and rigid or flexible. The older fusion design and ranking spec in README Sections 1 to 16 is historical context. Strand B is downstream of Strand A.

Consequence: precision on the mined subset matters more than headline APS, and abstention is a feature, not a failure.

DLD structural independent and dependent domain linker labels are a separate, secondary task. They must not silently become DFL ground truth.

---

## 2. Verified state

All numbers below come from committed files. Paths are relative to the repo root.

### 2.1 Frozen baseline on protected development data

- Matched SAM-like baseline, three seeds (7, 17, 31), evaluated on 7,427 known development residues (1,267 positive). Mean APS **0.7918**, mean AUROC 0.9076, mean Fmax 0.7022. Source: `reports/DFL_SAM_REPEAT_AND_LOCKED_CAID3_2026-10-17.md`.
- Two variants (multiscale features; boundary plus query) beat baseline mean APS by +0.0067 and +0.0037 but failed the preregistered no-positive-protein specificity limit (loss must stay within 0.02). The frozen choice is the baseline seed-7 checkpoint. Rule: `literature/DFL_SAM_REPEAT_SELECTION_RULE_2026-10-17.md`. Freeze manifest: `reports/dfl_sam_repeat_frozen_selection_2026-10-17.json`.
- Neither variant demonstrates superiority. Preserve these negative results.

### 2.2 One-time locked CAID3 common-27 comparison (SPENT)

27 proteins, 15,570 residues, 1,214 positive. Official CAID evaluator output in `reports/dfl_sam_repeat_locked_caid3_2026-10-17/`.

| Predictor | APS | AUROC | Fmax |
|---|---|---|---|
| IPA-AF2-Linker | 0.482 | 0.898 | 0.491 |
| LINKER-Pred-Lite | 0.429 | 0.872 | 0.473 |
| DisoFLAG-DFL | 0.409 | 0.846 | 0.464 |
| LINKER-Pred | 0.398 | 0.884 | 0.468 |
| LINKER-Pred2 | 0.386 | 0.889 | 0.464 |
| APOD | 0.369 | 0.852 | 0.438 |
| **LinkerForge frozen SAM baseline** | **0.241** | **0.695** | **0.345** |
| DFLpred | 0.097 | 0.491 | 0.148 |

The project baseline is second from last. CAID3 is public, was inspected earlier, and is now spent. It is retrospective. Nothing may be selected, tuned, or thresholded using CAID3 again.

### 2.3 Mask-policy result (the strongest asset)

Same model scores, two evaluation masks. The TE82 masked CNN scores 0.450 APS when unknown residues are excluded. The same scores give 0.115 APS when unknown residues are counted as negative. The audit table lists 13,694 residues (16.23%) under the original mask and 36,938 (6.02%) under the zero-as-negative mask. The column meanings should be checked against the script before quoting them. Source: `reports/DFL_POOR_TRANSFER_CAUSE_AUDIT_2026-10-12.md`, rows at lines 20 and 21.

- The 93.35% unknown-label inventory is the reason. Most residues in the DisProt-derived annotations are unannotated, and the pipeline treats them as negatives. Source: `reports/DFL_LABEL_MASK_AUDIT_2026-10-12.md`.
- The model did not change. Mask policy alone moves reported APS by about 0.33. So mask policy, not model quality, is the dominant term in reported performance on this benchmark.
- This result needs no GPU and no new split. It is the strongest thing the project owns.
- Caveat: the frozen score files that support the reproduction are gitignored under `artifacts/`. A clean clone cannot recompute 0.450 and 0.115 without them. The committed report and JSON are the evidence until the scores are reproduced.

### 2.4 Corrected numbers (both were wrong when relayed)

- The field's best published result on CAID3 common-27 is **APS 0.482, Fmax 0.491** (IPA-AF2-Linker). It is not "F1 0.38". Source: section 2.2 above.
- The chain extension bound is **3.8 Å per residue**, measured as the virtual C-alpha to C-alpha distance. It was verified over 412 residues of real structure in PR #3 (`linkerforge.geometry`). The 3.0 to 3.5 Å figure that two literature surveys gave is a **design convention**. Used as a bound it makes real protein segments look impossible. Helical rigid linkers have an axial rise of about 1.5 Å per residue.
- General lesson, recorded in the claims register: a number that is conventional in a design context but physical in form does not say which it is, and nobody downstream re-derives it. Re-derive before quoting.

### 2.5 Data state

- `data/processed/caid3_protected_dfl_v1` is **not eligible for new training**. A prospective sequence-only screen (DP94) overlaps the proposed holdout.
- `scripts/preflight_caid3_protected_training.py` and `scripts/check_detector_autoresearch.py` fail closed on this. Do not bypass them.
- DP94 is a candidate benchmark, not an independently validated one. Its label semantics and published-method exposure are unverified.
- Splits so far are by protein ID, not by family cluster. No family-disjoint split exists yet.
- `data/`, `artifacts/`, `.venv-gpu/` and checkpoints are gitignored. A fresh clone does not contain the data needed to train or rescore.

### 2.6 Hardware and rules

- Dorya's laptop: RTX 3050, 6 GB, no torch installed at last check. Long jobs belong on the college cloud. Free tiers are for feature extraction only.
- Per round: one local GPU job at a time, at most three architecture variants, six GPU-hours. Source: `RESEARCH_HANDOFF.md`, binding rule 3.
- Frozen encoders mean embeddings can be extracted once on free cloud GPUs and small heads trained locally.
- Benchmark sequences and IDs may be used for exclusion. Benchmark labels, per-protein errors, scores and thresholds may not be used for training, filtering, feature selection or model selection. Unknown residues are not negatives in any new claim.

---

## 3. Pull requests (drafts when this snapshot was written; merged 18 October)

Base branch for #1 to #3 is `codex/dfl-research-handoff`. Dorya is not a repo collaborator, so GitHub review requests fail. The repo is the system of record.

| PR | Head branch | What it contains | Tested | NOT verified |
|---|---|---|---|---|
| #1 | `claude/dfl-high-precision-toolkit` | nnPU loss (`class_prior` null, launcher refuses to run), ensemble scoring, Clopper-Pearson precision threshold (`detector/selective.py`), AlphaFold feature loader, span decoders. Status file: `reports/track-ml/2026-10-17-high-precision-toolkit-status.md` (on that branch, not on this one) | 43 CPU tests pass, 2 skipped | torch loss and gradients, PU launcher wiring (py_compile only), ensemble CLI (never run), any run on real scores or real AlphaFold data, full repo test suite. The Clopper-Pearson bound is residue-level and optimistic, and it is selective prediction with a binomial bound, not conformal risk control |
| #2 | `claude/linker-property-library` | 30-entry template library of engineered linkers from Chen et al. 2013 (doi:10.1016/j.addr.2012.09.039), two-axis property caller, coverage audit | 28 tests | Literature-reported use only, not function-tested. Source text has internal inconsistencies, flagged per entry. **Not a contribution:** SynLinker (Liu et al. 2015, doi:10.1093/bioinformatics/btv447) shipped equivalent curation in 2015 |
| #3 | `claude/linker-geometry-feasibility` | `linkerforge.geometry`: required linker length per attachment pair, multi-tether satisfiability (proof, witness or `undetermined`), chain topology | 80 tests, standard library only, runs on the laptop | Geometric admissibility is necessary, not sufficient. Says nothing about expression, folding, domain compatibility or which loops tolerate interruption. Novelty of 3+ simultaneous tethers is an unsearched-literature negative result |
| #4 | `claude/uniprot-mining-sizing` | `MINING_RUN_SIZING.md`: UniProt 2026_03 counts, storage arithmetic, Swiss-Prot first recommendation, leakage control design, resumable pipeline | Arithmetic only | AlphaFold DB coverage, current AFDB entry count, encoder throughput, TrEMBL and UniRef50 mean lengths |
| #5 | `claude/session-handoff-claims-register` | `SESSION_HANDOFF_2026-10-17.md`: claims register and publication plan, with links to the two artifact pages | Documentation only | Its own content is a summary of the reports above |

PR #1 and PR #5 both say that checkpoint and score artifacts are not in git. That is expected.

### CI gap (affects every PR into the active branch)

`.github/workflows/ci.yml` triggers `pull_request` only for `main` and `integration`. **PRs into `codex/dfl-research-handoff` run no CI at all.** An absence of red checks there is not a pass. There is also no Python job: `tests/` never runs in CI. Python changes must be verified locally and the exact command stated in the PR. Fix belongs in its own PR.

---

## 4. Publication plan

Two papers. Paper I first.

### Paper I: benchmark and methods analysis (UNBLOCKED)

Needs no laptop, no GPU, no new split and no new data. Target: a bioinformatics methods journal, not an ML main track. The contribution is the analysis, not a method. The abstract must say that no new method is introduced.

Contents:
1. Equal-coverage CAID3 reanalysis. The published table mixes predictors scored on 90.3% of residues with full-coverage ones. Recompute on a common coverage.
2. Mask-policy counterfactual (section 2.3), reproduced from frozen scores.
3. The 93.35% unknown-label inventory, and what it implies for any negative-labelled benchmark.
4. Input modality audit: the LINKER-Pred family is not sequence-only at inference. This needs verification against the tools' own documentation before it is stated.
5. About ten clean negative ablations (multiscale features, boundary plus query, context span, bounded span, and others under `reports/`).

Paper I also supplies the argument Paper II needs: that abstention is a conclusion, not an apology.

### Paper II: certified high-precision caller plus mined set (GATED)

Do not draft Paper II before the first prerequisite passes. Drafting creates pressure to produce the number the draft is shaped around. Prerequisites, in order:
1. A family-disjoint split that passes the hash-checked preflight against all nominated benchmark sequence registries.
2. Conformal risk control with abstention (Learn then Test), with thresholds calibrated on a small clean protein-split calibration set. Only conformal methods certify a precision level. PU, noisy-label methods and ensembles improve the model but still leave a hand-picked threshold.
3. An adjudicated no-linker set. Without it there is no specificity claim, and therefore no mining claim.
4. An untouched test set, opened once, after freezing model, weights, thresholds, evaluator and competitor coverage.

---

## 5. Integrity constraints (violating these turns a result into a retraction)

1. **The precision-0.90 operating point is a threshold selected on the split it is reported on.** It is not a guarantee. Never call it one, and never put "guarantee" or "certified" near it, especially in a title.
2. **0.7918 and 0.241 are different quantities.** 0.7918 is development APS of the SAM baseline. 0.241 is locked CAID3 APS of the same checkpoint on a different benchmark. Never place them in the same table, sentence or abstract without a clear label and a reason.
3. **CAID3 is spent.** Anything selected, tuned, thresholded or filtered using CAID3 labels or scores is not a result. This includes choosing the common-27 subset for convenience after seeing it.
4. **The template library is not a contribution.** SynLinker shipped 2,260 curated linkers with length, composition, solvent-accessibility and proteolytic-site filters in 2015, plus fusion structure prediction from two uploaded domains. Build the library to use, not to publish.
5. **Li et al. 2026 (doi:10.1016/j.synbio.2026.02.003) figures are permanently unverified.** The PMC copy is captcha-walled, ScienceDirect blocks robots, and there is no university library access. Do not use its numbers as a benchmark or comparison. This includes its structure-versus-sequence ablation.
6. **Do not use the word "first" near multi-tether work.** Two connection points between the same domains is established practice (GPCR fusion, Baker-lab two-rigid-junction work). No precedent was found for three or more simultaneous peptide linkers, but the crosslink and staple literature was never searched. Search it before any novelty word appears.
7. **pLDDT is not a validity score for linkers.** Low pLDDT is the disorder signal (CAID uses 1 minus pLDDT). Filtering generated linkers on structural confidence rejects correctly disordered ones. A flexible linker has no single reference backbone, so self-consistency RMSD does not apply.
8. **Never report a success probability for a fusion.** The best zero-shot viability predictor is ESMFold pLDDT at AUC 0.73 (all metrics combined 0.74), and no dataset maps fusion constructs to expression outcomes. Report flags, not probabilities.
9. **Breadth is not quality.** Mined natural linkers have no functional labels. Library size is not a quality claim.
10. **Unknown is not negative** in any new label, metric or claim.
11. **Do not manufacture records** from advertised database counts. SynLinker's home page gives no usable data.

---

## 6. The classifier target is not reachable

Dorya's original target was 90% recall at 99%+ precision. It is not reachable with the current labels. The ceiling is the labels, not the model. Unannotated residues count as negatives, so many apparent false positives are real linkers that nobody annotated. Mask policy (section 2.3) is the evidence. Architecture is not the bottleneck: the best linker model in the repo is a 3-layer CNN, and SAM and Mask2Former decoders were negative in-repo. Data is the lever.

Keep reporting real numbers and say plainly where the target sits. Changing experiment parameters to reach a number is not an option.

---

## 7. Design-model niche (Strand B), briefly

- Q1 (feasibility) mostly reduces to a calculator. PR #3 is that calculator.
- Q2 (attachment sites) is already published: ProDomino (Nature Methods 2025, public code at github.com/Niopek-Lab/ProDomino, frozen ESM-2 3B plus a 2-layer MLP, AUROC 0.84, 78% lab-confirmed). Its approach is the one to copy, not to rebuild.
- Q3 (linker sequence) is mostly a lookup: cleavable is a lookup against published protease consensus, length is a repeat count. Template retrieval beats generative. The four tick-boxes are not orthogonal.
- Cleavable plus rigid is not one segment. Cleavage sites sit in flexible loops and need 10+ residues of local unfolding. It must be a rigid spacer plus a short flexible protease cassette.
- The defensible design claim is the dataset: Li et al. state that public linker datasets are "virtually nonexistent". Natural inter-domain linkers number in the thousands (LinkerDB), not the low hundreds of engineered ones. LinkerDB is unreachable from this environment (see section 8).
- Cheapest decisive experiments are in `CHECKLIST.md`.

Dorya wants multi-point tethering: several linkers between the same two domains at once, to rigidify relative orientation. This is not junction count. Chain topology is the binding constraint: one polypeptide gives one connection for free, a second tether means leaving and re-entering domain A, and a third repeats it.

---

## 8. Environment facts

- LinkerDB (George and Heringa 2002, `https://ibi.vu.nl/programs/linkerdbwww`) redirects `ibi.vu.nl` to `www.ibi.vu.nl` and then to `zeus.few.vu.nl`. The last is refused with a 403 CONNECT failure. The allowlist must include `www.ibi.vu.nl` (or the `*.vu.nl` suffix). `ibi.vu.nl` alone is not enough. LinkerDB has not been updated since 2002, so its "100x more data" figure may not hold.
- Only WebSearch and WebFetch reach outside the allowlist from the cloud container. PMC serves a captcha. `www.ebi.ac.uk` is denied.
- Web archives and public DNS-over-HTTPS resolvers are denied. There is no mirror route.
- Running thread sessions keep their original sandbox. Only new sessions pick up a changed environment.

---

## 9. Open questions for Dorya

1. **End-to-end (tandem) fusion, internal insertion, or both?** Recommended: end-to-end. Scoping cannot finalize Q2, and the calculator models insertion only as a two-tether approximation, until this is answered. A decision card is open in the project chat.
2. Allowlist for LinkerDB (section 8).
3. Whether the frozen score artifacts under `artifacts/` can be provided, so the mask-policy counterfactual can be reproduced from a clean clone.

---

## 10. Honest limits of this document

- Everything in section 2 was read from committed files. The PR #1 status file was read from that branch. The `SESSION_HANDOFF_2026-10-17.md` on PR #5 was read from the PR body only.
- Not verified in this document: the DLD structural linker count (about 2,600 versus DisProt's about 250), the 370 Li et al. measured linkers, the modality audit claim in Paper I, the DP94 label semantics, and everything in the "not verified" column of section 3.
- Not run: any torch code, the ensemble CLI, any GPU job, any check of the full test suite.
