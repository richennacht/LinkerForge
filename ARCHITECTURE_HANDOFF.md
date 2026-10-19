# LinkerForge architecture handoff

Date: 17 October 2026. Written for someone picking this up cold from the repository alone.

Read order: this file, then `SPEC.md` and `CHECKLIST.md` (on PR #6), then `RESEARCH_HANDOFF.md` on the active branch. Those three are the binding rules. This file is the design map: what to build, how each part is shaped, which papers it rests on, what already exists, and what is not verified.

Status update (18 October 2026): Dorya explicitly resumed work. PRs #1–#7 are merged into `codex/dfl-research-handoff`. This is a 17 October architecture snapshot; read `CURRENT_DIRECTION_2026-10-18.md` for the current classifier-first implementation order. Its evidence and benchmark limitations still apply.

Every draft PR named here targets `codex/dfl-research-handoff`. As of this writing the CI workflow does not run on PRs into that branch (see section 11), so a missing red check is not a pass.

---

## 0. Two strands, one dependency chain

The project has two strands. They are connected by one dependency: the classifier mines the natural linkers that the design model is trained on.

```
Strand A: detector (active)
  classifier (component 1) ──> mining pipeline (component 2) ──> mined natural linkers
                                                                       │
Strand B: design model                                                 │
  calculator (component 3) ─────────────────────────────┐              │
  attachment-site scorer (component 4) <── labels from flanking residues┘
  linker sequence (component 5) <── mined linkers, template library
```

Component 3 (calculator) is the only piece that is already built. Components 1 and 2 are the most gated, because they need a family-disjoint split and a clean negative set before any claim can be made.

Dorya's stated intent, in their words relayed through the scoping thread: build a classifier with enough trustworthiness to run, use it on unannotated proteins to deduce whether a linker is present, then train a model that takes two unrelated domains and returns (a) whether they are compatible for linking, (b) the most suitable attachment sites, and (c) the linker type, learned from the mined natural linkers.

---

## 1. Verified state, in one place

Numbers below are from committed files or from named sources. Anything not re-checked is marked.

- Matched SAM-like baseline, three seeds on protected development data: mean APS **0.7918**, AUROC 0.9076, Fmax 0.7022 (7,427 known residues, 1,267 positive). Source: `reports/DFL_SAM_REPEAT_AND_LOCKED_CAID3_2026-10-17.md`.
- One-time locked CAID3 common-27 comparison, **spent**: LinkerForge baseline APS **0.241**, AUROC 0.695, Fmax 0.345. IPA-AF2-Linker is best at APS 0.482, Fmax 0.491. Source: SPEC.md section 2.2.
- Mask-policy result: with model scores frozen, APS moves from **0.450** (unknown residues excluded) to **0.115** (unknown residues counted as negative). 93.35% of labels are unknown. Source: `reports/DFL_POOR_TRANSFER_CAUSE_AUDIT_2026-10-12.md` lines 20 and 21, `reports/DFL_LABEL_MASK_AUDIT_2026-10-12.md`. **The frozen score files are gitignored under `artifacts/` and a clean clone cannot recompute this.** H3 and P0 in CHECKLIST.md cover that gap.
- No family-disjoint split exists. Splits are by protein ID. The current split `caid3_protected_dfl_v1` is not eligible for new training because it overlaps the DP94 candidate benchmark. Preflight fails closed.
- No conformal prediction exists in the repository.
- No adjudicated no-linker set exists, so no specificity claim and no mining claim can be made yet.
- Architecture is not the bottleneck. The best linker model in the repo is a 3-layer CNN (LINKER-Pred family). SAM and Mask2Former decoders were negative in-repo. Data is the lever.

---

## 2. Corrections log

Each row gives the original claim, the correction, and which one is right. The corrections are the most useful part of this record. The original text is kept so nobody re-introduces it from an old thread.

| # | Original claim (where it appeared) | Correction | Which is right, and source |
|---|---|---|---|
| C1 | Best published CAID3 linker predictor reaches "F1 0.38" (high-precision survey; relayed into scoping and the claims thread) | Wrong metric and wrong number. The field's best on common-27 is APS 0.482 and Fmax 0.491. | APS 0.482 / Fmax 0.491, IPA-AF2-Linker, SPEC.md section 2.2. The 0.386 in that table is LINKER-Pred2's **APS**, a plausible source of the mix-up, not confirmed. |
| C2 | Chain extension per residue is 3.0 to 3.5 Å (both architecture surveys; relayed to Dorya as fact) | 3.0 to 3.5 Å is a **design convention**. Used as a bound it makes real segments look impossible. | Bound is **3.8 Å per residue**, the virtual Cα–Cα distance. Verified over 412 residues of real structure in PR #3: no segment exceeds 3.8, and short extended segments routinely exceed 3.5. Helical rigid linkers: about 1.5 Å axial rise per residue (both surveys agree). |
| C3 | "How many linkers" means how many junctions a topology has (both architecture surveys) | Misread. Dorya means several linkers between the same two domains at once, at different contact points, to fix relative orientation. | Multi-point tethering, constraint satisfaction over domain B's placement. Built in PR #3. Chain topology is the binding constraint. |
| C4 | Compatibility is a binary classifier (the original framing) | Any two domains are linkable given a long enough linker. The question is conditional feasibility. | Conditional feasibility. Skip an end-to-end compatible/incompatible classifier: no labels, and it hides the failure modes. |
| C5 | Attachment count is a regression output | Counts are discrete topology choices. | One linker for tandem fusion, two for insertional (entry and exit), one per junction beyond that. Count falls out of the chosen topology. |
| C6 | Dorya's target of 90% recall at 99% precision (original brief) | Not reachable. The ceiling is the labels: unannotated residues count as negatives, so many apparent false positives are real linkers. | Keep reporting real numbers and state where the target sits. Changing experiment parameters to reach a number is not an option (SPEC section 6). |
| C7 | Template library is a contribution (library thread, before the competitor read) | Not a contribution. SynLinker (2015) shipped 2,260 curated linkers with length, composition, solvent-accessibility and proteolytic-site filters, plus fusion structure prediction from two domains. | SynLinker. Build the library to use, not to publish (PR #2, SPEC section 5 item 4). |
| C8 | The 2026 Li et al. paper is the competitor for linker design (early scoping) | It is a **scorer, not a generator**. One domain pair, one length, no property switches, no protease sites, no generation step. | The real limit on the design claim is SynLinker plus the dataset gap. Li et al.'s figures are unverified (see section 13). |
| C9 | Cleavable plus rigid is a single linker class (template library thread, untested prediction) | The prediction was that protease sites break helices. Now cited as confirmed: cleavage sites sit in flexible loops and rarely in or near helices, needing a local unfolding step of at least 10 residues. | Cleavable plus rigid must be a **rigid spacer plus a short flexible protease cassette**. The Hubbard proteolysis analysis is cited by the Q3 survey and was not re-read here. |
| C10 | The four property tick-boxes (cleavable/non-cleavable, rigid/flexible) are orthogonal (Dorya's framing) | Not fully orthogonal. Cleavable plus rigid is the awkward corner. | Treat the four combinations as dependent. The cleavable+rigid case uses the cassette design in C9. |
| C11 | Template retrieval beats generative (library thread prior) | Literature supports ranking retrieval first, but it is a **hypothesis until tested**. A natural linker count in the thousands (LinkerDB) is larger than the engineered low hundreds the prior assumed, which weakens the memorization argument somewhat. | Test with experiment X2 (section 12) before building a generator. |
| C12 | Low pLDDT is a quality signal to filter generated linkers on | Wrong for linkers. Low pLDDT is the disorder signal (CAID uses 1 minus pLDDT). Filtering on confidence rejects correctly disordered linkers. | No pLDDT filter on linkers. Verify only that a rigid linker is helical and that domains keep their folds. |
| C13 | Structure features help: AlphaFold pLDDT and RSA are worth adding (high-precision survey) | Conflict: IPA-AF2-Linker uses them, but IDPFunNet found them worse than ProtT5 alone. | Unresolved. Structure features must be settled by the project's own ablation, not by literature (SPEC section 4 and checklist X1). |
| C14 | Linker sequence size: about 2,600 DLD linkers vs about 250 DisProt (Strand A survey) | Unverified in this repository. | Treat as unverified (SPEC section 10). |
| C15 | Every Li et al. 2026 figure, including the structure-vs-sequence ablation | Unverified. PMC is captcha-walled, ScienceDirect blocks robots, Dorya has no library access. Only the author version up to part way through Results was readable. | Permanently unverified until a library copy is available. Do not use as benchmark or comparison. |
| C16 | Multi-tether connection is novel | Two connection points are established practice (GPCR crystallography fusion literature; Baker-lab two-rigid-junction work, 2021). No precedent found for three or more simultaneous peptide linkers. | Negative search result only. The crosslink and staple literature was **not searched**. Do not write "first" near multi-tether work until that search is done. |

General lesson recorded in the claims register: a number that is conventional in a design context but physical in form does not say which it is, and nobody downstream re-derives it. Re-derive before quoting (C2 is the example).

---

## 3. Component 1: high-precision classifier (Strand A)

**Job.** Given a protein sequence, label each residue as flexible linker (DFL), not linker, or abstain. Operate at a precision high enough that its positive calls can be mined at scale.

**Architecture (proposed, from the high-precision and linker-predictor surveys).**

- Input: sequence. A frozen protein language model encoder (ESM-2 is the repo default; ProtT5/ProstT5 is an alternative the surveys suggest combining). Optional AlphaFold features (pLDDT, relative solvent accessibility, a "between two confident domains" signal from PAE or a domain parser) only if the project's own ablation shows they help.
- Head: small 1D convolutional head (the 3-layer CNN family is the best linker model in the literature the surveys reviewed). No SAM or Mask2Former query decoder: no prompt at inference, and the repo's own negative results rule it out.
- Span decoding: a constraint that allows a linker only between two ordered regions (the CRF-style trick SignalP and DeepTMHMM use). Cuts false calls on terminal tails. Code exists in PR #1 (span decoders).
- Learning under unknown labels: positive, unknown and negative labels are kept separate. Unknown is not negative (SPEC section 5 item 10). PU loss (nnPU) is in PR #1, but its `class_prior` is null, so the launcher refuses to run until the prior is set from a source that is not CAID3.
- Abstention: three-way output (linker, not linker, abstain). Abstention is a feature here, not a shortfall.

**Calibration, the part that carries the claim.**

- Interim (already in PR #1): Clopper-Pearson precision threshold, `detector/selective.py`. It is residue-level and ignores within-protein correlation, so it is optimistic. It is selective prediction with a binomial bound, **not** conformal risk control. It is an input to the next step, not a substitute for it.
- Target: conformal risk control via Learn then Test (Angelopoulos et al., 2021; citation from memory, not checked here) on a small clean calibration set split by protein family. Only this route certifies a precision level.
- Hard rule: a threshold chosen on the split it is reported on is **not a guarantee** and must never be described as one, especially in a title (SPEC section 5 item 1). The current 0.90 operating point is of this kind.
- Exchangeability caveat: the guarantee holds only within the family-disjoint design. Proteins sharing family structure break the assumption, so B1 comes first.

**Metrics to report.** APS, AUROC, Fmax on equal-coverage sets; **precision at the coverage a large mining run would actually use**, with the abstention rate. Headline APS is not the decision metric (Dorya's instruction in the relay: precision over recall).

**Data.**

- DisProt-derived flexible linker spans are the positive source. DLD structural independent and dependent domain linkers are a separate, secondary task and must not silently become DFL ground truth.
- Benchmark sequences and IDs may be used for exclusion. Benchmark labels, per-protein errors, scores and thresholds may not be used for training, filtering, feature selection or model selection.
- DLD's structure-derived set is the proposed route to more training data (structural independent and dependent domain linkers, about 2,600 per the survey, unverified here). Train it as a separate masked head.

**Gates before any training (from RESEARCH_HANDOFF.md rule 3 and CHECKLIST B1 to B3).**

1. Family-disjoint train, development and auxiliary split with a hash-checked preflight against all nominated benchmark sequence registries, including `data/processed/prospective_dfl_holdout_dp94_v1/`. Current split fails. Proposed clustering: MMseqs2 at 30% identity and 50% coverage (from MINING_RUN_SIZING.md, a proposal, not tested).
2. Adjudicated no-linker set, disjoint from B1 and CAID3. Without it there is no specificity claim, so no mining claim.
3. Label-ceiling analysis: estimate how many apparent false positives are unannotated real linkers.
4. One local GPU job at a time. At most three architecture variants and six GPU-hours per round. Frozen encoders: extract embeddings once on a free cloud GPU, train small heads locally on the 6 GB RTX 3050.

**Existing code.** PR #1 (`claude/dfl-high-precision-toolkit`): nnPU loss, ensemble scoring, Clopper-Pearson threshold, AlphaFold feature loader, span decoders. 43 CPU tests pass, 2 skipped (numpy only). **Untested:** torch loss and gradients, PU launcher wiring (py_compile only), ensemble CLI (never run), any run on real scores or real AlphaFold data, and the full repo test suite.

**Not built.** Family-disjoint split (B1), conformal thresholds (B2), adjudicated no-linker set (B3), label-ceiling analysis (B4), any PU training run.

---

## 4. Component 2: mining pipeline over UniProt

**Job.** Run the classifier across unannotated proteins and output spans with confidence, flanking residues and provenance. The mined set trains components 4 and 5.

**Architecture (proposed).**

- Stream sequences from UniProt. Do not store embeddings for everything: Swiss-Prot embeddings alone would take about 535 GB (PR #4 arithmetic). Stream, score, write only the calls.
- Order: Swiss-Prot first (575,748 sequences, 209M residues, UniProt 2026_03). Then UniRef50 (38.8M clusters) on the college cloud, which the handoff names as the place for long jobs. Full TrEMBL (149.4M) is not feasible on the stated budget.
- Abstain-by-default: only the high-precision operating point writes mined spans. The threshold is taken from component 1's calibration set, never from the mined set.
- Leakage control: exclude any protein homologous to a benchmark sequence (benchmark IDs and sequences may be used for exclusion). Proposed MMseqs2 30%/50%, untested.
- Output record per span: source accession, UniProt release, residues, classifier score, threshold ID, abstention flag, flanking 10 residues, config hash, commit SHA. Label natural observations distinctly from predictions (README rule).

**Existing work.** PR #4 (`claude/uniprot-mining-sizing`, draft, `MINING_RUN_SIZING.md`): counts, storage arithmetic, Swiss-Prot first, leakage design, resumable pipeline design. **Arithmetic only.** Not verified: current AlphaFold DB coverage of UniProt (the container could not reach UniProt or EBI), AFDB entry count, encoder throughput, TrEMBL and UniRef50 mean lengths.

**Blockers.** No family split (component 1 gate). Environment access to UniProt and AFDB from the cloud container was denied during this work. Verify when the environment allows it.

**Not built.** The pipeline itself, AFDB coverage check, encoder throughput measurement.

---

## 5. Component 3: geometric feasibility calculator

**Job.** Given two domain structures and candidate attachment points, return the linker length each site pair requires, and whether a set of tethers is simultaneously satisfiable. Deterministic. Explainable. The calculator is the only part of the design model that is solid today.

**Architecture (built, PR #3).** Package `linkerforge.geometry`, standard library only, no torch, runs on the laptop. Modules:

- **Span window.** Required linker length per attachment pair, with both ends: a physical lower bound (the chain must reach) and a design figure, plus an upper bound for helical linkers (a rigid rod that is too long holds the domains apart at a wrong fixed angle). Feasibility is a window, not a ceiling.
- **Extension constant.** 3.8 Å per residue as the bound (virtual Cα–Cα). 3.0 to 3.5 Å is a design convention, not used as a bound (C2). Helical rise about 1.5 Å per residue.
- **Accessibility.** Shrake–Rupley solvent accessible surface at each anchor, checked analytically in tests. Radii and maximum-ASA values cited in PR #3 (van Rosmalen, Tien, Bondi, Sowdhamini). Not re-checked here.
- **Occlusion and detour.** Clearance along the path; a detour bound when the path crosses a domain; whether the chain exits toward the partner.
- **Multi-tether satisfiability.** Constraint satisfaction over the placement of domain B, not closed form. Verdicts are two-sided: a violated pairwise distance-geometry bound **proves** infeasibility; a returned placement is a **witness** that can be rechecked; anything else returns `undetermined`. This is the only provable, non-tuned claim in the project (register I1). Its boundary: geometric admissibility is necessary, not sufficient.
- **Tether types.** Peptide linkers and engineered disulfides. A disulfide has a fixed short distance between two cysteine side chains, not a length-dependent range, so it needs its own span constraint.
- **Chain topology report.** Every verdict reports the topology it implies. One polypeptide gives one connection free. A second tether needs leaving domain B and re-entering domain A (interleaved or insertional). A third repeats it. Alternatives (circularized construct, two chains, chemical crosslinks) are separate engineering problems.
- **Calibration.** Extension bound and clash tolerance calibrated on real structures: 4E43 (dimer) and 4AKE (open) fixtures, with provenance in the repo.

**Tests.** 80 pass. Real-structure calibration included. Run `python3 -m unittest` on the package tests (exact invocation is in PR #3).

**Not covered by geometry, so not claimed.** Whether a construct expresses or folds. Whether two domains are compatible as folds rather than shapes. Which loops tolerate interruption. Learned residual risk (steric clash beyond geometry, interface formation, aggregation) has no label source and is not built.

**Open decision (blocks the insertion case).** End-to-end (tandem) or internal insertion? The calculator models insertion only as a two-tether approximation. Recommended default: end-to-end. Decision card is open in the project chat. Not answered.

**Novelty of multi-tether.** Unverified (C16). Search the crosslink and staple literature before any novelty word.

**Papers this rests on.** ALBATROSS (end-to-end distance from sequence, R² 0.986, 30 to 60 sequences per second on CPU): https://www.nature.com/articles/s41592-023-02159-5. SynLinker, which already filters 2,260 linkers by end-to-end distance: https://academic.oup.com/bioinformatics/article/31/22/3700/241706 (doi 10.1093/bioinformatics/btv447). GPCR fusion crystallography (two connection points, established practice). Baker-lab two-rigid-junction work (2021, from the scoping thread). Constants and radii as cited in PR #3.

**Next steps (not started).** Insertion-specific loop geometry (blocked on E2). JSON schema and docs are present in PR #3; read them before extending.

---

## 6. Component 4: attachment-site scorer

**Job.** For each domain, score each residue on whether it tolerates a linker or insertion. Rank candidate site pairs, combining the scores with the calculator (component 3).

**Architecture (proposed, consensus of two independent surveys).**

- **Copy ProDomino.** Frozen ESM-2 3B encoder plus a two-layer MLP per residue. Trained on 174,872 natural domain-interruption cases from InterPro and UniRef50. AUROC 0.84 on held-out experimental data; about 78% of predictions confirmed in the lab. Code is public: https://github.com/Niopek-Lab/ProDomino. Paper: https://www.nature.com/articles/s41592-025-02741-z (Nature Methods 2025).
- **Two-tower.** Run the scorer independently on each domain. Attachment tolerance is a property of the domain, so the learned part factorizes over single domains. It trains on single-protein data, not fusion-pair data (which does not exist).
- **Pair ranking.** Closed-form geometric re-ranker over site pairs, using component 3's output.
- **Keep the tick-boxes out of the encoder.** Rigid, flexible and cleavable are constraints on the ranker and on component 5, not inputs to the tolerance model. One trained model then serves all combinations.
- **Fine-tuning set.** Domain insertion profiling of Kir2.1: 759 donor motifs at all 435 positions, over 300,000 variants, fitness scored against a flexible linker control. Local backbone flexibility predicts tolerance, stiffness predicts failure. Paper: https://www.nature.com/articles/s41467-021-27342-0. Too narrow to pretrain on; use to fine-tune.

**Labels.**

- Natural insertion events (ProDomino's source).
- Flanking residues of natural linkers mined by component 2. The surveys argue these are positions that tolerate a linker, since they are by definition where a linker sits. **This is an untested assumption.** It must be checked before it trains anything: a flanking residue may be a domain boundary residue for reasons unrelated to tolerance.

**Validation.** Hold out a domain family. Report AUROC and precision at a fixed coverage. Compare against ProDomino's own output as a baseline. The Kir2.1 data is the only experimental check the surveys identified.

**Caveats.** The surveys could not open the Kir2.1 paper's full text or the Gly-Ser flexibility paper in full; those citations come from search metadata. Verify before relying on them.

**Not built.** Everything. This is the first component to build after the classifier's gates pass, because it can run on single-protein data.

---

## 7. Component 5: linker sequence

**Job.** Given two domains, chosen attachment sites, a required length, and the cleavable and rigid/flexible choices, return a linker sequence.

**Decomposition (the most useful finding in the Q3 survey).** Only about a third of the task is a learning problem.

- **Cleavable.** A lookup against published protease consensus sequences. Generating a *novel* protease site means predicting protease specificity, a harder and different problem. Not learned.
- **Length.** A repeat count. Not learned.
- **Rigid versus flexible in context.** The only part where the specific domain pair changes the answer. This is where a model earns its place: context fit, and escaping repetitive (GGGGS)n, which recombines during cloning.

**Architecture ranking (from the Q3 survey; the ranking is the survey's view, not tested).**

1. **Template library plus substitution.** Ships first and becomes the baseline. Built in PR #2 (30 entries, each with a literature source, 28 tests, two-axis property caller, coverage audit). **Not a contribution** (C7). Build to use.
2. **Fill-in-middle protein language model** fine-tuned on mined linkers. Domain A as prefix, domain B as suffix, linker as the middle. ProtFIM is the closest published match: it beat ProGen2-large at about a thirtieth of the parameters. Paper name only, not read here.
3. **Adapters on a frozen language model** if per-combination label counts are small.
4. **RFdiffusion inpainting** for rigid linkers only.
5. **De novo sequence diffusion.** Not justified by the surveys.

**The four combinations are not orthogonal.** Cleavable plus rigid cannot be one segment. Cleavage sites sit in flexible loops and need local unfolding of 10 or more residues. The design is a **rigid spacer plus a short flexible protease cassette** (C9, C10).

**Validation, and what cannot be validated.**

- **Experiment X2 (the decisive, cheapest test, not yet run).** Hold out natural linkers. For each, retrieve the closest template from the library using flanking domain context only. Measure sequence recovery and length error. Include a retrieval baseline (nearest library entry) and a length-matched frequency baseline. A strong retrieval result means a generator has little room. Report either way.
- **In-silico checks.** A rigid linker is helical (secondary structure check). Both domains keep their folds. A flexible linker has no single reference conformation, so self-consistency RMSD against an intended backbone does not apply.
- **Not checkable without a bench.** That the fusion works, expresses or folds.
- **Do not filter on pLDDT** (C12).

**Data.** Natural inter-domain linkers number in the thousands per LinkerDB (George and Heringa 2002, doi 10.1093/protein/15.11.871; https://ibi.vu.nl/programs/linkerdbwww). LinkerDB was last updated in 2002, so its size figures are a snapshot. **Unreachable from this environment** (E1). Engineered fusion linkers number in the low hundreds.

**Competitor check.** Li et al. 2026 (doi 10.1016/j.synbio.2026.02.003, Synth Syst Biotechnol 13:362-375) is a scorer, not a generator, and does not take this niche. Its figures are unverified (C15). SynLinker (2015) covers curation (C7).

**Defensible claim, if any.** The dataset, not the generator. Li et al. state in print that public linker datasets are "virtually nonexistent". Mined natural linkers are the gap. A generated linker can be claimed as novel and plausible. It cannot be claimed as better than a SynLinker retrieval, since mined linkers have no functional labels (Li et al. have 370 measured ones).

**Not built.** Fill-in-middle model, adapters, validation harness, X2. The library is partly built (PR #2).

---

## 8. Cross-cutting: viability flags

No component may report a probability that a fusion works.

- Best zero-shot predictor of design success: ESMFold pLDDT, AUC 0.73 across 614 designs. All metrics combined reached 0.74.
- Solubility predictors top out around MCC 0.46 on natural *E. coli* proteins, not chimeras.
- No published dataset maps fusion constructs to expression outcomes.

Present any viability output as a **flag with its AUC stated**, never as a success probability. Geometry, the calculator and the library can say what is physically or compositionally admissible. They cannot say the construct folds or expresses.

---

## 9. Integrity constraints (restated from SPEC.md section 5)

1. The precision-0.90 operating point is a threshold selected on the split it is reported on. It is not a guarantee. Never use "guarantee" or "certified" near it.
2. 0.7918 and 0.241 are different quantities. Never put them in the same table, sentence or abstract without a clear label and a reason.
3. CAID3 is spent. Anything selected, tuned, thresholded or filtered using CAID3 is not a result, including choosing the common-27 subset after seeing it.
4. The template library is not a contribution (C7).
5. Li et al. 2026 figures are permanently unverified. Do not use them as a benchmark or comparison, including the structure-vs-sequence ablation.
6. Do not use "first" near multi-tether work until the crosslink and staple literature is searched (C16).
7. pLDDT is not a validity score for linkers (C12).
8. Never report a success probability for a fusion (section 8).
9. Breadth is not quality. Library size and mined-set size are not quality claims.
10. Unknown is not negative in any new label, metric or claim.
11. Do not manufacture records from advertised database counts. SynLinker's row-level data was not acquired (the publisher blocks automated access), so no SynLinker records exist in this repository.

---

## 10. Literature-survey citations kept with their findings

- The borrowable-architectures survey (Strand A): the best linker predictor is a small CNN ensemble on a protein language model. Architecture gains came from data, not structure.
- The linker predictor survey: DLD structural linkers are downloadable. Negatives can be handled with positive, unknown and negative labels (IPA-AF2-Linker). CAID3 target numbers must be checked against the same subset.
- The high-precision survey: conformal risk control and Learn then Test are the only routes to a guaranteed precision. Ranked: clean calibration set split by protein, Learn then Test, nnPU plus a five-model ensemble, cleanlab for label errors, three-way output with abstention.
- The domain-compatibility survey (Q1): three layers. Geometry is computed. Attachment tolerance is learned (ProDomino). Viability is a flag only.
- The attachment-site survey (Q2): two-tower scorer plus closed-form geometric re-ranker. Nothing published takes two domains plus a linker specification and returns attachment points. ProDomino is the nearest.
- The linker generation survey (Q3): template library first, then fill-in-middle. Cleavable plus rigid is confirmed awkward.
- The design-model scoping: Q1 is a calculator, Q2 is already published, a third of Q3 is a lookup. The project's actual contribution in the design half is narrow.

The surveys' own threads are the sources. Their artifact pages are linked in the claims register and publication plan (artifacts `DVuSZKvmddaZUyY9bvWY8z` and `NwQKGpHyCA1RseRQmRYzhh`, and the competitor page `12ac39dd-3288-4e55-8078-9e9e571fe0ae`). Those pages are not in git; the content they carry is restated above.

---

## 11. Branch and PR state

Base branch for all of these: `codex/dfl-research-handoff`. All draft. Dorya is not a repository collaborator, so review requests fail. The repository is the system of record.

| PR | Head branch | Contents | Tested | Not verified |
|---|---|---|---|---|
| #1 | `claude/dfl-high-precision-toolkit` | nnPU loss, ensemble scoring, Clopper-Pearson threshold, AlphaFold feature loader, span decoders | 43 CPU tests pass, 2 skipped | torch paths, PU launcher, ensemble CLI, real data, full test suite |
| #2 | `claude/linker-property-library` | 30-entry library from Chen et al. 2013 (doi 10.1016/j.addr.2012.09.039), two-axis property caller, coverage audit | 28 tests | Function (literature-reported only); some source inconsistencies flagged per entry; not a contribution (C7) |
| #3 | `claude/linker-geometry-feasibility` | `linkerforge.geometry`: span window, accessibility, occlusion, multi-tether satisfiability, topology | 80 tests, stdlib only | Insertion geometry; expression or folding; novelty of 3+ tethers |
| #4 | `claude/uniprot-mining-sizing` | `MINING_RUN_SIZING.md` | Arithmetic only | AFDB coverage, AFDB count, encoder throughput, mean lengths |
| #5 | `claude/session-handoff-claims-register` | Session handoff: claims register and publication plan | Documentation only | Summarizes the reports; does not add evidence |
| #6 | `claude/spec-and-checklist` | `SPEC.md` and `CHECKLIST.md` | Documentation only | Binding rules; read first |
| this | `claude/architecture-handoff` | This file | Documentation only | Every architectural claim not already in SPEC, the surveys or PRs above |

**CI gap.** `.github/workflows/ci.yml` triggers `pull_request` only for `main` and `integration`. PRs into `codex/dfl-research-handoff` run no CI. No Python job exists, so `tests/` never runs in CI. Verify Python changes locally and state the command in the PR. Fix in its own PR (CHECKLIST H1).

---

## 12. Open items, blockers and next actions

**Decisions only Dorya can make.**

1. End-to-end (tandem) fusion or internal insertion, or both? Recommended: end-to-end. Blocks the insertion geometry and the final Q2 scope.
2. Whether the frozen score files under `artifacts/` can be provided, so the mask-policy result can be reproduced from a clean clone (P0, H3).

**Environment blockers.**

3. **LinkerDB (E1).** The redirect chain is `ibi.vu.nl` to `www.ibi.vu.nl` to `zeus.few.vu.nl`. The last is refused with a 403 CONNECT failure. `ibi.vu.nl` alone is not enough, and the first hop was reachable in one session. The allowlist needs `www.ibi.vu.nl` or `*.vu.nl`. Changed in Project settings, environment. Only new sessions pick up the change.
4. **UniProt and AlphaFold DB.** Denied from the cloud container during the sizing work. Needed for the coverage check in component 2.
5. **Li et al. 2026 (E3).** Permanently unverified without library access.

**Gated work, in dependency order (from CHECKLIST.md; do not start a blocked item).**

- Paper I (no GPU, no new split, no new data): P0 to P7. Benchmark and methods analysis: equal-coverage CAID3 reanalysis, mask-policy reproduction, unknown-label inventory, modality audit (verify against tool documentation), about ten clean negatives. Target is a methods journal. The abstract must say no new method is introduced.
- Paper II (gated; do not draft until B1 passes): B1 family-disjoint split with passing preflight; B2 conformal thresholds with abstention; B3 adjudicated no-linker set; B4 label-ceiling analysis; B5 one-time untouched test.
- Cheap decisive experiments: X1 (do relative-position structural features transfer across domain pairs; test on natural linkers grouped by domain pair, sequence-only baseline); X2 (template retrieval reconstructs held-out linkers from flanking context, with baselines).

**What the record cannot yet support.** No family-disjoint split. No conformal prediction. No adjudicated no-linker set, so no specificity or mining claim. No design model. The pilot dataset fails its own label gate at 91 records with no resolved linker interiors. The sequence-pattern side study (PR-level survey in `literature/LINKER_MOTIF_DISCOVERY_SURVEY_2026-10-12.md`) is descriptive only and does not assign function.

---

## 13. Verification status

- **Read from committed files in this repository:** SPEC.md and CHECKLIST.md (PR #6), PR #1 status file, RESEARCH_HANDOFF.md, AGENTS.md, the motif survey, PR #5 body, PR #4 and PR #3 summaries as reported in threads.
- **Read from survey threads (sessions that produced them):** the surveys for components 1 to 5, the claims register, the competitor read, the scoping read, the calculator and sizing summaries. Not re-read line by line in this session.
- **Unverified in this document:** DLD structural linker count (about 2,600 vs about 250 DisProt); the 370 Li et al. measured linkers; the modality audit claim; the DP94 label semantics; every figure from Li et al. 2026; the Hubbard proteolysis result (cited, not read); the Kir2.1 and Gly-Ser papers (cited from metadata); the Learn then Test and nnPU citations (from memory); the PR #3 accessibility constants (cited, not re-checked); AFDB coverage; encoder throughput; the flanking-residue label assumption (component 4); IDPFunNet's negative pLDDT result.
- **Not run:** any training, any split, any benchmark evaluation, the full test suite, the torch paths, the ensemble CLI, any GPU job.

---

## 14. References

Papers and tools named above. Links are included where the source thread or repository recorded one. Entries without a link were named in a survey and not located here.

- ProDomino, Nature Methods 2025: https://www.nature.com/articles/s41592-025-02741-z ; code https://github.com/Niopek-Lab/ProDomino
- Domain insertion profiling, Kir2.1: https://www.nature.com/articles/s41467-021-27342-0
- ALBATROSS, Nature Methods 2023: https://www.nature.com/articles/s41592-023-02159-5
- SynLinker, Liu et al. 2015, Bioinformatics: https://academic.oup.com/bioinformatics/article/31/22/3700/241706 (doi 10.1093/bioinformatics/btv447)
- George and Heringa 2002, linker classification: doi 10.1093/protein/15.11.871. LinkerDB: https://ibi.vu.nl/programs/linkerdbwww (unreachable, E1)
- Chen et al. 2013, linker properties: doi 10.1016/j.addr.2012.09.039
- Li et al. 2026, Synth Syst Biotechnol 13:362-375: doi 10.1016/j.synbio.2026.02.003 (unverified, see C15)
- DLD source study 2026: https://www.sciencedirect.com/science/article/pii/S0022283626003670
- DFLpred 2016: https://pmc.ncbi.nlm.nih.gov/articles/PMC4908364/
- ELM resource 2022: https://pmc.ncbi.nlm.nih.gov/articles/PMC8728146/
- STREME, Bailey 2021: https://pmc.ncbi.nlm.nih.gov/articles/PMC8479671/
- DiMotif, Asgari et al. 2019: https://www.nature.com/articles/s41598-019-38746-w
- LinkCraft 2026: https://www.sciencedirect.com/science/article/pii/S0022283626001877
- Prolines in linkers 2026: doi 10.1371/journal.pcbi.1014588
- CAID3 (benchmark, spent), LINKER-Pred, LINKER-Pred2, IPA-AF2-Linker, DisoFLAG-DFL, APOD, DFLpred: reported in SPEC.md section 2.2. Primary links not recorded in this repository.
- Learn then Test, Angelopoulos et al. 2021 (conformal risk control; citation from memory)
- nnPU, Kiryo et al. 2017 (positive-unlabeled learning; citation from memory)
- ProtFIM (fill-in-middle protein language model): named in the Q3 survey, not located here
- Hubbard, proteolysis and local unfolding: named in the Q3 survey, not located here
- ESMFold pLDDT zero-shot design-success AUC 0.73; solubility MCC around 0.46: as cited by the Q1 and Q3 surveys, not located here
- MMseqs2 (family clustering proposal): named in MINING_RUN_SIZING.md
- UniProt release 2026_03 (Swiss-Prot 575,748; TrEMBL 149.4M; UniRef50 38.8M): PR #4
