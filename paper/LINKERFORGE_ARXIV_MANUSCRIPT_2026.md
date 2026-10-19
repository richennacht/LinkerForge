# LinkerForge: An Auditable Prototype for Linker-Candidate Mining and Geometry-Guided Fusion Design

**Status (19 October 2026): Corrected research draft; not submission-ready.** This revision supersedes the earlier manuscript and the previously generated arXiv archives. No independent high-precision mining, published-method superiority, or experimentally functional fusion has been demonstrated. See [the companion audit](../reports/PAPER_CLAIMS_AUDIT_2026-10-19.md).

## Abstract

We describe a computational prototype that combines a sequence-based disordered flexible linker (DFL) classifier, a domain-gap candidate filter, a provenance-stratified retrieval index, deterministic geometry checks, an attachment-site heuristic, and a pairwise ranker for fusion-linker exploration. A frozen three-layer CNN was trained on a family-disjoint B1 V2 split. Its reported three-seed development residue average precision was 0.7611 ± 0.0110, but this is not a proteome-mining precision estimate. A local, exact-sequence-matched InterPro/Pfam audit found just one eligible gap among 28 development proteins, covering 38 of 547 known positive residues and zero known negative residues. A separate repository builder reported 93 true-positive and one false-positive *known residues* after a different domain-annotation path and threshold search. These results cannot establish the precision of called spans in unlabelled proteomes, and their discrepancy requires reconciliation.

The builder's summary reports 27,595 model-scored candidate spans and 171 additional natural/template records across several source corpora. The underlying large dataset and checkpoint are not distributed in this Git branch; only source code and summary reports are tracked. Model-scored gaps are hypotheses, not verified biological linkers or training ground truth. Small retrieval, attachment-site and synthetic-corruption ranking evaluations are exploratory diagnostics. The project does not yet have an independently adjudicated eligible no-linker test, a calibrated span-level error bound, an untouched common-target comparison with published methods, or experimental fusion measurements.

## 1. Scope and evidence categories

A domain gap is a sequence interval between two annotated domains. It is not automatically disordered or flexible. DisProt DFL annotations, structural domain boundaries, engineered templates, and model-inferred gaps describe different evidence types and remain separate. Unannotated residues are unknown, not negative. The active detector and benchmark-firewall rules are in [RESEARCH_HANDOFF.md](../RESEARCH_HANDOFF.md) and the [untouched-comparison contract](../literature/DFL_UNTOUCHED_COMPARISON_CONTRACT_2026-10-17.md).

The present work implements prototype modules for candidate generation, geometric feasibility, retrieval, and ranking. It does not validate fusion expression, folding, function, cleavage, or biological compatibility. A geometric pass means only that the implemented necessary-rule checks did not reject a candidate; it does not prove a physically attainable conformation.

## 2. Methods actually implemented

The detector is a three-layer 1D CNN using frozen protein-sequence embeddings and masked training labels. The domain filter admits 5–60-residue gaps between adjacent, non-overlapping, distinct domain annotations when the sequence and coordinates pass its checks. The implementation in [domain_gate.py](../src/linkerforge/detector/domain_gate.py) does not itself verify that the flanking domains are folded or that the gap is a DFL.

The threshold routine in [conformal_mining.py](../src/linkerforge/detector/conformal_mining.py) searches residue-score cutoffs and computes an ordinary one-sided Clopper–Pearson bound for each cutoff. It does **not** implement a valid, prespecified split-conformal or learn-then-test guarantee: threshold selection reuses the same labels without multiplicity correction, residues within proteins are correlated, and the mining caller classifies whole spans using a different mean/max-score rule. Accordingly, its bound is a descriptive calculation for the selected development data, not a 95% guarantee for new proteins or mined spans.

The geometry module calculates necessary length and accessibility heuristics from specified inputs. Retrieval evaluates natural-linker reconstruction against fixed and other baseline arms. The attachment-site evaluation labels domain termini as positive and interior positions as negative, while its features include terminal position; therefore near-perfect scores on that constructed task cannot establish general attachment-site prediction. The pairwise ranker is tested primarily against programmatically generated corruptions; its scores do not estimate experimental fusion success.

## 3. Observed results and their limits

| Diagnostic | Recorded result | Supported interpretation |
| --- | --- | --- |
| B1 V2 development CNN, three seeds | Residue APS 0.7611 ± 0.0110; AUROC 0.8583 ± 0.0080 | Development discrimination on known residues only |
| Builder's gated development cutoff search | 93 TP, 1 FP among called *known residues*; threshold 0.639120; nominal Clopper–Pearson calculation 4.95% | Descriptive, adaptively selected residue result; not span precision or a prospective error bound |
| Independent exact-sequence InterPro audit | 1 eligible gap / 28 proteins; 38 known-positive, 0 known-negative, 3 unknown residues inside it | Severe coverage loss with this stricter annotation path; precision not estimable |
| Builder's 40-protein control panel | 39/40 proteins had no high-tier call | Mixed auxiliary/annotation-absence proxy; not an independently adjudicated no-linker panel |
| Builder's candidate inventory | 6,556 high-tier plus 21,039 provisional model-scored gaps; 91 CATH-derived, 57 DisProt-curated, 23 templates | Reported inventory, not 27,766 verified DFLs; the data files are not tracked in Git |
| X2 retrieval, CATH cohort | Hybrid arm: 43/91 within ±2 residues; geometry-only: 31/91 | Small retrospective natural-linker reconstruction, not engineered-fusion performance |
| Attachment-site constructed task | Reported AUROC 1.0000 and 289/289 top-1 on natural cohort | Terminal-vs-interior task with an explicit terminal feature; not an independent biological site test |
| Pairwise ranker | Reported 147/148 top-1 against six synthetic corruptions | Distinguishes chosen perturbations; not superiority over published predictors or experimental designs |

The builder used a 0.10 target for its *residue-level* nominal FDR calculation; the 0.639120 cutoff is not a 0.92285 span threshold, and the paper must not interchange the two. The builder's high-tier span calling rule uses a frozen high residue threshold of 0.735198 and a provisional mean threshold of 0.50; it is not the calibrated rule evaluated in the table. The earlier manuscript's 44.74%-to-1.06% **proteome** FDR comparison, 95% guaranteed precision, and 97.50% adjudicated no-linker specificity were unsupported.

The reported 184,554 unreviewed reference-proteome domain-gap observations are a separate stage-one audit and were not all CNN-scored. Adding them to 161,317 scanned proteins as though they were the same kind of record is invalid. Source corpora may also overlap, so these counts are not a unique-protein total without deduplication evidence.

## 4. Reproducibility and release status

Tracked source includes [the dataset builder](../scripts/build_mined_and_natural_linker_dataset.py), [the summary](../reports/mined_and_natural_linkers_v1_summary_2026-10-19.json), [the X2 diagnostic](../reports/x2_retrieval_and_part2_design_2026-10-19.json), [the attachment diagnostic](../reports/part3_attachment_scorer_2026-10-19.json), and [the ranking diagnostic](../reports/part4_end_to_end_ranker_and_sfgfp_mcherry_2026-10-19.json). The claimed large dataset files and frozen CNN weights are ignored local artifacts, not included in this Git branch. The hashes in reports are provenance assertions that cannot be independently checked from this branch alone. The existing PDB exports are geometric illustration artifacts, not AlphaFold, Rosetta, MD, or experimental structures of the proposed fusions.

Prior to a submission or mining release, the team must (1) reconcile the two domain-annotation routes and rerun the exact production caller on a sealed span-level set, (2) collect independent adjudication of eligible no-linker and positive gaps, (3) preregister a cluster-aware calibration procedure with multiplicity handling and evaluate it once on an untouched family-disjoint set, (4) publish or reproducibly regenerate source data and model weights subject to licensing, and (5) compare against matched published methods on identical targets and coverage. CAID3 has already been used retrospectively and cannot serve as a new blind result. No claim here requires that future tests succeed.

## 5. Conclusion

LinkerForge now has an integrated *prototype* for proposing and ranking linker candidates. Its current numbers support software and small retrospective diagnostic claims only. High-precision proteome mining, biological attachment-site generalization, and functional fusion-linker improvement remain open hypotheses.
