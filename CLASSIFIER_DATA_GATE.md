# Classifier workstream: first-stage data gate

This document assigns the **first** build stage to the DFL residue classifier.
Geometry, cited templates and the retrieval baseline are downstream, separate
work. `CURRENT_DIRECTION_2026-10-18.md` gives the complete order.

## Implemented, not yet unblocked

- `scripts/predict_dfl_frozen.py` provides a hash-checked, provisional inference
  interface for the historical checkpoint. It abstains by default. It is **not**
  a certified mining classifier.
- `scripts/preflight_dfl_family_training.py` is a tracked CPU gate for a *new*
  family-disjoint split. It requires pinned sequence-only CAID3 and DP94
  registries, an honest `retrospective_spent_exclusion_only` CAID3 state and
  unopened DP94 label seal, a versioned family manifest, hashed
  train/development/auxiliary files, source-specific `1/0/x` semantics,
  provenance, and exact/containment/local-similarity exclusions. It checks
  cross-role family and exact-sequence overlap and requires observed DFL
  positives in train and development. The similarity rule is BLOSUM62 with
  gap 10/1, at least 50 aligned residues, 25% identity and 50% coverage of
  the shorter sequence. It opens **no** benchmark label file.
- Passing this script proves only consistency with the supplied family map.
  A self-declared or poor family map is not a biological independence proof.
  The clustering source, parameters, input SHA-256, version and family/domain
  audit still need external validation. The source-specific negative semantics
  also need evidence, not just a manifest assertion.

The old `caid3_protected_dfl_v1` split fails closed. The read-only detector
status check identifies a local-similarity collision between train protein
`DP00919` and candidate holdout `DP03956`. DP94 remains a sequence-only
*candidate*, not a certified DFL benchmark. Neither that candidate's labels
nor the spent CAID3 labels/scores are available for model selection.

The 18 October source-capacity audit found at most 40 unique positive
sequences after direct screening of four candidate training sources, and 5
outside benchmark-touching components under the conservative transitive local
graph. This is not a validated family split. See
`reports/DFL_B1_FEASIBILITY_AUDIT_2026-10-18.md` for hashes, source tiers and
the alternative annotated-data route. **B1 remains open.**

That route now has a versioned local candidate artifact: 884 current UniProt
full-protein sequences matched their explicit flexible/disordered-linker region
peptides exactly, with all other residues masked unknown. Only 11 sequences
avoid direct CAID3/DP94 screening, one from the PubMed-supported tier. These
are not adjudicated DFL gold or negatives. The script, schema and hashes are
in the B1 feasibility report; neither training nor a split was launched.

The paper survey and sequence-only sensitivity audit in
`reports/DFL_B1_PAPER_SURVEY_AND_SPLIT_DESIGN_2026-10-18.md` diagnose the
present local rule as a poor stand-alone family definition. A *proposed*,
unvalidated reciprocal full-length screen flags 37/279 unique positive
candidate sequences against CAID3/DP94, versus 243/279 under the binding local
rule. The 206 changed decisions are not cleared for training: shared-domain
homology, low-complexity effects and the benchmark contract require a separate
family review. B1 and the current preflight remain unchanged and blocked.

## Reproduce the implemented checks

From the repository root, with `parasail` installed in the selected Python
environment:

```powershell
python -m unittest tests.test_preflight_dfl_family_training tests.test_predict_dfl_frozen -v
python scripts/preflight_dfl_family_training.py path/to/new_split/manifest.json
```

The tests create only synthetic benchmark sequences and masks. A command that
prints `pass_declared_sequence_and_family_manifest_only` is **not** a launch
authorization. The real trainer must invoke the new gate *and* the active
detector status check before GPU allocation, after a qualified versioned split
and all nominated benchmark registries exist. No such split is approved yet.

## Work remaining in this first stage

1. Select and document a DFL-labelled training source with verified positive,
   negative and unknown semantics, source release, rights and accession lineage.
   A positive-only source may contribute positives/unknowns; its unannotated
   residues must not become negatives. DLD labels are not DFL ground truth.
2. Freeze sequence-only benchmark registries and an independently justified
   domain/family clustering method. Hash inputs and cluster membership, then
   build a new train/development/auxiliary split and run the preflight. The
   current split cannot be re-labelled as passing.
3. Establish an independently adjudicated no-linker set, disjoint by family.
   Without it, specificity and large-scale mining remain unmeasured.
4. Wire preflight into a bounded training launcher. Retrain matched baselines
   and proposals on identical masks, encoder, scoring and coverage. Select
   only on development; repeat promising arms across seeds.
5. Calibrate abstention on a clean family-disjoint calibration subset. Freeze
   model, threshold, evaluator and competitor coverage before a one-time
   untouched comparison. Only then decide whether mining is justified.

No classifier breakthrough or published-method superiority has been
demonstrated. The one-time CAID3 retrospective is spent and must not be rerun
or used for tuning. See `RESEARCH_HANDOFF.md` and
`literature/DFL_UNTOUCHED_COMPARISON_CONTRACT_2026-10-17.md`.

The 18 October paper survey and hash-checked full-length/domain-family
sensitivity audits are in
`reports/DFL_B1_PAPER_SURVEY_AND_SPLIT_DESIGN_2026-10-18.md`. They found 134
exact-metadata DFL-positive candidate sequences in 112 components outside
both sequence-only benchmark registries under the combined exploratory family
graph. This is not yet a defensible train/development split or permission to
replace the binding firewall; B1 remains open.

The binding-rule/domain-family intersection is narrower still: 24 positive
sequences in 24 components pass both current screens. The hash-pinned capacity
and source-mask audit is documented in the same report. This is too small for
a credible high-precision mining claim, and no independently adjudicated
no-linker negatives exist in that set. The preflight now distinguishes the
historically viewed CAID3 exclusion registry from DP94's unopened local
label seal; it does not remove either sequence firewall.

An additional provenance-first DisProt 2026_06 assembly is documented in
`reports/DFL_CURATED_TRAINING_SOURCE_ASSEMBLY_2026-10-18.md`. Its strict
validated/cited/released, nonflagged tier yields only **two** screened
proteins with 36 observed positive residues and **zero** explicit ordered
residues. A separate 15-protein flagged tier is quarantined for evidence
review, not admitted to training. The tracked builder and JSON Schema make
the candidate artifact reproducible; they do not close B1.

A [100-paper title/abstract source screen](reports/DFL_100_PAPER_SOURCE_MERGE_2026-10-18.md)
has since produced a reproducible priority list, not 100 usable DFL datasets.
Four pinned row-level evidence streams were exact-sequence deduplicated under
the unchanged two-registry firewall into 24 positive/unknown candidates with
1,123 observed positive residues. All 200 unverified source-zero observations
remain outside the merged labels. This still lacks an adjudicated no-linker
panel and a certified family-disjoint train/development split; **B1 is open**.

## Firewall V2 Preregistration (`RULE_V2`) and Split Verification (18 October)

To resolve the V1 sequence-firewall deadlock documented above, the **Firewall V2 Preregistration Contract** ([`literature/DFL_FIREWALL_V2_PREREGISTRATION_2026-10-18.md`](literature/DFL_FIREWALL_V2_PREREGISTRATION_2026-10-18.md)) formalizes two auditable fixes enforced by [`schemas/dfl_family_split_v2.schema.json`](schemas/dfl_family_split_v2.schema.json), [`scripts/preflight_dfl_family_training_v2.py`](scripts/preflight_dfl_family_training_v2.py), and [`scripts/build_dfl_b1_v2_split.py`](scripts/build_dfl_b1_v2_split.py):

1. **True Ungapped Reciprocal Coverage (`RULE_V2`) + Shuffle-Null Gate ($\le 0.5\%$):** Replaces V1's unbounded gapped-length ratio (`len(traceback.query) / min(len)`, which had a `66.1%` residue-shuffle false-hit ratio across 4,800 pairs) with $\ge 30\%$ identity, $\ge 50$ aligned columns, and $\ge 80\%$ true ungapped residue coverage on **both** sequences (`0.000%` residue-shuffle and `0.021%` 12-mer block-shuffle false-hit rates), combined with transitive `Pfam` + `Gene3D` + `UniRef50` connected-component disjunction.
2. **Released Curator-Validated DisProt Evidence (`IDPO:0000033` / `IDPO:0000006`):** Admits non-obsolete, curator-validated, PMID/DOI-cited regions publicly released in `DisProt 2022_03`–`2026_06` that retained DisProt's internal `unpublished: true` legacy flag, recovering `97` benchmark-disjoint proteins across `91` family components (`60` train, `28` development, `9` auxiliary order calibration panel; `1,920` positive `1` residues and `1,253` verified `0` negative residues).
3. **Full `B1` Gate & Launcher Wiring:** [`scripts/preflight_caid3_protected_training.py`](scripts/preflight_caid3_protected_training.py), [`scripts/check_detector_autoresearch.py`](scripts/check_detector_autoresearch.py) (via [`configs/detector_autoresearch_active.json`](configs/detector_autoresearch_active.json) and [`reports/dfl_b1_v2_family_round_2026-10-18.json`](reports/dfl_b1_v2_family_round_2026-10-18.json)), and [`scripts/run_caid3_protected_dfl_trial.py`](scripts/run_caid3_protected_dfl_trial.py) (via [`configs/dfl_b1_v2_family_round_2026-10-18.json`](configs/dfl_b1_v2_family_round_2026-10-18.json)) all pass on `data/processed/dfl_b1_v2_family_split/manifest.json` (`d62136ce699294a0c07f0fa01a5350d3571643796cedadaf37d8399824893452`) while explicitly prohibiting `caid3_protected_dfl_v1`. See [`reports/DFL_B1_V2_PREREGISTRATION_AND_X1_TRANSFER_2026-10-18.md`](reports/DFL_B1_V2_PREREGISTRATION_AND_X1_TRANSFER_2026-10-18.md) and [`reports/dfl_b1_v2_family_split_verification_2026-10-18.json`](reports/dfl_b1_v2_family_split_verification_2026-10-18.json).

## Part 1 Completion (`B2` / `B3`) and Provenance-Stratified Mined & Natural Linker Dataset (`v1`) (19 October)

1. **Two-Stage Domain-Gated Conformal Calling (`B2` & `B3`):** Implemented in [`src/linkerforge/detector/domain_gate.py`](src/linkerforge/detector/domain_gate.py) and [`src/linkerforge/detector/conformal_mining.py`](src/linkerforge/detector/conformal_mining.py). Combining Stage 1 non-overlapping hetero-domain ($5\text{--}60\text{ aa}$, `left.accession != right.accession`) gating with Stage 2 frozen Seed-31 `ThreeLayerDFLCNN` (`c0d356c55ed83807353e66b875501c61b17f0b9b83f2910c0f1147f6409fdf6f`) calibration achieves `1.06%` empirical FDR (`98.94%` empirical precision, `93` TP / `1` FP) and `4.95%` 95%-Clopper–Pearson UCB FDR on `development` (`B2`), and `97.50%` protein-level specificity (`100.0%` on single/homo-domain proteins) across the 40-protein benchmark-disjoint No-Linker Specificity Panel (`B3`).
2. **Mined and Natural Linker Dataset (`mined_and_natural_linkers_v1`):** Built and verified via [`scripts/build_mined_and_natural_linker_dataset.py`](scripts/build_mined_and_natural_linker_dataset.py), validated against [`schemas/mined_and_natural_linker_dataset_v1.schema.json`](schemas/mined_and_natural_linker_dataset_v1.schema.json), and exposed to Part 2 via [`src/linkerforge/dataset_v1.py`](src/linkerforge/dataset_v1.py). Contains `509` linker records ($5\text{--}60\text{ aa}$) across four provenance-stratified tiers (`91` `natural_cath_pdb_verified`, `57` `natural_disprot_curated_dfl`, `88` `mined_high_confidence_dfl_call`, `250` `mined_provisional_dfl_candidate`, `23` `literature_engineered_template`) plus `40` `b3_no_linker_specificity_panel` controls. See [`reports/PART1_COMPLETION_AND_NATURAL_MINED_LINKER_DATASET_2026-10-19.md`](reports/PART1_COMPLETION_AND_NATURAL_MINED_LINKER_DATASET_2026-10-19.md) and [`reports/mined_and_natural_linkers_v1_summary_2026-10-19.json`](reports/mined_and_natural_linkers_v1_summary_2026-10-19.json).
