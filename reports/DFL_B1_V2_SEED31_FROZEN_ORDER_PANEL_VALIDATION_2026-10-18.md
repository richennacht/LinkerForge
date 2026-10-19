# Frozen seed-31 validation: the 99% mining claim is not supported

The seed-31 CNN checkpoint was selected by the highest B1 V2 development APS
among seeds 7, 17, and 31. A deployment-equivalent scorer then processed **all**
sequence windows one at a time, independent of label locations, rather than
only windows containing known labels. On the 28 development proteins, the
known-label result remained 153 TP, 1 FP, 394 FN at threshold
`0.7351981997489929` (99.35% *empirical annotated-residue* precision and
27.97% known-positive recall). It also called 1,126 of 13,713 unknown
residues. These are not counted as correct or incorrect without adjudication.

Before scoring the auxiliary panel, the model checkpoint, run manifest,
split, threshold, development audit, and evaluator source were hash-frozen
in `reports/dfl_b1_v2_seed31_order_panel_selection_2026-10-18.json`, pushed
as commit `d73acf6`. The evaluator refuses a changed selection, checkpoint,
config, split, development audit, or source hash. Selection SHA-256:
`aeac492ca9dd7c50662c8a2e3b790f3cc2708b46c869f065afc81a5019a3ee1e`.
The auxiliary output SHA-256 is
`cbc770fda9d5fdcb83997d702ef31fa1e68c95c4d232f027a85e513463c268c1`;
raw output remains in ignored local artifacts.

The separate nine-protein, family-disjoint auxiliary panel contains 194
curator-labelled experimental-order residues, **no verified DFL-positive
residues**, and 3,347 unknown residues. At the frozen threshold, the model
called **14 of 194 known-order residues** (7.22% false-positive rate on those
residues; specificity 92.78%) across three proteins. It also called 364
unknown residues (10.88% of that unknown mask). No threshold or model was
changed after observing this result; this panel is now spent for selection.

This does **not** estimate precision: the panel has no verified positives,
unknown is not negative, and order-only proteins are not independently
adjudicated whole-protein no-linker examples. It does show that the
development-selected 99% operating point is not robust even on this limited
negative challenge. Thus the seed-31 checkpoint is **not certified for
proteome-scale mining**, and the three-seed annotated-development percentages
must not be advertised as validated 99% precision. The public CAID3 result
was not rerun, and DP94 labels were not opened.

Next validation gate: recruit or curate a prospectively sealed, family-audited
set with both DFL-positive spans and explicitly adjudicated no-linker
proteins/residues; freeze a future checkpoint, scorer, threshold, coverage
rules, and metrics before unsealing. Report raw TP/FP/FN, protein-level false
calls, precision-recall with uncertainty, and abstention/coverage. Until then,
the classifier can be used only as a provisional candidate generator for
manual review, not as a high-precision natural-linker database builder.
