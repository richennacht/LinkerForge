# B1 V2 three-layer CNN: local reproduction and launch record (18 Oct 2026)

Status at initial commit: **preflight passed; one bounded development-only GPU
trial prepared**. This is not a certified high-precision natural-linker miner
and no protected benchmark labels or scores are used.

The ignored V2 split described in PR #8 was absent from this checkout. It was
regenerated with `scripts/build_dfl_b1_v2_split.py` from the pinned local
DisProt 2026_06, CAID3 *sequence-only*, and DP94 *sequence-only* sources.
The reconstructed split has the same 60/28/9 protein counts, 91 family
components, 1,920 DFL-positive residues and 1,253 curated-order negative
residues reported in PR #8. Its **manifest hash differs** from the PR report
because each row records the generating checkout commit and Python version.
This reproduction was generated at commit
`353bbcab83ae5148bf8ee7835544869876d201d9` with Python 3.12.14; the
reproduced split manifest SHA-256 is
`d62136ce699294a0c07f0fa01a5350d3571643796cedadaf37d8399824893452`.
The family manifest SHA-256 remains
`db484272985535f08cf182166503973e77753abd00c7a4df1a197cc9db971dd2`.
No original split file was overwritten. Its data remain ignored by Git.

The reproduced root CAID3 and DP94 registry copies were checked by
`scripts/preflight_caid3_protected_training.py`; result:
`pass_v2_preregistered_sequence_and_family_firewall` with 60 train, 28
development, and 9 auxiliary proteins. `scripts/check_detector_autoresearch.py`
now returns `launch_or_audit` for the reproduced V2 split. The old
`caid3_protected_dfl_v1` split remains prohibited. The auxiliary order panel
is hash-checked but not trained on or used for model selection in this trial.

The generic SAM/YOLO launcher in PR #8 is **not** used here: its inherited
auxiliary-data check expects the older fully observed disorder-transfer masks
and would reject this V2 order panel. Instead, the dedicated
`scripts/train_dfl_b1_v2_three_layer_cnn.py` runs the requested 1/5/1-kernel,
50/5-channel three-convolution head on a frozen local ESM-2 encoder. It trains
only on verified `1` and `0` residues, ignores `x` in loss and metrics,
selects checkpoints by masked development APS, records AUROC/Fmax and raw
label counts, and limits the run to 12 epochs or 1,800 seconds. Its config
pins the regenerated split/role hashes and uses seed 7. It refuses an existing
output or concurrent internal B1 lock and loads ESM weights offline.

Validation before launch: 5 existing CNN tests, 5 V2 firewall tests, and 2
new window-mask tests passed. The CLI `--preflight-only` passed with 60 train
and 28 development proteins. Local CUDA is available and the frozen ESM-2
checkpoint loads from cache. `nvidia-smi` reported no other GPU Python job.

This is one baseline development trial, not a published-method comparison.
The auxiliary order panel is not a full-proteome no-linker calibration set;
conformal abstention and an independently adjudicated no-linker specificity
panel remain required before a precision or mining claim.

## Seed-7 completion

The hidden independent GPU process completed all 12 epochs in 26.312 seconds
(peak PyTorch allocation 186,616,832 bytes) and exited. Its best checkpoint
was selected at epoch 2 by masked development APS. On the 547 annotated
positive and 518 annotated-order negative development residues, the best
checkpoint recorded APS **0.770304**, AUROC **0.706613**, and 101-threshold-grid
Fmax **0.688020**. The isolated checkpoint SHA-256 is
`a72e2f00f0a0becf5460a24366c668bf877d916fff25f29f33bd500c1f2359a5`.
The checkpoint and raw epoch log remain ignored local artifacts; the tracked
round ledger records their identity. This is one-seed development evidence,
not a protected-benchmark score, calibrated 90%-precision result, or
permission to run proteome-scale mining.
