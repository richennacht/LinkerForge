# B1 V2 three-layer CNN: three-seed development and precision audit

Status: **development-only, not a certified 99%-precision linker miner**.
The family-disjoint B1 V2 split was rebuilt from the recorded sources, and the
preflight passed before each of the two new GPU runs. All three seeds use the
same 60 training and 28 development proteins, frozen ESM-2 35M encoder,
1/5/1-kernel 50/5-channel CNN, 12 epochs, masked `1`/`0` training with `x`
ignored, and development-APS epoch selection. The nine-protein auxiliary
order panel and CAID3/DP94 labels were not used for fitting or selection.

Split manifest SHA-256:
`d62136ce699294a0c07f0fa01a5350d3571643796cedadaf37d8399824893452`.
Family manifest SHA-256:
`db484272985535f08cf182166503973e77753abd00c7a4df1a197cc9db971dd2`.
The seed-repeat trainer source SHA-256 was
`e8c6c516e8d722d2d3bd561994fcbedd54ea8533ea9a4fbf592b8e45f8d50e05`;
its training commit was `92eeb948b6682a74eca3731f188b758364e21c1a`.
The seed-7 runner hash/commit are retained in the per-run manifest and round
ledger. Raw epoch logs, checkpoint hashes, config hashes, runtime, and peak
PyTorch VRAM are recorded in `reports/dfl_b1_v2_cnn_round_2026-10-18.json`.

| Seed | Best epoch | Dev APS | Dev AUROC | Dev Fmax | TP/FP at >=99% observed precision | Known-positive recall | Threshold |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 2 | 0.770304 | 0.706613 | 0.688020 | 100/1 | 0.182815 | 0.648354 |
| 17 | 5 | 0.791283 | 0.730852 | 0.693603 | 139/1 | 0.254113 | 0.563638 |
| 31 | 4 | 0.817759 | 0.781327 | 0.725806 | 153/1 | 0.279707 | 0.735190 |

Each high-precision threshold was chosen to maximize *known-positive* recall
on the same 28-protein development set, with score ties kept intact and at
least 50 true-positive residues required. There are 547 labelled positives
and 518 labelled negatives in that set. The empirical observed-label
precisions are 99.010%, 99.286%, and 99.351% respectively, but this is
selection and measurement on the same development sample. Even an exact
binomial lower bound treating each threshold as fixed (which it was **not**)
would only be 95.39%, 96.66%, and 96.96% one-sided at 95% confidence.
The 9,193 unknown residues in the scored development windows are excluded;
they may include false positives. Curated order residues are not a
representative sample of all non-linker residues, and no independently
adjudicated no-linker-protein set or sealed external DFL test was scored.
Therefore **none of these numbers validates 99% precision for proteome-scale
mining or superiority over published methods**. The three different cutoffs
also show that a portable threshold has not yet been calibrated.

Reproduction: `scripts/audit_dfl_b1_v2_precision.py` rechecks the split and
checkpoint hashes, recomputes development APS from the selected checkpoint,
and writes isolated audit JSON without opening benchmark labels. CPU
recomputation differed from CUDA training APS by at most about 2.2e-5; the
audit allows 1e-4 for this backend numerical effect. The three precision
audit JSON files remain ignored local artifacts, while this report and the
round ledger carry the important counts and hashes. The precision-frontier
unit tests cover ties, minimum TP support, and the zero-FP confidence limit.

Next gate: obtain a separately adjudicated, family-audited DFL-positive and
no-linker evaluation set; freeze a single model, threshold, and evaluator
before opening its labels. If that set does not support 99% precision with
useful recall and uncertainty, do not mine proteins under a 99% claim.
