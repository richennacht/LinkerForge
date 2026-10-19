# B1 classifier training-start gate — 18 October 2026

Decision: **do not allocate GPU or claim that B1 is train-ready.** The
three-layer CNN implementation exists, but the current real DFL data do not
support a defensible high-precision binary training/calibration run. This is
a data and evaluation gate, not an architecture failure. No benchmark labels,
per-protein benchmark errors, or public CAID3 scores were used.

## Reproducible checks

`scripts/check_dfl_b1_training_readiness.py` re-hashes and validates the
deduplicated candidate pool, counts known and unknown masks, and exits 2 on
the actual local data. The input manifest SHA-256 is
`9e0def03d4d3d95f9cb4ddea6841b9e8973070868e5e052cc594f278069c1299`;
the records SHA-256 is
`0f107d08fde1563fc2dd7a523f583eedde31939e1e81f5d5b74c44e4f63a5121`.
It finds 24 unique candidate sequences in 24 exploratory components, 1,123
observed DFL-positive residues, **zero verified negative residues**, and
12,721 unknown residues. It rejects tampered records and a 1/x candidate
whose unknowns have been rewritten as zeros; three unit tests pass with the
local Python environment. The pool has no frozen family-disjoint
train/development allocation, independently adjudicated no-linker panel, or
passing split-specific benchmark firewall. Component IDs alone do not
certify biological/domain-family independence.

An independent read-only source-capacity run of
`scripts/audit_dfl_b1_source_capacity.py --include-idpfunnet` used the pinned
IDPFunNet author training file SHA-256
`89e374e30924a277fb4bdcf1b7a82776017e4ab70397757b185fa9d0b16954b2`.
Its 552 records include 124 with an observed DFL-positive residue, of which
17 positive sequences are direct non-hits against the pinned CAID3 and DP94
sequence registries under the declared local rule. The combined source audit
has 40 unique direct non-hit positive sequences; only **five** positive
sequences lie outside benchmark-connected components in its conservative
transitive local-similarity graph. This graph is not a validated family
partition and is a sensitivity bound, not a replacement for the existing
Pfam/Gene3D/UniRef50 audit. It cannot establish a safe training split. The
new source's redistribution rights and negative-mask semantics remain
unverified, so its zeros are not training negatives. The full source audit
read benchmark IDs and sequences only.

The separate published-negative-capacity audit was rerun. Of 84 TR166
nonpositive records, 83 directly collide with the protected registries under
the binding local rule; the one direct survivor has no independently verified
negative semantics and has not passed a transitive family audit. Its 16,415
published `0` observations are therefore **not** a reservoir of 16,415
verified negative residues. This closes the tempting fallback of importing
the published non-DFL track wholesale.

The existing `scripts/preflight_dfl_family_training.py` is the required
sequence/family manifest gate for any future versioned split. Passing that
gate alone would still be insufficient for a high-precision caller trained
only on positives: precision and false-positive rates need independently
adjudicated negatives and a no-linker calibration panel. The old protected
split fails the DP94 sequence firewall and was not retried for training.

## Smallest valid route to a first training run

1. Secure source-level evidence for `0` meaning verified non-DFL (rather than
   merely unannotated or a different function), with rights and accession
   provenance. Preserve unresolved residues as `x`.
2. Build a separately adjudicated no-linker panel to set a high-precision
   operating threshold on development data; never mine negatives from the
   candidate pool's `x` mask.
3. Create a versioned family/domain-disjoint train/development split after
   excluding every nominated benchmark sequence and family. Keep benchmark
   labels sealed.
4. Pass a hash-pinned full-registry preflight before allocating the single
   local GPU; then run the existing three-convolution baseline with bounded
   runtime, isolated outputs, multiple seeds, and development-only selection.

Until these conditions hold, beginning a real B1 training run would produce
misleading probabilities, not a high-precision natural-linker miner. This
check deliberately did **not** launch a synthetic or positive-only run and
does not imply a publishable model exists.
