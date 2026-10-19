# Frozen DFL classifier: provisional inference

`scripts/predict_dfl_frozen.py` is the first callable classifier interface for
the sequence-to-linker-to-design workflow. It loads the previously selected
SAM-like DFL checkpoint only after verifying the checkpoint and manifest hashes
in `reports/dfl_sam_repeat_frozen_selection_2026-10-17.json`. It does not train,
use CAID3 labels, or evaluate a benchmark. The checkpoint is a historical
development selection, **not** a certified mining model.

Example, on the explicitly synthetic fixture (requires a local checkpoint and
cached encoder; these files are not part of a clean clone):

```powershell
.venv-gpu\Scripts\python.exe scripts\predict_dfl_frozen.py `
  --fasta tests\fixtures\synthetic_dfl_inference_smoke.fasta `
  --selection reports\dfl_sam_repeat_frozen_selection_2026-10-17.json `
  --encoder-revision 6fbf070e65b0b7291e7bbcd451118c216cff79d8 `
  --device cpu
```

The default returns one score per residue and **abstains from calling spans**.
An explicit `--threshold` returns 1-based inclusive *candidate* spans; it is
not a precision guarantee or a verified annotation. No command in this document
may be run on CAID3 or another nominated benchmark. The encoder revision was
observed in the local cache; the original training manifest did not pin it, so
the report records the revision used at inference and does not assert exact
training equivalence. Inputs and outputs are hash-labelled. The script refuses
overwriting an existing JSON output.

This interface is intentionally upstream of the design workflow:

1. A future, family-disjoint DFL classifier must pass the full benchmark
   firewall and an independently adjudicated no-linker test before it may mine
   unannotated proteins. The current checkpoint may not be promoted for mining.
2. Only provenance-bearing, reviewed mined linkers can enter a natural-linker
   retrieval set. Predicted spans are not function-tested linkers.
3. The deterministic geometry calculator from PR #3 supplies attachment-site
   feasibility and length windows. The PR #2 template library supplies an
   auditable candidate source. A matched retrieval baseline must select a
   training-set linker using **flanking context only**, with whole family/domain
   pairs held out; report length and sequence-recovery metrics alongside the
   template and geometry approaches. Do not claim retrieval beats generation
   until this test exists.

Current gate: the old protected split overlaps the proposed DP94 holdout and
cannot be used for new training. The CAID3 comparison is spent retrospective.
This entry point changes neither fact. See
`literature/DFL_UNTOUCHED_COMPARISON_CONTRACT_2026-10-17.md`.
