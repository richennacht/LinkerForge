# First-stage three-layer DFL CNN

`src/linkerforge/detector/cnn.py` now implements the small, three-convolution
residue head specified in `ARCHITECTURE_HANDOFF.md` section 3. It accepts
frozen per-residue protein-language-model embeddings in `(batch, length,
channels)` order and returns one **logit**, not a calibrated probability, per
residue. The 1/5/1 convolution kernel and 50/5 hidden-channel layout match the
released LINKER-Pred-Lite CNN family, but this is a new untrained head with no
published weights. It can accept ESM-2 embeddings by setting their actual
channel dimension; no encoder or feature archive is silently downloaded.
`src/linkerforge/detector/sequence_cnn.py` connects an already loaded frozen
ESM-2 encoder and tokenizer to the head and returns one logit per input residue.
It bounds input to 510 residues and checks encoder/head channel alignment.

The same module provides strict `1/0/x` mask encoding and a masked BCE loss.
Only verified `1` and `0` residues contribute to training; `x` and right
padding contribute neither loss nor gradient. A batch with no known labels
fails closed. The head handles variable sequence lengths and exposes logits
for later, independently calibrated abstention and span decoding.

Run a CPU-only, synthetic end-to-end training smoke:

```powershell
.venv-gpu\Scripts\python.exe scripts\smoke_dfl_three_layer_cnn.py
.venv-gpu\Scripts\python.exe scripts\smoke_dfl_three_layer_cnn.py --with-cached-esm
.venv-gpu\Scripts\python.exe -m unittest tests.test_dfl_three_layer_cnn -v
```

These commands validate shape, forward/backward optimization and mask
semantics. They do **not** train on proteins, produce a checkpoint, calibrate a
threshold or validate biological performance. The cached-ESM option is offline
and fails if that model is not present locally;
it verified score lengths 20 and 11 for two synthetic amino-acid sequences on
the configured Windows environment. The historical
`scripts/train_dfl_masked_cnn.py` is pinned to the old TR166/TE82 split and is
**not** an eligible route to a new mining model. A real bounded trainer must
first invoke the versioned family/benchmark firewall and the active detector
integrity check, and needs an independently justified family map, verified
negative semantics and a defensible no-linker calibration set. The current
source pool fails those requirements; B1 in `CHECKLIST.md` remains open.

The executable published-negative capacity audit is
`scripts/audit_dfl_published_negative_capacity.py`. Against the pinned local
sources on 18 October it found 84 published nonpositive records but only **one**
outside the binding direct CAID3/DP94 sequence screen. The source contains
16,415 published zero residues, yet those zeros are not independently
adjudicated negatives and cannot be relabelled as a no-linker gold panel. This
confirms why simply training the now-working CNN on the available old records
would not complete B1 or support a mining precision claim. The audit parses
benchmark sequences only, not benchmark labels.
