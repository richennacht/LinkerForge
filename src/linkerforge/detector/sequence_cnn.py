"""Sequence-to-logit adapter for a frozen ESM-2 encoder and the DFL CNN.

No threshold, span call, or mining claim is made by this uncalibrated adapter.
The encoder and tokenizer are supplied by the caller and must be version-pinned.
"""

from __future__ import annotations

import torch

from .cnn import ThreeLayerDFLCNN


def sequence_logits(backbone, tokenizer, head: ThreeLayerDFLCNN,
                    sequences: list[str], *, device: torch.device | str = "cpu") -> list[torch.Tensor]:
    if not sequences or any(not sequence or len(sequence) > 510 for sequence in sequences):
        raise ValueError("Provide nonempty sequences of at most 510 residues")
    if any(set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO") for sequence in sequences):
        raise ValueError("Sequences must be uppercase amino-acid strings")
    tokens = tokenizer(sequences, return_tensors="pt", padding=True,
                       add_special_tokens=True, truncation=False)
    tokens = {name: value.to(device) for name, value in tokens.items()}
    backbone.eval()
    with torch.no_grad():
        hidden = backbone(**tokens).last_hidden_state
    length = max(map(len, sequences))
    if hidden.ndim != 3 or hidden.shape[0] != len(sequences) or hidden.shape[1] < length + 2:
        raise ValueError("ESM-2 output does not align one-to-one with residues plus BOS/EOS")
    if hidden.shape[2] != head.embedding_dim:
        raise ValueError("Encoder and CNN embedding dimensions differ")
    logits = head(hidden[:, 1:1 + length])
    return [logits[index, :len(sequence)] for index, sequence in enumerate(sequences)]
