"""Three-convolution DFL residue head for frozen per-residue embeddings.

This is a trainable architecture, not a trained or mining-qualified model.
Inputs are (batch, residues, embedding channels); outputs are residue logits.
The 1/5/1 kernel and 50/5 hidden-channel layout follow the released
LINKER-Pred-Lite CNN family. No published weights are loaded here.
"""

from __future__ import annotations

import torch
from torch import nn
from torch.nn import functional as F


class ThreeLayerDFLCNN(nn.Module):
    def __init__(self, embedding_dim: int, first_width: int = 50, second_width: int = 5) -> None:
        super().__init__()
        if min(embedding_dim, first_width, second_width) < 1:
            raise ValueError("All channel widths must be positive")
        self.embedding_dim = embedding_dim
        self.conv1 = nn.Conv1d(embedding_dim, first_width, kernel_size=1)
        self.conv2 = nn.Conv1d(first_width, second_width, kernel_size=5, padding=2)
        self.conv3 = nn.Conv1d(second_width, 1, kernel_size=1)

    def forward(self, embeddings: torch.Tensor) -> torch.Tensor:
        if embeddings.ndim != 3 or embeddings.shape[-1] != self.embedding_dim:
            raise ValueError("Expected embeddings shaped (batch, residues, embedding_dim)")
        if embeddings.shape[0] < 1 or embeddings.shape[1] < 1:
            raise ValueError("Empty embedding batch or sequence")
        hidden = F.relu(self.conv1(embeddings.transpose(1, 2)))
        hidden = F.relu(self.conv2(hidden))
        return self.conv3(hidden).squeeze(1)


def labels_to_tensor(masks: list[str], width: int, *, device: torch.device | str = "cpu") -> torch.Tensor:
    """Encode 1/0/x masks as 1/0/-1; right padding is always unknown."""
    if not masks or width < 1:
        raise ValueError("Masks and width must be nonempty")
    labels = torch.full((len(masks), width), -1, dtype=torch.int8, device=device)
    for index, mask in enumerate(masks):
        if not mask or len(mask) > width or set(mask) - set("01x"):
            raise ValueError(f"Invalid 1/0/x mask at batch index {index}")
        labels[index, :len(mask)] = torch.tensor(
            [{"x": -1, "0": 0, "1": 1}[value] for value in mask],
            dtype=torch.int8, device=device,
        )
    return labels


def masked_bce_with_logits(logits: torch.Tensor, labels: torch.Tensor) -> torch.Tensor:
    """Average BCE over *verified* 0/1 residues; never train on x or padding."""
    if logits.ndim != 2 or labels.shape != logits.shape:
        raise ValueError("Logits and labels must have matching (batch, residues) shapes")
    if not torch.all((labels == -1) | (labels == 0) | (labels == 1)):
        raise ValueError("Labels must be -1 (unknown), 0, or 1")
    known = labels != -1
    if not torch.any(known):
        raise ValueError("No verified labels in batch")
    return F.binary_cross_entropy_with_logits(logits[known], labels[known].to(logits.dtype))
