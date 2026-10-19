"""CPU-only synthetic forward/backward check of the first-stage DFL CNN.

This command uses no proteins, benchmark labels, pretrained encoder or GPU.
It cannot create a mining-qualified checkpoint.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.linkerforge.detector.cnn import ThreeLayerDFLCNN, labels_to_tensor, masked_bce_with_logits
from src.linkerforge.detector.sequence_cnn import sequence_logits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--with-cached-esm", action="store_true",
                        help="Also verify sequence-to-logit alignment with locally cached ESM-2; no download")
    args = parser.parse_args()
    torch.manual_seed(20261018)
    model = ThreeLayerDFLCNN(8)
    features = torch.randn(2, 16, 8)
    labels = labels_to_tensor(["xx111000xxxxxxxx", "x0011xxxxxxxxxxx"], 16)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    before = float(masked_bce_with_logits(model(features), labels).detach())
    for _ in range(12):
        optimizer.zero_grad(set_to_none=True)
        loss = masked_bce_with_logits(model(features), labels)
        loss.backward()
        optimizer.step()
    after = float(masked_bce_with_logits(model(features), labels).detach())
    if not after < before:
        raise RuntimeError("Synthetic loss did not decrease")
    result = {"status": "synthetic_cpu_smoke_pass", "architecture": "three_layer_1d_cnn_1_5_1",
                      "scores_per_sequence": list(model(features).shape),
                      "observed_positives": int((labels == 1).sum()),
                      "observed_negatives": int((labels == 0).sum()),
                      "ignored_unknown_or_padding": int((labels == -1).sum()),
                      "loss_before": before, "loss_after": after,
                      "scientific_status": "untrained_on_proteins_no_precision_claim"}
    if args.with_cached_esm:
        from transformers import AutoModel, AutoTokenizer
        encoder_id = "facebook/esm2_t12_35M_UR50D"
        tokenizer = AutoTokenizer.from_pretrained(encoder_id, local_files_only=True)
        backbone = AutoModel.from_pretrained(encoder_id, local_files_only=True).cpu()
        esm_head = ThreeLayerDFLCNN(backbone.config.hidden_size)
        sequences = ["ACDEFGHIKLMNPQRSTVWY", "MGSAGGSGGSA"]
        logits = sequence_logits(backbone, tokenizer, esm_head, sequences)
        result["cached_esm_sequence_lengths"] = [len(row) for row in logits]
        result["encoder_id"] = encoder_id
        result["resolved_encoder_revision"] = getattr(backbone.config, "_commit_hash", None)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
