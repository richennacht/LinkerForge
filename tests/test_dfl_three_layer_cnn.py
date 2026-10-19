"""Synthetic-only architecture and mask tests. No biological labels are used."""

import unittest

import torch

from src.linkerforge.detector.cnn import ThreeLayerDFLCNN, labels_to_tensor, masked_bce_with_logits
from src.linkerforge.detector.sequence_cnn import sequence_logits


class ThreeLayerCNNTests(unittest.TestCase):
    def test_variable_lengths_and_three_convolutions(self):
        model = ThreeLayerDFLCNN(8)
        self.assertEqual(sum(isinstance(layer, torch.nn.Conv1d) for layer in model.modules()), 3)
        self.assertEqual(tuple(model(torch.randn(2, 13, 8)).shape), (2, 13))
        self.assertEqual(tuple(model(torch.randn(1, 1, 8)).shape), (1, 1))
        with self.assertRaises(ValueError):
            model(torch.randn(2, 13, 7))

    def test_unknown_and_padding_have_no_loss_or_gradient(self):
        labels = labels_to_tensor(["10xx", "x01"], 5)
        self.assertEqual(labels.tolist(), [[1, 0, -1, -1, -1], [-1, 0, 1, -1, -1]])
        logits = torch.zeros(2, 5, requires_grad=True)
        loss = masked_bce_with_logits(logits, labels)
        loss.backward()
        self.assertEqual(logits.grad[0, 2:].tolist(), [0.0, 0.0, 0.0])
        self.assertEqual(logits.grad[1, [0, 3, 4]].tolist(), [0.0, 0.0, 0.0])
        changed = logits.detach().clone()
        changed[labels == -1] = 1000
        self.assertEqual(float(loss.detach()), float(masked_bce_with_logits(changed, labels)))

    def test_invalid_or_only_unknown_masks_fail_closed(self):
        with self.assertRaises(ValueError):
            labels_to_tensor(["10?"], 3)
        with self.assertRaises(ValueError):
            masked_bce_with_logits(torch.zeros(1, 3), labels_to_tensor(["xxx"], 3))
        with self.assertRaises(ValueError):
            masked_bce_with_logits(torch.zeros(1, 3), torch.tensor([[1, 2, -1]]))

    def test_cpu_training_step_reduces_synthetic_loss(self):
        torch.manual_seed(17)
        model = ThreeLayerDFLCNN(4)
        embeddings = torch.randn(4, 12, 4)
        masks = ["1" * 6 + "0" * 6] * 4
        labels = labels_to_tensor(masks, 12)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.02)
        first = float(masked_bce_with_logits(model(embeddings), labels).detach())
        for _ in range(35):
            optimizer.zero_grad()
            loss = masked_bce_with_logits(model(embeddings), labels)
            loss.backward()
            optimizer.step()
        last = float(masked_bce_with_logits(model(embeddings), labels).detach())
        self.assertLess(last, first)

    def test_sequence_adapter_preserves_residue_alignment(self):
        class Tokenizer:
            def __call__(self, sequences, **kwargs):
                length = max(map(len, sequences)) + 2
                return {"input_ids": torch.zeros(len(sequences), length, dtype=torch.long)}

        class Encoder(torch.nn.Module):
            def forward(self, input_ids):
                return type("Output", (), {"last_hidden_state":
                            torch.ones(*input_ids.shape, 8)})()

        head = ThreeLayerDFLCNN(8)
        scores = sequence_logits(Encoder(), Tokenizer(), head, ["ACDEF", "MGS"])
        self.assertEqual([len(row) for row in scores], [5, 3])
        with self.assertRaises(ValueError):
            sequence_logits(Encoder(), Tokenizer(), head, ["ACD?F"])


if __name__ == "__main__":
    unittest.main()
