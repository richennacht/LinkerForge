"""Synthetic checks for fail-closed source-mask capacity accounting."""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.audit_dfl_binding_survivor_masks import audit


class BindingSurvivorMaskTests(unittest.TestCase):
    def test_keeps_unverified_zero_and_positive_conflict_visible(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "capacity.json"
            sequence = "ACD"
            artifact = {
                "status": "screened_capacity_not_training_split",
                "candidates": [{"record_id": "C0001", "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                                "source_ids": [{"source": "published_dfl_train", "protein_id": "DP1"},
                                               {"source": "disprot_positive_only", "protein_id": "DP1"}]}],
            }
            path.write_text(json.dumps(artifact), encoding="utf-8")
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            rows = {("published_dfl_train", "DP1"): {"sequence": sequence, "labels": "1x0"},
                    ("disprot_positive_only", "DP1"): {"sequence": sequence, "labels": "xx1"}}
            with patch("scripts.audit_dfl_binding_survivor_masks.source_rows", return_value=(rows, {})):
                result = audit(Path(directory), path, digest)
            self.assertEqual(result["counts"]["survivor_sequences"], 1)
            self.assertEqual(result["counts"]["published_unverified_zero_residues"], 1)
            self.assertEqual(result["counts"]["cross_source_positive_zero_conflicts"], 1)

    def test_rejects_changed_capacity_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "capacity.json"
            path.write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                audit(Path(directory), path, "0" * 64)


if __name__ == "__main__":
    unittest.main()
