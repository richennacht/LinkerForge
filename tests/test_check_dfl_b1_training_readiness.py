import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from scripts.check_dfl_b1_training_readiness import audit


def _pool(tmp_path, labels="1xx"):
    record = {
        "record_type": "deduplicated_dfl_positive_candidate",
        "training_eligibility": "candidate_not_B1_no_adjudicated_negatives_or_frozen_split",
        "sequence": "AGS",
        "sequence_sha256": hashlib.sha256(b"AGS").hexdigest(),
        "labels": labels,
        "source_observations": [{"source": "curated", "source_id": "example"}],
        "family_component_id": "candidate-component-1",
    }
    raw = json.dumps(record).encode() + b"\n"
    (tmp_path / "positive_unknown_candidates.jsonl").write_bytes(raw)
    manifest = {
        "status": "deduplicated_positive_unknown_candidates_not_B1",
        "records_sha256": hashlib.sha256(raw).hexdigest(),
        "counts": {"unique_positive_sequences": 1,
                   "merged_positive_residues": labels.count("1"),
                   "merged_unknown_residues": labels.count("x")},
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


class ReadinessTests(unittest.TestCase):
    def test_positive_unknown_pool_does_not_become_training_ready(self):
        with tempfile.TemporaryDirectory() as directory:
            pool = Path(directory)
            _pool(pool)
            result = audit(pool)
            self.assertEqual(result["status"], "blocked_not_training_ready")
            self.assertEqual(result["residue_counts"], {"1": 1, "0": 0, "x": 2})
            self.assertFalse(result["gates"]["verified_negative_residues"])

    def test_tampered_pool_fails_before_census(self):
        with tempfile.TemporaryDirectory() as directory:
            pool = Path(directory)
            _pool(pool)
            (pool / "positive_unknown_candidates.jsonl").write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                audit(pool)

    def test_unknown_cannot_be_relabelled_zero_in_candidate_pool(self):
        with tempfile.TemporaryDirectory() as directory:
            pool = Path(directory)
            _pool(pool, "100")
            with self.assertRaisesRegex(ValueError, "Invalid or duplicate candidate record"):
                audit(pool)


if __name__ == "__main__":
    unittest.main()
