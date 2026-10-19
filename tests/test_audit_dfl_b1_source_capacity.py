"""Synthetic B1-capacity audit tests; no real benchmark labels are opened."""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from scripts.audit_dfl_b1_source_capacity import audit


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value.encode() if isinstance(value, str) else json.dumps(value).encode()
    path.write_bytes(raw)
    return hashlib.sha256(raw).hexdigest()


class CapacityAuditTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        caid = "".join(f">C{i}\n{'M' * 80}\n{'1' * 80}\n" for i in range(31))
        caid_hash = write(self.root / "artifacts/caid3-linker.fasta", caid)
        pool = self.root / "data/processed/linker_pred2_style_pool_v1"
        files = {}
        for name, key, sequence, task in (("published_dfl_train", "T", "A" * 80, "DFL"),
                                          ("disprot_positive_only", "U", "C" * 80, "DFL_positive_unlabeled")):
            row = {"protein_id": key, "sequence": sequence, "labels": "1" + "x" * 79,
                   "task": task}
            files[name] = {"sha256": write(pool / f"{name}.jsonl", json.dumps(row) + "\n"), "records": 1}
        write(pool / "manifest.json", {
            "status": "prepared_not_linker_pred2_exact_split_not_training_cleared",
            "source_sha256": {"caid3": caid_hash}, "files": files})
        catalog = self.root / "data/processed/dfl_positive_candidates_2026_10_12.jsonl"
        row = {"protein_id": "V", "sequence": "G" * 80, "labels": "1" + "x" * 79,
               "record_type": "disprot_dfl_positive_candidate"}
        catalog_hash = write(catalog, json.dumps(row) + "\n")
        write(catalog.with_name("dfl_positive_candidates_2026_10_12.manifest.json"), {
            "status": "positive_unknown_candidates_not_trainable", "output_sha256": catalog_hash,
            "records": 1})
        dp = self.root / "data/processed/prospective_dfl_holdout_dp94_v1"
        fasta = ">D1\n" + "W" * 80 + "\n"
        fasta_hash = write(dp / "sequences.fasta", fasta)
        write(dp / "manifest.json", {"status": "candidate_not_certified_independent",
                                    "sequence_fasta_sha256": fasta_hash,
                                    "holdout_proteins": 1, "holdout_ids": ["D1"]})

    def test_counts_unique_unflagged_synthetic_positives(self):
        result = audit(self.root)
        self.assertEqual(result["status"], "sequence_screen_only_not_B1")
        self.assertEqual(result["combined_unique_surviving_sequences"], 3)
        self.assertEqual(result["declared_local_graph"]["positive_sequences_outside_benchmark_components"], 3)

    def test_tampered_source_hash_rejected(self):
        path = self.root / "data/processed/linker_pred2_style_pool_v1/published_dfl_train.jsonl"
        path.write_text("tampered")
        with self.assertRaisesRegex(ValueError, "Source hash changed"):
            audit(self.root)

    def test_dld_task_rejected(self):
        pool = self.root / "data/processed/linker_pred2_style_pool_v1"
        path = pool / "published_dfl_train.jsonl"
        row = json.loads(path.read_text())
        row["task"] = "DLD"
        manifest_path = pool / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["files"]["published_dfl_train"]["sha256"] = write(path, json.dumps(row) + "\n")
        write(manifest_path, manifest)
        with self.assertRaisesRegex(ValueError, "Non-DFL source"):
            audit(self.root)

    def test_weak_uniprot_candidate_is_reported_separately(self):
        weak = self.root / "data/processed/weak"
        row = {"source_id": "WEAK", "sequence": "G" * 80, "labels": "1" + "x" * 79,
               "record_type": "uniprot_explicit_flexible_linker_candidate",
               "training_eligibility": "candidate_not_DFL_gold_or_negative_adjudicated",
               "annotation": {"evidence_tier": "explicit_region_weak"}}
        path = weak / "candidates.jsonl"
        source_hash = write(path, json.dumps(row) + "\n")
        write(weak / "manifest.json", {"status": "explicit_flexible_linker_candidates_not_DFL_gold",
                                       "validated": 1, "output_sha256": {"candidates.jsonl": source_hash}})
        result = audit(self.root, weak_candidates=weak)
        self.assertEqual(result["weak_uniprot_candidates"]["counts"]["explicit_region_weak_total"], 1)
        self.assertEqual(result["combined_unique_surviving_sequences"], 3)


if __name__ == "__main__":
    unittest.main()
