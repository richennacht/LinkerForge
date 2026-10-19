import json
import tempfile
import unittest
from pathlib import Path

from scripts.summarize_dfl_mmseqs_audit import sha, summarize


class MMseqsSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        records = [
            {"record_id": "C0000", "role": "candidate_positive_not_training_cleared"},
            {"record_id": "C0001", "role": "candidate_positive_not_training_cleared"},
            {"record_id": "DP94_000", "role": "sealed_benchmark_sequence_only", "registry": "DP94"},
        ]
        self.records = self.root / "records.jsonl"
        self.records.write_text("".join(json.dumps(row) + "\n" for row in records), encoding="utf-8")
        self.fasta = self.root / "sequences.fasta"
        self.fasta.write_text(">C0000\nACDE\n>C0001\nFGHI\n>DP94_000\nKLMN\n", encoding="utf-8")
        manifest = {"status": "sequence_only_family_audit_not_split", "candidate_sequences": 2,
                    "benchmark_sequences": 1,
                    "files_sha256": {self.records.name: sha(self.records), self.fasta.name: sha(self.fasta)}}
        (self.root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        self.tsv = self.root / "clusters.tsv"
        self.tsv.write_text("C0000\tC0000\nC0000\tDP94_000\nC0001\tC0001\n", encoding="utf-8")

    def test_counts_touching_and_candidate_only_clusters(self):
        result = summarize(self.root, self.tsv, sha(self.tsv))
        self.assertEqual(result["counts"]["candidate_sequences_in_benchmark_clusters"], 1)
        self.assertEqual(result["counts"]["candidate_sequences_outside_benchmark_clusters"], 1)

    def test_rejects_hash_change(self):
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            summarize(self.root, self.tsv, "0" * 64)

    def test_rejects_missing_member(self):
        self.tsv.write_text("C0000\tC0000\nC0001\tC0001\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "omits"):
            summarize(self.root, self.tsv, sha(self.tsv))


if __name__ == "__main__":
    unittest.main()
