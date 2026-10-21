from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from jsonschema import Draft202012Validator

from scripts.build_linker_curation_review_packet import build_packet


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(index: int, tier: str, score: float, source: str) -> dict:
    return {
        "record_id": f"record-{index}",
        "protein_id": f"protein-{index}",
        "authority_tier": tier,
        "linker_sequence": "GGGGS",
        "linker_start_1based": 11,
        "linker_end_1based": 15,
        "linker_length": 5,
        "homology_cluster_id": f"cluster-{index}",
        "evidence_grade": "test-evidence",
        "mining_source_corpus": source,
        "benchmark_firewall_verified": True,
        "classifier_scores": {"cnn_mean_score": score},
        "flanking_context": {
            "left_domain_id": "PF00001",
            "right_domain_id": "PF00002",
            "left_flank_sequence": "AAAAA",
            "right_flank_sequence": "CCCCC",
        },
        "geometry": None,
        "biophysical_properties": {},
    }


class BuildLinkerCurationReviewPacketTest(unittest.TestCase):
    def _fixture(self, root: Path) -> tuple[Path, Path]:
        dataset = root / "dataset"
        dataset.mkdir()
        records = [
            _record(1, "mined_high_confidence_dfl_call", 0.80, "source-a"),
            _record(2, "mined_high_confidence_dfl_call", 0.78, "source-b"),
            _record(3, "mined_provisional_dfl_candidate", 0.61, "source-a"),
            _record(4, "mined_provisional_dfl_candidate", 0.55, "source-b"),
            _record(5, "natural_cath_pdb_verified", 0.70, "natural"),
        ]
        all_records = dataset / "all_records.jsonl"
        all_records.write_text(
            "".join(json.dumps(row) + "\n" for row in records), encoding="utf-8", newline="\n"
        )
        manifest = {
            "dataset_id": "mined_and_natural_linkers_v1",
            "files": {
                "all_records_jsonl": {
                    "path": "all_records.jsonl",
                    "records": len(records),
                    "sha256": _sha256(all_records),
                }
            },
        }
        (dataset / "manifest.json").write_text(
            json.dumps(manifest) + "\n", encoding="utf-8", newline="\n"
        )
        config = root / "config.json"
        config.write_text(
            json.dumps(
                {
                    "source_dataset_id": "mined_and_natural_linkers_v1",
                    "selection_seed": "test-seed",
                    "minimum_total_records": 1,
                    "maximum_total_records": 10,
                    "tier_quotas": {
                        "mined_high_confidence_dfl_call": 1,
                        "mined_provisional_dfl_candidate": 1,
                        "natural_cath_pdb_verified": 1,
                    },
                    "length_bins": [10, 20, 30, 45, 60],
                    "claim_boundary": "test boundary",
                }
            )
            + "\n",
            encoding="utf-8",
            newline="\n",
        )
        return dataset, config

    def test_builds_deterministic_independent_review_packets(self) -> None:
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            dataset, config = self._fixture(root)
            first = root / "first"
            second = root / "second"
            first_manifest = build_packet(dataset, config, first)
            second_manifest = build_packet(dataset, config, second)

            self.assertEqual(3, first_manifest["selection"]["selected_total"])
            self.assertEqual(
                (first / "selection_index.jsonl").read_bytes(),
                (second / "selection_index.jsonl").read_bytes(),
            )
            curator_a = [json.loads(line) for line in (first / "curator_a.jsonl").read_text().splitlines()]
            curator_b = [json.loads(line) for line in (first / "curator_b.jsonl").read_text().splitlines()]
            self.assertEqual(
                [row["record_id"] for row in curator_a],
                [row["record_id"] for row in curator_b],
            )
            self.assertTrue(all(row["compatibility_outcome"] == "not_assessed" for row in curator_a))
            self.assertTrue(all(row["review"]["status"] == "pending" for row in curator_a))
            self.assertTrue(all(row["review_id"].startswith("curator_a:") for row in curator_a))

            schema_path = Path(__file__).resolve().parents[1] / "schemas" / "linker-curation-review-record.v1.schema.json"
            validator = Draft202012Validator(json.loads(schema_path.read_text(encoding="utf-8")))
            for row in curator_a + curator_b:
                validator.validate(row)

    def test_refuses_to_overwrite_existing_packet(self) -> None:
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            dataset, config = self._fixture(root)
            output = root / "packet"
            build_packet(dataset, config, output)
            with self.assertRaises(FileExistsError):
                build_packet(dataset, config, output)


if __name__ == "__main__":
    unittest.main()
