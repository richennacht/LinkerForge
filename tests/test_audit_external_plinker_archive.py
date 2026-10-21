from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile

from scripts.audit_external_plinker_archive import build_audit


class ExternalPLinkerArchiveAuditTest(unittest.TestCase):
    def test_audits_without_extracting_or_executing_archive(self) -> None:
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            archive_path = root / "Linker.zip"
            silver = {
                "protein_accession": "P1", "protein_length": 11,
                "domain_A_acc": "PF00001", "domain_A_start": 1, "domain_A_end": 3,
                "domain_B_acc": "PF00001", "domain_B_start": 9, "domain_B_end": 11,
                "linker_start": 4, "linker_end": 8, "linker_length": 5,
                "sequences": {"domain_A": "AAA", "linker": "GGGGS", "domain_B": "CCC", "full": "AAAGGGGSCCC"},
                "confidence": "silver",
            }
            gold = {
                "pdb_id": "1ABC",
                "domains": {"domain_A": {"end": 3}, "domain_B": {"start": 9}},
                "linker": {"start": 4, "end": 8},
            }
            with ZipFile(archive_path, "w") as archive:
                archive.writestr("Linker - Copy/data/processed/dataset_v2.jsonl", json.dumps(silver) + "\n")
                archive.writestr("Linker - Copy/data/processed/dataset_v1.jsonl", json.dumps(gold) + "\n")
                archive.writestr("Linker - Copy/data/processed/mined_linkers.jsonl", "{}\n")
                archive.writestr("Linker - Copy/docs/data_card.md", "Total 7,903\n| Gold | 11 |\n")
                archive.writestr("Linker - Copy/.venv/ignored.py", "raise RuntimeError('must not run')")
                archive.writestr("Linker - Copy/.git/config", "untrusted")
            existing = root / "existing.jsonl"
            existing.write_text(json.dumps({"protein_id": "P1", "linker_sequence": "GGGGS"}) + "\n", encoding="utf-8")

            result = build_audit(archive_path, existing)

            self.assertEqual("quarantined_not_training_authorized", result["status"])
            self.assertEqual(1, result["silver_profile"]["within_linkerforge_5_to_30_residues"])
            self.assertEqual(1, result["overlap_with_existing_27766_corpus"]["exact_protein_identifier_overlap"])
            self.assertEqual(0, result["gold_profile"]["declared_boundary_mismatch_records"])
            self.assertFalse(result["fitness"]["compatibility_ranker_training"])


if __name__ == "__main__":
    unittest.main()
