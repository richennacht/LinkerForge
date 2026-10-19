from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from scripts.preflight_dfl_family_training_v2 import (
    RULE_V2,
    v2_alignment_metrics,
    v2_homolog_hit,
    verify,
)


def _sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _sha256_path(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


class PreflightDflFamilyTrainingV2Test(unittest.TestCase):
    def test_v2_alignment_metrics_bounds_coverage_to_one(self) -> None:
        left = "GS" * 35
        right = "G" * 70
        metrics = v2_alignment_metrics(left, right)
        self.assertLessEqual(float(metrics["left_residue_coverage"]), 1.0)
        self.assertLessEqual(float(metrics["right_residue_coverage"]), 1.0)
        self.assertGreaterEqual(float(metrics["left_residue_coverage"]), 0.0)
        self.assertGreaterEqual(float(metrics["right_residue_coverage"]), 0.0)

    def test_v2_homolog_hit_detects_true_homolog_and_rejects_unrelated(self) -> None:
        seq_a = "ACDEFGHIKLMNPQRSTVWY" * 4
        seq_b = "ACDEFGHIKLMNPQRSTVWY" * 4
        seq_unrelated = "MKTIIALSYIFCLVFADYKDDDDK" * 3
        self.assertTrue(v2_homolog_hit(seq_a, seq_b))
        self.assertFalse(v2_homolog_hit(seq_a, seq_unrelated))

    def _write_valid_v2_fixture(self, root: Path) -> Path:
        reg_dir = root / "registries"
        reg_dir.mkdir(parents=True, exist_ok=True)

        caid_fasta_text = ">DP09999\n" + ("MKTIIALSYIFCLVFADYKDDDDK" * 3) + "\n"
        (reg_dir / "caid3.fasta").write_text(caid_fasta_text, encoding="utf-8", newline="\n")
        caid_reg = {
            "schema_version": "1.0.0",
            "name": "caid3",
            "status": "sequence_only_protected",
            "label_seal_state": "retrospective_spent_exclusion_only",
            "source_release": "CAID3_linker_reference",
            "fasta_path": "caid3.fasta",
            "fasta_sha256_canonical_lf": _sha256_text(caid_fasta_text),
            "ids": ["DP09999"],
        }
        caid_reg_path = reg_dir / "caid3.registry.json"
        caid_reg_path.write_text(json.dumps(caid_reg, indent=2) + "\n", encoding="utf-8", newline="\n")

        dp94_fasta_text = ">DP08888\n" + ("WCHMYWVCHMYWVCHMYWVCHMYW" * 3) + "\n"
        (reg_dir / "dp94.fasta").write_text(dp94_fasta_text, encoding="utf-8", newline="\n")
        dp94_reg = {
            "schema_version": "1.0.0",
            "name": "dp94",
            "status": "sequence_only_labels_sealed",
            "label_seal_state": "unopened",
            "source_release": "DisoFLAG_DP94_sequence_only",
            "fasta_path": "dp94.fasta",
            "fasta_sha256_canonical_lf": _sha256_text(dp94_fasta_text),
            "ids": ["DP08888"],
        }
        dp94_reg_path = reg_dir / "dp94.registry.json"
        dp94_reg_path.write_text(json.dumps(dp94_reg, indent=2) + "\n", encoding="utf-8", newline="\n")

        train_rows = [
            {
                "protein_id": "DP00001",
                "sequence": "AGGGSAGGGSAGGGSAGGGSAGGGSAGGGSAGGGSAGGGSAGGGSAGGGSAGGGSAGGGSAAAAA",
                "labels": "1" * 20 + "x" * 45,
                "source": "DisProt_2026_06",
                "task": "DFL_curated_positive_and_experimental_order",
                "family_component_id": "comp_train_pos",
                "provenance": {"source_id": "DP00001", "input_sha256": "a" * 64},
            },
            {
                "protein_id": "DP00002",
                "sequence": "DEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKRCCCCC",
                "labels": "0" * 25 + "x" * 40,
                "source": "DisProt_2026_06",
                "task": "DFL_curated_positive_and_experimental_order",
                "family_component_id": "comp_train_neg",
                "provenance": {"source_id": "DP00002", "input_sha256": "a" * 64},
            },
        ]
        dev_rows = [
            {
                "protein_id": "DP01001",
                "sequence": "LVIFLVIFLVIFLVIFLVIFLVIFLVIFLVIFLVIFLVIFLVIFLVIFLVIFLVIFLVIFDDDDD",
                "labels": "1" * 20 + "x" * 45,
                "source": "DisProt_2026_06",
                "task": "DFL_curated_positive_and_experimental_order",
                "family_component_id": "comp_dev_pos",
                "provenance": {"source_id": "DP01001", "input_sha256": "a" * 64},
            },
            {
                "protein_id": "DP01002",
                "sequence": "NQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQNQEEEEE",
                "labels": "0" * 25 + "x" * 40,
                "source": "DisProt_2026_06",
                "task": "DFL_curated_positive_and_experimental_order",
                "family_component_id": "comp_dev_neg",
                "provenance": {"source_id": "DP01002", "input_sha256": "a" * 64},
            },
        ]
        aux_rows = [
            {
                "protein_id": "DP02001",
                "sequence": "PTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTPTGGGGG",
                "labels": "0" * 30 + "x" * 35,
                "source": "DisProt_2026_06_order_panel",
                "task": "DFL_experimental_order_negative_panel",
                "family_component_id": "comp_aux_neg",
                "provenance": {"source_id": "DP02001", "input_sha256": "a" * 64},
            }
        ]

        membership = {
            "DP09999": "comp_bench_caid3",
            "DP08888": "comp_bench_dp94",
            "DP00001": "comp_train_pos",
            "DP00002": "comp_train_neg",
            "DP01001": "comp_dev_pos",
            "DP01002": "comp_dev_neg",
            "DP02001": "comp_aux_neg",
        }
        family_doc = {
            "schema_version": "2.0.0",
            "status": "family_clustered_v2",
            "method": "transitive_Pfam_Gene3D_UniRef50_plus_RULE_V2",
            "version": "2.0.0",
            "source_sha256": "a" * 64,
            "membership": membership,
        }
        family_path = root / "family_manifest.json"
        family_path.write_text(json.dumps(family_doc, indent=2) + "\n", encoding="utf-8", newline="\n")

        files_info = {}
        counts = {}
        mask_counts = {}
        for role, rows in [
            ("train", train_rows),
            ("development", dev_rows),
            ("auxiliary", aux_rows),
        ]:
            rel = f"{role}.jsonl"
            p = root / rel
            p.write_text("".join(json.dumps(r) + "\n" for r in rows), encoding="utf-8", newline="\n")
            files_info[role] = {"path": rel, "records": len(rows), "sha256": _sha256_path(p)}
            counts[role] = len(rows)
            mask_counts[role] = {ch: sum(r["labels"].count(ch) for r in rows) for ch in "10x"}

        manifest = {
            "schema_version": "2.0.0",
            "status": "preregistered_v2_family_disjoint_split",
            "preregistration_doc": "literature/DFL_FIREWALL_V2_PREREGISTRATION_2026-10-18.md",
            "similarity_rule": RULE_V2,
            "shuffle_null_verification": {
                "seeds": [173, 271, 811],
                "sample_size_each_side": 40,
                "total_pairs": 4800,
                "observed_hits": 11,
                "residue_shuffle_hits": 0,
                "block12_shuffle_hits": 1,
                "residue_shuffle_rate": 0.0,
                "block12_shuffle_rate": 0.00020833333333333335,
            },
            "source_policies": {
                "DisProt_2026_06": {
                    "task": "DFL_curated_positive_and_experimental_order",
                    "zero_semantics": "explicit_negative_curated",
                },
                "DisProt_2026_06_order_panel": {
                    "task": "DFL_experimental_order_negative_panel",
                    "zero_semantics": "explicit_negative_curated",
                },
            },
            "benchmark_registries": [
                {"name": "caid3", "path": "registries/caid3.registry.json", "sha256": _sha256_path(caid_reg_path)},
                {"name": "dp94", "path": "registries/dp94.registry.json", "sha256": _sha256_path(dp94_reg_path)},
            ],
            "family_manifest": {"path": "family_manifest.json", "sha256": _sha256_path(family_path)},
            "files": files_info,
            "counts": counts,
            "mask_counts": mask_counts,
            "limitations": ["Synthetic test fixture for V2 firewall verification."],
        }
        manifest_path = root / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")
        return manifest_path

    def test_verify_passes_valid_v2_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest_path = self._write_valid_v2_fixture(Path(tmp))
            result = verify(manifest_path)
            self.assertEqual(result["status"], "pass_v2_preregistered_sequence_and_family_firewall")

    def test_verify_rejects_excessive_shuffle_hit_rate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest_path = self._write_valid_v2_fixture(Path(tmp))
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            payload["shuffle_null_verification"]["residue_shuffle_rate"] = 0.012
            manifest_path.write_text(json.dumps(payload), encoding="utf-8", newline="\n")
            with self.assertRaises(Exception):
                verify(manifest_path)

    def test_verify_rejects_benchmark_family_collision(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest_path = self._write_valid_v2_fixture(root)
            family_path = root / "family_manifest.json"
            fam_doc = json.loads(family_path.read_text(encoding="utf-8"))
            fam_doc["membership"]["DP08888"] = "comp_train_pos"
            family_path.write_text(json.dumps(fam_doc, indent=2) + "\n", encoding="utf-8", newline="\n")
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            payload["family_manifest"]["sha256"] = _sha256_path(family_path)
            manifest_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
            with self.assertRaisesRegex(ValueError, "Benchmark-family collision"):
                verify(manifest_path)


if __name__ == "__main__":
    unittest.main()
