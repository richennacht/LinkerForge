"""Unit tests for B2/B3 two-stage domain-gated conformal mining and dataset_v1 loader."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.dataset_v1 import load_dataset_manifest, load_linker_records
from linkerforge.detector.conformal_mining import (
    calibrate_conformal_threshold,
    call_two_stage_interdomain_spans,
)
from linkerforge.detector.domain_gate import DomainMatch


class TestConformalMining(unittest.TestCase):
    def test_calibrate_conformal_threshold_attains_low_fdr(self) -> None:
        # 40 clean true positives at high scores (>= 0.80) and 50 negatives at low scores (<= 0.40)
        labels = [1] * 40 + [0] * 50
        scores = [0.85 + 0.002 * i for i in range(40)] + [0.05 + 0.005 * i for i in range(50)]
        res = calibrate_conformal_threshold(
            labels,
            scores,
            target_fdr=0.10,
            alpha=0.05,
            min_true_positives=20,
        )
        self.assertTrue(res.attained)
        self.assertIsNotNone(res.calibrated_threshold)
        self.assertEqual(res.false_positives, 0)
        self.assertEqual(res.true_positives, 40)
        self.assertLessEqual(res.ucb_fdr, 0.10)
        self.assertAlmostEqual(res.empirical_recall, 1.0)

    def test_calibrate_conformal_threshold_abstains_when_unattainable(self) -> None:
        # Inseparable scores where positives have lower scores than negatives
        labels = [1] * 10 + [0] * 20
        scores = [0.20] * 10 + [0.90] * 20
        res = calibrate_conformal_threshold(
            labels,
            scores,
            target_fdr=0.10,
            alpha=0.05,
            min_true_positives=5,
        )
        self.assertFalse(res.attained)
        self.assertIsNone(res.calibrated_threshold)
        self.assertAlmostEqual(res.abstention_rate, 1.0)

    def test_two_stage_caller_rejects_single_domain_and_homo_repeat(self) -> None:
        seq = "A" * 100
        high_scores = [0.95] * 100
        # Single domain -> Stage 1 abstains even if all CNN scores are 0.95
        single = [DomainMatch(10, 80, "PF00001", "domain")]
        out_single = call_two_stage_interdomain_spans(seq, single, high_scores)
        self.assertFalse(out_single["stage1_passed"])
        self.assertEqual(out_single["called_high_confidence_spans"], 0)

        # Homo-repeat domain pair -> Stage 1 abstains
        homo = [
            DomainMatch(1, 30, "PF00001", "domain"),
            DomainMatch(45, 90, "PF00001", "domain"),
        ]
        out_homo = call_two_stage_interdomain_spans(seq, homo, high_scores)
        self.assertFalse(out_homo["stage1_passed"])
        self.assertEqual(out_homo["called_high_confidence_spans"], 0)

    def test_two_stage_caller_stratifies_hetero_domain_spans(self) -> None:
        seq = "A" * 120
        domains = [
            DomainMatch(1, 30, "PF00001", "domain"),
            DomainMatch(45, 75, "PF00002", "domain"),
            DomainMatch(90, 120, "PF00003", "domain"),
        ]
        # Gap 1: 31..44 (14 aa) -> high score 0.85
        # Gap 2: 76..89 (14 aa) -> provisional score 0.60
        scores = [0.10] * 120
        for i in range(30, 44):
            scores[i] = 0.85
        for i in range(75, 89):
            scores[i] = 0.60
        out = call_two_stage_interdomain_spans(
            seq,
            domains,
            scores,
            high_confidence_threshold=0.735198,
            provisional_threshold=0.50,
        )
        self.assertTrue(out["stage1_passed"])
        self.assertEqual(out["called_high_confidence_spans"], 1)
        self.assertEqual(out["called_provisional_spans"], 1)
        self.assertEqual(out["candidate_spans"][0]["calling_tier"], "mined_high_confidence_dfl_call")
        self.assertEqual(out["candidate_spans"][1]["calling_tier"], "mined_provisional_dfl_candidate")

    def test_dataset_v1_loader_verifies_sha256_and_filters(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            rec_path = tmp_dir / "all_records.jsonl"
            row = {
                "record_id": "cath:test",
                "authority_tier": "natural_cath_pdb_verified",
                "linker_sequence": "GGGGSGGGGS",
                "linker_length": 10,
                "geometry": {"physical_3p8A_feasible": True},
            }
            rec_path.write_text(json.dumps(row) + "\n", encoding="utf-8")
            import hashlib

            digest = hashlib.sha256(rec_path.read_bytes()).hexdigest()
            manifest = {
                "files": {
                    "all_records_jsonl": {
                        "path": "all_records.jsonl",
                        "records": 1,
                        "sha256": digest,
                    }
                }
            }
            (tmp_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
            loaded = load_linker_records(
                tmp_dir,
                tiers=["natural_cath_pdb_verified"],
                only_feasible_3p8A=True,
            )
            self.assertEqual(len(loaded), 1)
            self.assertEqual(loaded[0]["record_id"], "cath:test")

    def test_real_dataset_v1_manifest_if_present(self) -> None:
        from linkerforge.dataset_v1 import DEFAULT_DATASET_DIR

        if not (DEFAULT_DATASET_DIR / "manifest.json").exists():
            self.skipTest("mined_and_natural_linkers_v1 not built in this checkout")
        manifest = load_dataset_manifest(DEFAULT_DATASET_DIR)
        self.assertEqual(manifest["dataset_id"], "mined_and_natural_linkers_v1")
        natural = load_linker_records(
            DEFAULT_DATASET_DIR,
            tiers=["natural_cath_pdb_verified", "natural_disprot_curated_dfl"],
        )
        self.assertEqual(len(natural), 148)
        mined_hi = load_linker_records(
            DEFAULT_DATASET_DIR,
            tiers=["mined_high_confidence_dfl_call"],
        )
        self.assertEqual(len(mined_hi), 6556)
        disprot_hi = [
            r for r in mined_hi if r.get("mining_source_corpus") == "DisProt_2026_06_Pfam_Gaps"
        ]
        self.assertEqual(len(disprot_hi), 88)


if __name__ == "__main__":
    unittest.main()
