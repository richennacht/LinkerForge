"""Fail-closed checks for the CAID3/DP94 sequence and V2 family firewall (`B1`)."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from build_caid3_protected_dfl_splits import OUTPUT, REFERENCE, ROOT, match, reference_sequences
from check_detector_autoresearch import inspect as inspect_autoresearch
from preflight_caid3_protected_training import V2_DEFAULT_OUTPUT, verify
import run_caid3_protected_dfl_trial as trial_launcher


class Caid3BlacklistTests(unittest.TestCase):
    def test_exact_identifier_is_excluded(self):
        identifier, sequence = reference_sequences(REFERENCE)[0]
        self.assertEqual(match({"protein_id": identifier, "sequence": sequence},
                               reference_sequences(REFERENCE))["reason"],
                         "exact_id_sequence_or_containment")

    def test_exact_sequence_with_new_identifier_is_excluded(self):
        _, sequence = reference_sequences(REFERENCE)[0]
        self.assertIsNotNone(match({"protein_id": "synthetic_test_id", "sequence": sequence},
                                   reference_sequences(REFERENCE)))

    def test_old_versioned_split_fails_prospective_holdout_gate(self):
        with self.assertRaisesRegex(ValueError, "Prospective DP94 benchmark sequence overlap"):
            verify(OUTPUT / "manifest.json")

    def test_v2_family_disjoint_split_passes_preflight_and_launcher(self):
        v2_manifest = V2_DEFAULT_OUTPUT / "manifest.json"
        result = verify(v2_manifest)
        self.assertEqual(result["status"], "pass_v2_preregistered_sequence_and_family_firewall")
        self.assertEqual(result["counts"], {"train": 60, "development": 28, "auxiliary": 9})

        status = inspect_autoresearch(ROOT / "reports/dfl_b1_v2_family_round_2026-10-18.json")
        self.assertEqual(status["action"], "launch_or_audit")
        self.assertEqual(status["preflight_status"], "pass_v2_preregistered_sequence_and_family_firewall")

        train, dev = trial_launcher.preflight()
        self.assertEqual(len(train), 60)
        self.assertEqual(len(dev), 28)


if __name__ == "__main__":
    unittest.main()
