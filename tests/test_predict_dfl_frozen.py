"""No benchmark labels or real model scores are used in these contract tests."""

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/predict_dfl_frozen.py"
spec = importlib.util.spec_from_file_location("predict_dfl_frozen", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class FrozenDFLInferenceTests(unittest.TestCase):
    def test_fasta_rejects_multiple_entries_and_noncanonical_residues(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "input.fasta"
            fasta.write_text(">p\nACD\n>q\nEFG\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                module.read_fasta(fasta)
            fasta.write_text(">p\nACDX\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                module.read_fasta(fasta)
            fasta.write_text(">p note\nacde\n", encoding="utf-8")
            self.assertEqual(module.read_fasta(fasta), ("p", "ACDE"))

    def test_abstain_is_default_and_explicit_threshold_is_provisional(self):
        common = dict(identifier="synthetic", sequence="ACDE", scores=[0.1, 0.8, 0.9, 0.2],
                      checkpoint_sha256="0" * 64, model_id="test-encoder", revision="test-rev")
        abstain = module.make_report(**common, threshold=None)
        self.assertTrue(abstain["abstained_from_span_calls"])
        self.assertEqual(abstain["candidate_spans"], [])
        self.assertEqual(abstain["claim_status"], "provisional_inference_not_certified_for_mining")
        called = module.make_report(**common, threshold=0.5)
        self.assertFalse(called["abstained_from_span_calls"])
        self.assertEqual([(r["start"], r["end"]) for r in called["candidate_spans"]], [(2, 3)])
        with self.assertRaises(ValueError):
            module.make_report(**{**common, "scores": [0.1, float("nan"), 0.9, 0.2]}, threshold=None)

    def test_frozen_selection_hashes_and_path_gate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = root / "artifacts" / "run"
            run.mkdir(parents=True)
            checkpoint = run / "best.pt"
            checkpoint.write_bytes(b"synthetic checkpoint bytes; never loaded")
            manifest = run / "manifest.json"
            manifest.write_text(json.dumps({"status": "completed", "arch": "samlike",
                                            "blacklist_status": "passed_before_training",
                                            "model_id": "facebook/esm2_t12_35M_UR50D"}), encoding="utf-8")
            selection = root / "selection.json"
            record = {"frozen_before_caid3": True,
                      "selection_basis": "filtered_development_only", "arch": "samlike",
                      "variant": None, "checkpoint": "artifacts/run/best.pt",
                      "checkpoint_sha256": module.sha256(checkpoint),
                      "checkpoint_manifest_sha256": module.sha256(manifest)}
            selection.write_text(json.dumps(record), encoding="utf-8")
            actual, _, _ = module.verified_checkpoint(selection, root)
            self.assertEqual(actual, checkpoint.resolve())
            checkpoint.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                module.verified_checkpoint(selection, root)
            checkpoint.write_bytes(b"synthetic checkpoint bytes; never loaded")
            record["checkpoint"] = "../outside.pt"
            selection.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "inside artifact root"):
                module.verified_checkpoint(selection, root)


if __name__ == "__main__":
    unittest.main()
