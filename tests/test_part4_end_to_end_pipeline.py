"""Tests for Part 4 (`linkerforge.pipeline`): scorecard semantics, context extraction, evaluation."""

from __future__ import annotations

from pathlib import Path
import random
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.pipeline import (  # noqa: E402
    SYNTHETIC_NEGATIVE_TYPES,
    decoy_ranking_benchmark,
    evaluate_hard_rejection_scorecard,
    generate_synthetic_negatives,
    record_query_context,
)
from linkerforge.pipeline.fusion_demo import linker_ca_path, place_domain_b  # noqa: E402


class TestScorecard(unittest.TestCase):
    def test_rules_without_inputs_are_reported_not_evaluated(self) -> None:
        card = evaluate_hard_rejection_scorecard("GGGGSGGGGS")
        self.assertTrue(card.hard_rejection_pass)
        self.assertIn("physical_span_closure", card.not_evaluated)
        self.assertIn("anchor_accessibility", card.not_evaluated)

    def test_span_closure_uses_l_plus_one_bonds(self) -> None:
        # (L + 1) * 3.8 A: 7 residues reach 30.4 A, 6 residues reach 26.6 A.
        self.assertTrue(evaluate_hard_rejection_scorecard("G" * 7, required_path_length_A=30.0).physical_span_closure_pass)
        self.assertFalse(evaluate_hard_rejection_scorecard("G" * 6, required_path_length_A=30.0).physical_span_closure_pass)
        self.assertTrue(evaluate_hard_rejection_scorecard("", required_path_length_A=3.0).physical_span_closure_pass)
        self.assertFalse(evaluate_hard_rejection_scorecard("", required_path_length_A=10.0).hard_rejection_pass)

    def test_constructed_corruptions_trip_their_own_rules(self) -> None:
        negs = generate_synthetic_negatives("GGGGSGGGGSGGGGS", target_length=15)
        self.assertEqual(set(negs), set(SYNTHETIC_NEGATIVE_TYPES))
        self.assertFalse(evaluate_hard_rejection_scorecard(negs["all_hydrophobic_core"]).hydrophobic_stretch_pass)
        self.assertFalse(
            evaluate_hard_rejection_scorecard(negs["protease_site_buried_in_rigid_helix"]).composite_rigid_cleavable_pass
        )

    def test_buried_port_is_rejected_only_when_accessibility_is_supplied(self) -> None:
        self.assertFalse(
            evaluate_hard_rejection_scorecard("GGGGS", port_a_sasa=0.05, port_b_sasa=0.5).anchor_accessibility_pass
        )


class TestContext(unittest.TestCase):
    def test_no_geometry_is_not_invented(self) -> None:
        ctx = record_query_context({"linker_sequence": "GGGGS", "linker_length": 5})
        self.assertIsNone(ctx["path_A"])
        self.assertIsNone(ctx["sasa_a"])

    def test_real_geometry_is_used_when_recorded(self) -> None:
        ctx = record_query_context(
            {
                "linker_sequence": "GGGGS",
                "geometry": {"inter_port_ca_distance_A": 12.5, "port_a_accessibility_proxy": 0.7, "port_b_accessibility_proxy": 0.4},
                "biophysical_properties": {"curated_rigidity_label": "rigid"},
            }
        )
        self.assertEqual(ctx["path_A"], 12.5)
        self.assertEqual(ctx["rigidity"], "rigid")


class TestDecoyBenchmark(unittest.TestCase):
    def test_structure_and_chance_levels(self) -> None:
        rng = random.Random(1)
        recs = []
        for i in range(80):
            seq = "".join(rng.choice("GSAPEK") for _ in range(10 + (i % 3)))
            recs.append(
                {
                    "record_id": f"r{i}", "homology_cluster_id": f"c{i}",
                    "linker_sequence": seq, "linker_length": len(seq),
                    "biophysical_properties": {"rigidity_predicted_label": "flexible"},
                    "flanking_context": {"left_flank_sequence": "A" * 30, "right_flank_sequence": "K" * 30},
                }
            )
        res = decoy_ranking_benchmark(recs[:40], recs[40:], n_decoys=9, min_decoys=5)
        self.assertGreater(res["n_queries"], 10)
        self.assertAlmostEqual(res["chance"]["hit_at_1"], 0.1, delta=0.02)
        for arm in ("random", "training_composition_prior", "flank_conditioned_nn", "linkerforge_ranker"):
            self.assertIn("mrr", res["arms"][arm])
        # a random scorer must sit near chance
        self.assertLess(abs(res["arms"]["random"]["hit_at_1"]["mean"] - 0.1), 0.2)


class TestPlacement(unittest.TestCase):
    def test_placed_domain_keeps_minimum_gap(self) -> None:
        rng = np.random.default_rng(0)
        a = rng.normal(size=(300, 3)) * 8.0
        b = rng.normal(size=(300, 3)) * 8.0
        placement = place_domain_b(a, b, a[0], b[0], n_rotations=12, seed=3)
        placed = placement["transform"](b)
        d = np.sqrt(((a[:, None, :] - placed[None, :, :]) ** 2).sum(-1)).min()
        self.assertGreaterEqual(d, 3.99)

    def test_linker_trace_reports_reachability(self) -> None:
        pa, pb = np.zeros(3), np.array([30.0, 0.0, 0.0])
        _, short = linker_ca_path(pa, pb, 3, np.array([0.0, 1.0, 0.0]))
        self.assertFalse(short["reachable_at_3p8A_per_bond"])
        pts, longer = linker_ca_path(pa, pb, 12, np.array([0.0, 1.0, 0.0]))
        self.assertTrue(longer["reachable_at_3p8A_per_bond"])
        self.assertEqual(len(pts), 12)
        self.assertLessEqual(longer["ca_ca_spacing_max_A"], 3.85)


if __name__ == "__main__":
    unittest.main()
