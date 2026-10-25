"""Tests for Part 3 (`linkerforge.attachment`): position-controlled scorer and site-pair ranking."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.attachment import (  # noqa: E402
    diagnose_original_task_leakage,
    evaluate_attachment_scorer_benchmark,
    fit_two_tower_attachment_scorer,
    matched_offset_positions,
    rank_attachment_site_pairs,
    score_domain_attachment_sites,
)
from linkerforge.geometry.structure import read_chains  # noqa: E402

PATTERN = "ACDEFGHIKLMNPQRSTVWY" * 3


def _record(i: int, left: str, right: str) -> dict:
    return {
        "record_id": f"rec_{i}",
        "homology_cluster_id": f"cluster_{i % 5}",
        "linker_sequence": "GGGGSGGGGS",
        "linker_length": 10,
        "flanking_context": {"left_flank_sequence": left, "right_flank_sequence": right},
    }


class TestPart3AttachmentScorer(unittest.TestCase):
    def setUp(self) -> None:
        # Real composition signal: the linker-adjacent ends are Gly/Ser/Pro rich.
        core = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNAL"
        self.signal = [_record(i, core + "GGSGPGS", "GSGPGSG" + core) for i in range(30)]
        # No signal: each window is a palindrome, so both edges have identical composition.
        pal = PATTERN + PATTERN[::-1]
        self.null = [_record(i, pal, pal) for i in range(30)]

    def test_positional_feature_weights_are_zero(self) -> None:
        scorer = fit_two_tower_attachment_scorer(self.signal, epochs=40)
        self.assertEqual(len(scorer.weights), 12)
        self.assertEqual(scorer.weights[9], 0.0)
        self.assertEqual(scorer.weights[11], 0.0)

    def test_matched_offsets_are_mirror_images(self) -> None:
        pos, neg = matched_offset_positions(40, linker_side="c_terminal")
        self.assertEqual(pos, [39, 38, 37])
        self.assertEqual(neg, [0, 1, 2])
        pos, neg = matched_offset_positions(40, linker_side="n_terminal")
        self.assertEqual(pos, [0, 1, 2])
        self.assertEqual(neg, [39, 38, 37])

    def test_benchmark_is_at_chance_when_only_position_differs(self) -> None:
        res = evaluate_attachment_scorer_benchmark(self.null, n_folds=5, seed=3, bootstrap_reps=50)
        auroc = res["arms"]["sequence_context_logistic"]["auroc_pooled"]
        self.assertAlmostEqual(auroc, 0.5, delta=0.02)

    def test_benchmark_detects_real_composition_signal(self) -> None:
        res = evaluate_attachment_scorer_benchmark(self.signal, n_folds=5, seed=3, bootstrap_reps=50)
        self.assertGreater(res["arms"]["sequence_context_logistic"]["auroc_pooled"], 0.8)
        self.assertAlmostEqual(res["arms"]["random_control"]["auroc_pooled"], 0.5, delta=0.1)

    def test_original_design_is_solved_by_position_alone(self) -> None:
        diag = diagnose_original_task_leakage(self.null)
        self.assertGreater(diag["auroc_position_only_term_prox"], 0.95)

    def test_rank_attachment_site_pairs_on_real_pdb_fixture(self) -> None:
        scorer = fit_two_tower_attachment_scorer(self.signal, epochs=30)
        chains = read_chains(ROOT / "tests" / "fixtures" / "geometry" / "4e43_protease_dimer.pdb")
        dom_a, dom_b = chains["A"], chains["B"]
        seq_a, seq_b = "G" * len(dom_a.residues), "G" * len(dom_b.residues)
        ranked = rank_attachment_site_pairs(
            dom_a, seq_a, dom_b, seq_b, scorer,
            linker_type="flexible", requested_residues=15, topology="end_to_end_tandem",
            candidate_sites_per_domain=4, top_k=3,
        )
        self.assertEqual(len(ranked), 3)
        self.assertEqual(ranked[0].rank, 1)
        self.assertGreaterEqual(ranked[0].combined_site_pair_score, ranked[1].combined_site_pair_score)
        self.assertGreater(ranked[0].inter_port_ca_distance_A, 0.0)
        sites = score_domain_attachment_sites(seq_a, scorer, chain=dom_a, top_k=5)
        self.assertEqual(len(sites), 5)
        self.assertIsNotNone(sites[0].relative_sasa)


if __name__ == "__main__":
    unittest.main()
