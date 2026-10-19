"""Unit and integration tests for Part 3 (`linkerforge.attachment`)."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.attachment import (  # noqa: E402
    TwoTowerAttachmentScorer,
    evaluate_attachment_scorer_benchmark,
    fit_two_tower_attachment_scorer,
    rank_attachment_site_pairs,
    score_domain_attachment_sites,
)
from linkerforge.geometry.structure import read_chains  # noqa: E402


class TestPart3AttachmentScorer(unittest.TestCase):
    def setUp(self) -> None:
        self.records = [
            {
                "record_id": f"rec_{i}",
                "homology_cluster_id": f"cluster_{i % 5}",
                "linker_sequence": "GGGGSGGGGS",
                "linker_length": 10,
                "flanking_context": {
                    "left_flank_sequence": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYRGGSGPGS",
                    "right_flank_sequence": "GSGPGSGMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
                },
                "geometry": {
                    "port_a_accessibility_proxy": 0.48,
                    "port_b_accessibility_proxy": 0.51,
                },
            }
            for i in range(15)
        ]

    def test_fit_two_tower_scorer_keeps_tick_boxes_out_of_encoder(self) -> None:
        scorer = fit_two_tower_attachment_scorer(self.records, epochs=40)
        self.assertIsInstance(scorer, TwoTowerAttachmentScorer)
        self.assertEqual(len(scorer.weights), 12)
        seq = self.records[0]["flanking_context"]["left_flank_sequence"]
        c_term_score = scorer.score_sequence_position(seq, len(seq) - 1, attachment_mode="c_terminal")
        core_score = scorer.score_sequence_position(seq, 25, attachment_mode="c_terminal")
        self.assertGreater(c_term_score, core_score)

    def test_rank_attachment_site_pairs_on_real_pdb_fixture(self) -> None:
        scorer = fit_two_tower_attachment_scorer(self.records, epochs=30)
        pdb_path = ROOT / "tests" / "fixtures" / "geometry" / "4e43_protease_dimer.pdb"
        chains = read_chains(pdb_path)
        dom_a = chains["A"]
        dom_b = chains["B"]
        seq_a = "G" * len(dom_a.residues)
        seq_b = "G" * len(dom_b.residues)

        ranked = rank_attachment_site_pairs(
            dom_a,
            seq_a,
            dom_b,
            seq_b,
            scorer,
            linker_type="flexible",
            requested_residues=15,
            topology="end_to_end_tandem",
            candidate_sites_per_domain=4,
            top_k=3,
        )
        self.assertEqual(len(ranked), 3)
        self.assertEqual(ranked[0].rank, 1)
        self.assertGreaterEqual(ranked[0].combined_site_pair_score, ranked[1].combined_site_pair_score)
        self.assertGreater(ranked[0].inter_port_ca_distance_A, 0.0)

        sites = score_domain_attachment_sites(seq_a, scorer, chain=dom_a, top_k=5)
        self.assertEqual(len(sites), 5)
        self.assertIsNotNone(sites[0].relative_sasa)

    def test_evaluate_attachment_scorer_benchmark_cluster_disjoint(self) -> None:
        res = evaluate_attachment_scorer_benchmark(self.records, n_folds=5, seed=17)
        arms = res["arms"]
        self.assertIn("two_tower_structure_plus_sequence", arms)
        self.assertIn("random_control", arms)
        self.assertGreater(
            arms["two_tower_structure_plus_sequence"]["auroc"],
            arms["random_control"]["auroc"],
        )


if __name__ == "__main__":
    unittest.main()
