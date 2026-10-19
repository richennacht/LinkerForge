"""Unit and integration tests for Part 4 (`linkerforge.pipeline`)."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.pipeline import (  # noqa: E402
    AuthoritySequenceIndex,
    BASELINE_NAMES,
    SYNTHETIC_NEGATIVE_TYPES,
    build_fusion_3d_pdb_and_pymol_script,
    evaluate_hard_rejection_scorecard,
    evaluate_six_baselines_benchmark,
    fit_pairwise_compatibility_ranker,
    generate_synthetic_negatives,
    run_sfgfp_mcherry_fusion_demo,
)


class TestPart4EndToEndPipeline(unittest.TestCase):
    def setUp(self) -> None:
        self.records = [
            {
                "record_id": f"rec_{i}",
                "homology_cluster_id": f"cluster_{i % 5}",
                "source_dataset": "unit_test",
                "evidence_tier": "natural_cath_pdb_verified",
                "linker_sequence": "GGGGSGGGGSGGGGS" if i % 2 == 0 else "GSAGSAAGSGEF",
                "linker_length": 15 if i % 2 == 0 else 12,
                "geometry": {
                    "estimated_end_to_end_ca_distance_A": 26.0,
                    "port_a_accessibility_proxy": 0.55,
                    "port_b_accessibility_proxy": 0.52,
                },
                "flanking_context": {
                    "left_flank_sequence": "KRDHMVLLEFVTAAGITHGMDELYK",
                    "right_flank_sequence": "MVSKGEEDNMAIIKEFMRFKVHMEG",
                },
            }
            for i in range(15)
        ]

    def test_synthetic_negatives_and_hard_rejection_scorecard(self) -> None:
        negs = generate_synthetic_negatives(
            "GGGGSGGGGSGGGGS",
            target_length=15,
            requested_rigidity="flexible",
            left_flank="THGMDELYK",
            right_flank="MVSKGEEDN",
        )
        self.assertEqual(set(negs.keys()), set(SYNTHETIC_NEGATIVE_TYPES))

        pos_card = evaluate_hard_rejection_scorecard(
            "GGGGSGGGGSGGGGS",
            required_path_length_A=28.4,
            target_length=15,
            requested_rigidity="flexible",
            requested_cleavable=False,
        )
        self.assertTrue(pos_card.hard_rejection_pass)

        # Direct 0-aa fusion must fail physical span closure when required_path_length_A = 28.4 A
        zero_card = evaluate_hard_rejection_scorecard(
            "",
            required_path_length_A=28.4,
            target_length=15,
        )
        self.assertFalse(zero_card.hard_rejection_pass)
        self.assertFalse(zero_card.physical_span_closure_pass)

        # All-hydrophobic negative must fail hydrophobic stretch check
        hydro_card = evaluate_hard_rejection_scorecard(
            negs["all_hydrophobic_core"],
            required_path_length_A=28.4,
            target_length=15,
        )
        self.assertFalse(hydro_card.hard_rejection_pass)
        self.assertFalse(hydro_card.hydrophobic_stretch_pass)

        # Buried protease inside rigid helix must fail composite rigid+cleavable check
        buried_card = evaluate_hard_rejection_scorecard(
            negs["protease_site_buried_in_rigid_helix"],
            required_path_length_A=28.4,
            target_length=15,
            requested_rigidity="rigid",
            requested_cleavable=True,
        )
        self.assertFalse(buried_card.composite_rigid_cleavable_pass)

    def test_pairwise_ranker_and_six_baselines_benchmark(self) -> None:
        res = evaluate_six_baselines_benchmark(self.records, n_folds=5, seed=23)
        self.assertEqual(set(res["six_baselines_comparison"].keys()), set(BASELINE_NAMES))
        lf_stats = res["six_baselines_comparison"]["linkerforge_ranker"]
        direct_stats = res["six_baselines_comparison"]["direct_0aa_fusion"]
        self.assertGreater(lf_stats["mean_compatibility_score"], direct_stats["mean_compatibility_score"])
        self.assertGreater(lf_stats["mean_hard_rejection_pass"], direct_stats["mean_hard_rejection_pass"])
        self.assertGreaterEqual(
            res["natural_vs_6_synthetic_corruptions"]["top1_natural_rank_accuracy"],
            0.85,
        )

    def test_sfgfp_mcherry_fusion_demo_and_3d_pdb_pymol_generation(self) -> None:
        ranker = fit_pairwise_compatibility_ranker(self.records, epochs=40)
        auth_index = AuthoritySequenceIndex(self.records)
        demo = run_sfgfp_mcherry_fusion_demo(
            ranker,
            auth_index,
            requested_rigidity="flexible",
            requested_cleavable=False,
            target_length=15,
            required_path_length_A=28.4,
        )
        self.assertEqual(demo["domain_a"]["pdb_id"], "2B3P")
        self.assertEqual(demo["domain_b"]["pdb_id"], "2H5Q")
        self.assertTrue(demo["top_selected_linker"]["hard_rejection_pass"])
        self.assertIn("ATOM      1  N   MET A   1", demo["fusion_3d_structure"]["pdb_text"])
        self.assertIn("select designed_linker", demo["fusion_3d_structure"]["pymol_pml_script"])

        bundle = build_fusion_3d_pdb_and_pymol_script(
            domain_a_name="sfGFP",
            domain_a_sequence="MSKGEELFTG",
            linker_sequence="GGGGS",
            domain_b_name="mCherry",
            domain_b_sequence="MVSKGEEDNM",
        )
        self.assertEqual(bundle["total_residues"], 25)
        self.assertEqual(bundle["total_atoms"], 100)


if __name__ == "__main__":
    unittest.main()
