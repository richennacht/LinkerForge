"""Unit tests for Part 2 Geometry-Plus-Template Design and X2 Retrieval Baseline."""

from __future__ import annotations

import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.dataset_v1 import DEFAULT_DATASET_DIR
from linkerforge.design import (
    DesignRequest,
    build_rigid_cleavable_composite,
    design_linker_candidates,
)
from linkerforge.geometry.structure import read_chain
from linkerforge.properties import load_library
from linkerforge.retrieval import (
    evaluate_x2_retrieval_benchmark,
    flanking_pair_similarity,
    sequence_nw_identity,
)


class TestPart2DesignAndRetrieval(unittest.TestCase):
    def test_rigid_cleavable_composite_enforces_two_part_construction(self) -> None:
        lib = load_library()
        composites = build_rigid_cleavable_composite(25.0, lib, protease="MMP")
        self.assertGreater(len(composites), 0)
        top = composites[0]
        self.assertEqual(
            top["construction_mode"],
            "two_part_rigid_spacer_plus_flexible_protease_cassette",
        )
        self.assertIn("GGGGSPLGLWAGGGGS", top["sequence"])
        self.assertTrue(
            top["sequence"].startswith("AEAAAK") or top["sequence"].startswith("PAPAP") or top["sequence"].startswith("APAP")
        )

    def test_design_engine_supports_tandem_and_internal_insertion_topologies(self) -> None:
        if not (DEFAULT_DATASET_DIR / "manifest.json").exists():
            self.skipTest("mined_and_natural_linkers_v1 dataset not present")
        # 1. Tandem flexible non-cleavable with 3D structure fixture
        pdb_path = ROOT / "tests" / "fixtures" / "geometry" / "4e43_protease_dimer.pdb"
        chain_a = read_chain(pdb_path, "A")
        chain_b = read_chain(pdb_path, "B")
        report_tandem = design_linker_candidates(
            DesignRequest(
                distance_angstrom=18.0,
                rigidity="flexible",
                cleavability="non_cleavable",
                topology="end_to_end_tandem",
            ),
            domain_a=chain_a,
            anchor_a=99,
            domain_b=chain_b,
            anchor_b=1,
            top_k=3,
        )
        self.assertEqual(report_tandem.q1_feasibility["topology"], "end_to_end_tandem")
        self.assertIsNotNone(report_tandem.q2_attachment_and_context["structure_junction_report"])
        self.assertEqual(len(report_tandem.q3_candidates), 3)
        for cand in report_tandem.q3_candidates:
            self.assertEqual(cand["geometric_verdict"], "bridgeable")
            self.assertEqual(cand["rigidity_predicted_label"], "flexible")
            self.assertEqual(cand["cleavability_label"], "no_known_motif")

        # 2. Internal insertion two-tether topology
        report_insert = design_linker_candidates(
            DesignRequest(
                distance_angstrom=14.0,
                return_distance_angstrom=16.0,
                rigidity="flexible",
                cleavability="non_cleavable",
                topology="internal_insertion_two_tether",
            ),
            top_k=2,
        )
        self.assertEqual(
            report_insert.q1_feasibility["topology"], "internal_insertion_two_tether"
        )
        self.assertIsNotNone(report_insert.q1_feasibility["return_tether_span_window"])
        self.assertIn(
            "requires_second_return_tether_for_internal_insertion",
            report_insert.q3_candidates[0]["flags"],
        )

        # 3. Rigid + cleavable returns two-part composite constructs
        report_rc = design_linker_candidates(
            DesignRequest(
                distance_angstrom=26.0,
                rigidity="rigid",
                cleavability="cleavable",
                topology="end_to_end_tandem",
            ),
            top_k=3,
        )
        self.assertGreater(len(report_rc.q3_candidates), 0)
        self.assertEqual(
            report_rc.q3_candidates[0]["construction_mode"],
            "two_part_rigid_spacer_plus_flexible_protease_cassette",
        )
        self.assertEqual(report_rc.q3_candidates[0]["cleavability_label"], "cleavable")

    def test_flanking_pair_similarity_and_cluster_disjoint_x2_benchmark(self) -> None:
        flank_1 = {
            "left_flank_sequence": "AEDWCAVCQNGGELLCCEKCPKVFHLSCHVPTLTNFPSGEWICTFCRDLSK",
            "right_flank_sequence": "KKTEGLVKLTPIDKRKCERLLLFLYCHEMSLAFQDPVPLTVPDYYKIIKNP",
            "left_domain_id": "3.30.40.10",
            "right_domain_id": "1.20.920.10",
        }
        flank_2 = {
            "left_flank_sequence": "PKVFHLSCHVPTLTNFPSGEWICTFCRDLSK",
            "right_flank_sequence": "KKTEGLVKLTPIDKRKCERLLLFLYCHEMSLA",
            "left_domain_id": "3.30.40.10",
            "right_domain_id": "1.20.920.10",
        }
        flank_diff = {
            "left_flank_sequence": "GGGGGGGGGGGGGGGGGGGG",
            "right_flank_sequence": "PPPPPPPPPPPPPPPPPPPP",
            "left_domain_id": "2.60.40.10",
            "right_domain_id": "2.40.10.10",
        }
        sim_close = flanking_pair_similarity(flank_1, flank_2)
        sim_far = flanking_pair_similarity(flank_1, flank_diff)
        self.assertGreater(sim_close, sim_far)
        self.assertGreater(sequence_nw_identity("GGGGSGGGGS", "GGGGSGGGGS"), 0.99)

        # Synthetic 6-cluster benchmark to verify zero cluster leakage and all arms run
        synth_natural = []
        for i in range(6):
            synth_natural.append(
                {
                    "record_id": f"nat:{i}",
                    "homology_cluster_id": f"cluster_{i}",
                    "linker_sequence": "GGGGSGGGGS" if i % 2 == 0 else "AEAAAKEAAAKA",
                    "linker_length": 10 if i % 2 == 0 else 12,
                    "flanking_context": flank_1 if i % 2 == 0 else flank_diff,
                    "geometry": {"inter_port_ca_distance_A": 18.0},
                    "biophysical_properties": {
                        "rigidity_predicted_label": "flexible" if i % 2 == 0 else "rigid"
                    },
                }
            )
        synth_templates = [
            {
                "record_id": "tmpl:flex2",
                "homology_cluster_id": "tmpl:1",
                "linker_sequence": "GGGGSGGGGS",
                "linker_length": 10,
                "biophysical_properties": {"rigidity_predicted_label": "flexible"},
            },
            {
                "record_id": "tmpl:rigid2",
                "homology_cluster_id": "tmpl:2",
                "linker_sequence": "AEAAAKEAAAKA",
                "linker_length": 12,
                "biophysical_properties": {"rigidity_predicted_label": "rigid"},
            },
        ]
        res = evaluate_x2_retrieval_benchmark(
            synth_natural,
            [],
            synth_templates,
            has_3d_geometry=True,
            n_folds=3,
            seed=7,
        )
        self.assertEqual(res["total_evaluated_natural_linkers"], 6)
        self.assertIn("flanking_context_template_retrieval", res["arms"])
        self.assertIn("geometry_plus_flanking_hybrid_retrieval", res["arms"])


if __name__ == "__main__":
    unittest.main()
