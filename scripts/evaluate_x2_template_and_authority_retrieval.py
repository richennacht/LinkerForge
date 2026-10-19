"""Evaluate Part 2 Geometry-Plus-Template Design and Experiment X2 Retrieval Baseline.

Runs:
  1. Experiment X2 (`CHECKLIST.md` item X2 & `CURRENT_DIRECTION_2026-10-18.md` Step 4):
     5-fold cluster-disjoint (`homology_cluster_id`-disjoint) cross-validation
     reconstructing held-out natural linkers from flanking domain context (and 3D
     attachment geometry on Cohort A), comparing:
       - `fixed_ggggs3_control`
       - `training_median_length_gs_control`
       - `flanking_context_template_retrieval`
       - `flanking_context_nn_natural_authority`
       - `flanking_context_nn_natural_plus_mined`
       - `geometry_only_control` (Cohort A)
       - `geometry_plus_flanking_hybrid_retrieval` (Cohort A)
       - `length_matched_template_control` (Oracle-length control)
       - `length_matched_frequent_natural_control` (Oracle-length control)
  2. Part 2 Geometry-Plus-Template Design Verification (`CURRENT_DIRECTION_2026-10-18.md`
     Step 3 & `SPEC.md` section 7):
     Evaluates Q1 (feasibility), Q2 (attachment/flanking context), and Q3 (candidate
     selection) across all four `(rigidity, cleavability)` regimes (including the
     mandatory two-part composite construction for `rigid + cleavable`), both
     topologies (`end_to_end_tandem` and `internal_insertion_two_tether`), and the
     real 3D PDB fixture `tests/fixtures/geometry/4e43_protease_dimer.pdb`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.dataset_v1 import (  # noqa: E402
    DEFAULT_DATASET_DIR,
    load_dataset_manifest,
    load_linker_records,
)
from linkerforge.design import DesignRequest, design_linker_candidates  # noqa: E402
from linkerforge.geometry.structure import read_chain  # noqa: E402
from linkerforge.retrieval import evaluate_x2_retrieval_benchmark  # noqa: E402


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_evaluation(dataset_dir: Path, output_path: Path) -> dict:
    manifest = load_dataset_manifest(dataset_dir)
    manifest_sha = sha256_file(dataset_dir / "manifest.json")

    cath_records = load_linker_records(
        dataset_dir,
        tiers=["natural_cath_pdb_verified"],
        min_length=5,
        max_length=60,
    )
    disprot_records = load_linker_records(
        dataset_dir,
        tiers=["natural_disprot_curated_dfl"],
        min_length=5,
        max_length=60,
    )
    natural_records = cath_records + disprot_records
    mined_records = load_linker_records(
        dataset_dir,
        tiers=["mined_high_confidence_dfl_call", "mined_provisional_dfl_candidate"],
        min_length=5,
        max_length=60,
    )
    mined_high_only = load_linker_records(
        dataset_dir,
        tiers=["mined_high_confidence_dfl_call"],
        min_length=5,
        max_length=60,
    )
    template_records = load_linker_records(
        dataset_dir,
        tiers=["literature_engineered_template"],
        min_length=5,
        max_length=60,
    )

    # Cohort A: 91 coordinate-complete CATH v4.4.0 natural two-domain linkers (57 domain-pair clusters)
    cohort_a = evaluate_x2_retrieval_benchmark(
        cath_records,
        mined_records,
        template_records,
        has_3d_geometry=True,
        n_folds=5,
        seed=17,
    )

    # Cohort B: All 148 natural authorities (91 CATH + 57 DisProt curated IDPO:0000033 across 113 clusters)
    cohort_b = evaluate_x2_retrieval_benchmark(
        natural_records,
        mined_records,
        template_records,
        has_3d_geometry=False,
        n_folds=5,
        seed=17,
    )

    # Design Engine Verification Suite across all 4 (rigidity, cleavability) regimes & 2 topologies
    pdb_fixture = ROOT / "tests" / "fixtures" / "geometry" / "4e43_protease_dimer.pdb"
    chain_a = read_chain(pdb_fixture, "A")
    chain_b = read_chain(pdb_fixture, "B")

    sample_cath = cath_records[0]
    sample_flank = sample_cath["flanking_context"]

    design_cases = {
        "flexible_non_cleavable_tandem_3d": design_linker_candidates(
            DesignRequest(
                distance_angstrom=18.5,
                rigidity="flexible",
                cleavability="non_cleavable",
                topology="end_to_end_tandem",
                left_flank_sequence=sample_flank["left_flank_sequence"],
                right_flank_sequence=sample_flank["right_flank_sequence"],
                left_domain_id=sample_flank["left_domain_id"],
                right_domain_id=sample_flank["right_domain_id"],
            ),
            dataset_dir=dataset_dir,
            domain_a=chain_a,
            anchor_a=99,
            domain_b=chain_b,
            anchor_b=1,
            top_k=5,
        ).to_dict(),
        "rigid_non_cleavable_tandem": design_linker_candidates(
            DesignRequest(
                distance_angstrom=22.0,
                rigidity="rigid",
                cleavability="non_cleavable",
                topology="end_to_end_tandem",
            ),
            dataset_dir=dataset_dir,
            top_k=5,
        ).to_dict(),
        "flexible_cleavable_tandem": design_linker_candidates(
            DesignRequest(
                distance_angstrom=20.0,
                rigidity="flexible",
                cleavability="cleavable",
                topology="end_to_end_tandem",
            ),
            dataset_dir=dataset_dir,
            top_k=5,
        ).to_dict(),
        "rigid_cleavable_composite_tandem": design_linker_candidates(
            DesignRequest(
                distance_angstrom=28.0,
                rigidity="rigid",
                cleavability="cleavable",
                topology="end_to_end_tandem",
            ),
            dataset_dir=dataset_dir,
            top_k=5,
        ).to_dict(),
        "flexible_non_cleavable_internal_insertion": design_linker_candidates(
            DesignRequest(
                distance_angstrom=15.0,
                return_distance_angstrom=16.5,
                rigidity="flexible",
                cleavability="non_cleavable",
                topology="internal_insertion_two_tether",
            ),
            dataset_dir=dataset_dir,
            top_k=3,
        ).to_dict(),
    }

    payload = {
        "schema_version": "1.0.0",
        "experiment_id": "X2_retrieval_baseline_and_part2_geometry_template_design",
        "dataset_provenance": {
            "dataset_id": manifest["dataset_id"],
            "manifest_sha256": manifest_sha,
            "natural_cath_count": len(cath_records),
            "natural_disprot_count": len(disprot_records),
            "mined_total_count": len(mined_records),
            "mined_high_confidence_count": len(mined_high_only),
            "template_count": len(template_records),
        },
        "x2_retrieval_benchmark": {
            "cohort_a_cath_3d_natural_91": cohort_a,
            "cohort_b_all_natural_authorities_148": cohort_b,
        },
        "part2_design_engine_verification": {
            "fixture_pdb_sha256": sha256_file(pdb_fixture),
            "cases": design_cases,
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        default=DEFAULT_DATASET_DIR,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "reports" / "x2_retrieval_and_part2_design_2026-10-19.json",
    )
    args = parser.parse_args()
    result = run_evaluation(args.dataset_dir, args.output)
    print("Cohort A (CATH 3D 91) arms:")
    for arm, stats in result["x2_retrieval_benchmark"]["cohort_a_cath_3d_natural_91"]["arms"].items():
        print(
            f"  {arm:42s} | MAE={stats['length_mae_residues']:.2f} aa | "
            f"<=5aa={stats['length_within_5aa_count']}/{stats['evaluated_records']} | "
            f"NW_id={stats['sequence_nw_identity_mean']:.4f} | "
            f"CompCos={stats['composition_cosine_similarity_mean']:.4f} | "
            f"3.8A_feas={stats['physical_3p8A_feasible_count']}/{stats['evaluated_records']}"
        )
    print("\nCohort B (All Natural 148) arms:")
    for arm, stats in result["x2_retrieval_benchmark"]["cohort_b_all_natural_authorities_148"]["arms"].items():
        print(
            f"  {arm:42s} | MAE={stats['length_mae_residues']:.2f} aa | "
            f"<=5aa={stats['length_within_5aa_count']}/{stats['evaluated_records']} | "
            f"NW_id={stats['sequence_nw_identity_mean']:.4f} | "
            f"CompCos={stats['composition_cosine_similarity_mean']:.4f}"
        )


if __name__ == "__main__":
    main()
