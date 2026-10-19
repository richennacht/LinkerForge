"""Evaluate Part 4 (Component 5): Pairwise Compatibility Ranker, Hard Structural Rejection
Scorecard, 6 Required Baselines, and the sfGFP (2B3P) — linker — mCherry (2H5Q) Fusion Demo.

Runs:
  1. 5-fold cluster-disjoint (`homology_cluster_id`) cross-validation of
     `PairwiseCompatibilityRanker` and the 6 required baselines (`direct_0aa_fusion`,
     `fixed_ggggs3`, `random_length_matched_gs`, `geometry_only`,
     `nearest_neighbour_retrieval`, `linkerforge_ranker`) on the 148 natural linker
     authorities (`91` CATH + `57` DisProt) and on 600 high-confidence proteome-mined
     linkers (`mined_high_confidence_dfl_call`) from `mined_and_natural_linkers_v1`.
  2. End-to-end `sfGFP (PDB 2B3P) — linker — mCherry (PDB 2H5Q)` fusion demonstration
     across `flexible`, `rigid`, and `flexible+cleavable (TEV)` user modes, exporting
     the ranked candidate table, deterministic hard-rejection scorecard, and 3D `.pdb` +
     PyMOL `.pml` scripts under `reports/structures/`.
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
from linkerforge.pipeline import (  # noqa: E402
    AuthoritySequenceIndex,
    evaluate_six_baselines_benchmark,
    fit_pairwise_compatibility_ranker,
    run_sfgfp_mcherry_fusion_demo,
)


def run_evaluation(dataset_dir: Path, out_path: Path) -> dict:
    manifest = load_dataset_manifest(dataset_dir)
    manifest_sha = hashlib.sha256((dataset_dir / "manifest.json").read_bytes()).hexdigest()

    natural_records = load_linker_records(
        dataset_dir,
        tiers=["natural_cath_pdb_verified", "natural_disprot_curated_dfl"],
    )
    mined_hi = load_linker_records(
        dataset_dir,
        tiers=["mined_high_confidence_dfl_call"],
    )
    templates = load_linker_records(
        dataset_dir,
        tiers=["literature_engineered_template"],
    )

    # 1. Held-out cluster-disjoint evaluation on 148 natural authorities
    bench_natural = evaluate_six_baselines_benchmark(
        natural_records,
        n_folds=5,
        seed=23,
    )

    # 2. Held-out cluster-disjoint evaluation on 600 proteome-mined high-confidence linkers
    bench_mined = evaluate_six_baselines_benchmark(
        mined_hi[:600],
        n_folds=5,
        seed=23,
    )

    # Fit full production ranker and authority index on natural + high-confidence mined authorities
    full_train = list(natural_records) + list(mined_hi[:1200]) + list(templates)
    ranker = fit_pairwise_compatibility_ranker(full_train, seed=23)
    full_auth_index = AuthoritySequenceIndex(
        list(natural_records) + list(mined_hi) + list(templates)
    )

    # 3. sfGFP (2B3P) — linker — mCherry (2H5Q) end-to-end demonstrations
    demo_flexible = run_sfgfp_mcherry_fusion_demo(
        ranker,
        full_auth_index,
        requested_rigidity="flexible",
        requested_cleavable=False,
        target_length=15,
        required_path_length_A=28.4,
    )
    demo_rigid = run_sfgfp_mcherry_fusion_demo(
        ranker,
        full_auth_index,
        requested_rigidity="rigid",
        requested_cleavable=False,
        target_length=18,
        required_path_length_A=28.4,
    )
    demo_cleavable = run_sfgfp_mcherry_fusion_demo(
        ranker,
        full_auth_index,
        requested_rigidity="flexible",
        requested_cleavable=True,
        protease_family="TEV",
        target_length=19,
        required_path_length_A=28.4,
    )

    # Write 3D PDB and PyMOL .pml files for the top sfGFP-mCherry constructs
    struct_dir = out_path.parent / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)
    for mode_key, demo_obj in [
        ("sfgfp_mcherry_flexible", demo_flexible),
        ("sfgfp_mcherry_rigid", demo_rigid),
        ("sfgfp_mcherry_tev_cleavable", demo_cleavable),
    ]:
        pdb_file = struct_dir / f"{mode_key}.pdb"
        pml_file = struct_dir / f"{mode_key}.pml"
        pdb_file.write_text(demo_obj["fusion_3d_structure"]["pdb_text"], encoding="utf-8")
        pml_file.write_text(demo_obj["fusion_3d_structure"]["pymol_pml_script"], encoding="utf-8")
        # Keep JSON concise by storing line counts and relative artifact paths alongside preview
        demo_obj["fusion_3d_structure"]["pdb_artifact_relpath"] = f"reports/structures/{mode_key}.pdb"
        demo_obj["fusion_3d_structure"]["pml_artifact_relpath"] = f"reports/structures/{mode_key}.pml"
        del demo_obj["fusion_3d_structure"]["pdb_text"]

    payload = {
        "schema_version": "1.0.0",
        "dataset_id": manifest["dataset_id"],
        "dataset_manifest_sha256": manifest_sha,
        "authority_index_size": len(natural_records) + len(mined_hi) + len(templates),
        "ranker_parameters": {
            "feature_names": list(ranker.feature_names),
            "weights": list(ranker.weights),
            "outputs_fusion_success_probability": False,
        },
        "benchmark_natural_authorities_148": bench_natural,
        "benchmark_proteome_mined_high_confidence_600": bench_mined,
        "sfgfp_2b3p_mcherry_2h5q_demos": {
            "flexible_15aa": demo_flexible,
            "rigid_18aa": demo_rigid,
            "flexible_tev_cleavable_19aa": demo_cleavable,
        },
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR)
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "reports" / "part4_end_to_end_ranker_and_sfgfp_mcherry_2026-10-19.json",
    )
    args = parser.parse_args()
    res = run_evaluation(args.dataset_dir, args.out)
    b_nat = res["benchmark_natural_authorities_148"]
    print("Part 4 Natural Authorities (148) vs 6 Synthetic Corruptions:")
    print(
        f"  Top-1 Natural Rank = {b_nat['natural_vs_6_synthetic_corruptions']['top1_natural_rank_accuracy']:.4f} "
        f"| MRR = {b_nat['natural_vs_6_synthetic_corruptions']['mrr_natural_rank']:.4f}"
    )
    print("Part 4 Six Baselines Comparison (Natural Authorities 148):")
    for m, stats in b_nat["six_baselines_comparison"].items():
        print(
            f"  {m:30s} | HardPass={stats['mean_hard_rejection_pass']:.4f} "
            f"| Score={stats['mean_compatibility_score']:.4f} "
            f"| WinVsNegs={stats['mean_pairwise_win_rate_vs_synthetic_negatives']:.4f}"
        )


if __name__ == "__main__":
    main()
