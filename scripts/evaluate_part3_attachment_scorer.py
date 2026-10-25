"""Evaluate Part 3: position-controlled attachment-site scorer and geometric site-pair ranker.

Outputs `reports/part3_attachment_scorer_2026-10-19.json` containing:
  1. `position_matched_benchmark_*`: cluster-disjoint 5-fold evaluation where positives (residues
     adjacent to a natural linker) and negatives (same offsets from the opposite window edge)
     have identical position statistics.  Chance AUROC = 0.5.
  2. `original_task_leakage_diagnostic`: the earlier benchmark design, solved by position alone.
  3. `structure_site_pair_reranking_demo`: closed-form 3D site-pair ranking on a real PDB
     fixture.  This is a worked example, not a validated benchmark.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.attachment import (  # noqa: E402
    diagnose_original_task_leakage,
    evaluate_attachment_scorer_benchmark,
    fit_two_tower_attachment_scorer,
    rank_attachment_site_pairs,
)
from linkerforge.dataset_v1 import (  # noqa: E402
    DEFAULT_DATASET_DIR,
    load_dataset_manifest,
    load_linker_records,
)
from linkerforge.geometry.structure import Chain, read_chains  # noqa: E402

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H",
    "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q",
    "ARG": "R", "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}
MINED_SUBSAMPLE = 2000
SUBSAMPLE_SEED = 20261019


def chain_to_sequence(chain: Chain) -> str:
    return "".join(THREE_TO_ONE.get(r.name.upper(), "X") for r in chain.residues)


def run_evaluation(dataset_dir: Path, out_path: Path) -> dict:
    manifest = load_dataset_manifest(dataset_dir)
    manifest_sha = hashlib.sha256((dataset_dir / "manifest.json").read_bytes()).hexdigest()

    natural = load_linker_records(
        dataset_dir, tiers=["natural_cath_pdb_verified", "natural_disprot_curated_dfl"]
    )
    mined_hi_all = list(load_linker_records(dataset_dir, tiers=["mined_high_confidence_dfl_call"]))
    rng = random.Random(SUBSAMPLE_SEED)
    mined_hi = rng.sample(mined_hi_all, min(MINED_SUBSAMPLE, len(mined_hi_all)))

    bench_natural = evaluate_attachment_scorer_benchmark(natural, n_folds=5, seed=17)
    bench_mined = evaluate_attachment_scorer_benchmark(mined_hi, n_folds=5, seed=17)
    leakage = diagnose_original_task_leakage(list(natural) + mined_hi)

    scorer = fit_two_tower_attachment_scorer(list(natural) + mined_hi)

    pdb_4e43 = ROOT / "tests" / "fixtures" / "geometry" / "4e43_protease_dimer.pdb"
    chains = read_chains(pdb_4e43)
    dom_a, dom_b = chains["A"], chains["B"]
    seq_a, seq_b = chain_to_sequence(dom_a), chain_to_sequence(dom_b)
    demo = {}
    for topology in ("end_to_end_tandem", "internal_insertion_two_tether"):
        ranked = rank_attachment_site_pairs(
            dom_a, seq_a, dom_b, seq_b, scorer,
            linker_type="flexible", requested_residues=15, topology=topology,
            candidate_sites_per_domain=6, top_k=5,
        )
        demo[f"4e43_{topology}_top5"] = [p.to_dict() for p in ranked]

    payload = {
        "schema_version": "2.0.0",
        "dataset_id": manifest["dataset_id"],
        "dataset_manifest_sha256": manifest_sha,
        "scorer_parameters": {
            "weights": list(scorer.weights),
            "bias": scorer.bias,
            "positional_feature_weights_fixed_to_zero": True,
            "tick_boxes_excluded_from_encoder": True,
        },
        "mined_subsample": {
            "tier": "mined_high_confidence_dfl_call",
            "sampled": len(mined_hi),
            "available": len(mined_hi_all),
            "seed": SUBSAMPLE_SEED,
        },
        "position_matched_benchmark_natural_authorities": bench_natural,
        "position_matched_benchmark_mined_high_confidence_subsample": bench_mined,
        "original_task_leakage_diagnostic": leakage,
        "structure_site_pair_reranking_demo": {
            "status": "worked example on one PDB fixture; no held-out validation",
            "fixture": "tests/fixtures/geometry/4e43_protease_dimer.pdb",
            "results": demo,
        },
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR)
    parser.add_argument(
        "--out", type=Path, default=ROOT / "reports" / "part3_attachment_scorer_2026-10-19.json"
    )
    args = parser.parse_args()
    res = run_evaluation(args.dataset_dir, args.out)
    for key in (
        "position_matched_benchmark_natural_authorities",
        "position_matched_benchmark_mined_high_confidence_subsample",
    ):
        b = res[key]
        print(f"{key}: windows={b['n_windows']} clusters={b['n_homology_clusters']}")
        for arm, m in b["arms"].items():
            print(f"  {arm:40s} AUROC={m['auroc_pooled']:.4f} CI={m['auroc_cluster_bootstrap_95ci']}")
    d = res["original_task_leakage_diagnostic"]
    print(
        f"original-design position-only AUROC: term_prox={d['auroc_position_only_term_prox']} "
        f"rel_interior={d['auroc_position_only_rel_interior']}"
    )


if __name__ == "__main__":
    main()
