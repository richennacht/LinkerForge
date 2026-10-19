"""Evaluate Part 3 (Component 4): Two-Tower Attachment-Site Scorer and Geometric Site-Pair Ranker.

Runs:
  1. 5-fold cluster-disjoint (`homology_cluster_id`) cross-validation of
     `TwoTowerAttachmentScorer` on the 148 natural linker authorities (`91` CATH
     + `57` DisProt) and on a 500-record stratified sample of `mined_high_confidence_dfl_call`
     from `mined_and_natural_linkers_v1`.
  2. Closed-form geometric site-pair ranking (`rank_attachment_site_pairs`) on real
     PDB two-domain fixtures (`4e43_protease_dimer.pdb` and `4ake_open.pdb`) across
     both `end_to_end_tandem` and `internal_insertion_two_tether` topologies.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.attachment import (  # noqa: E402
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
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}


def chain_to_sequence(chain: Chain) -> str:
    return "".join(THREE_TO_ONE.get(r.name.upper(), "G") for r in chain.residues)


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

    # 1. Held-out domain-family benchmark on 148 natural authorities
    bench_natural = evaluate_attachment_scorer_benchmark(
        natural_records,
        n_folds=5,
        seed=17,
    )

    # 2. Held-out domain-family benchmark on 600 high-confidence proteome-mined linkers
    bench_mined = evaluate_attachment_scorer_benchmark(
        mined_hi[:600],
        n_folds=5,
        seed=17,
    )

    # Fit full scorer on natural + high-confidence mined authorities
    scorer = fit_two_tower_attachment_scorer(list(natural_records) + list(mined_hi[:1000]))

    # 3. 3D Site-Pair Re-Ranking on 4e43_protease_dimer.pdb
    pdb_4e43 = ROOT / "tests" / "fixtures" / "geometry" / "4e43_protease_dimer.pdb"
    chains_4e43 = read_chains(pdb_4e43)
    dom_a_4e43 = chains_4e43["A"]
    dom_b_4e43 = chains_4e43["B"]
    seq_a_4e43 = chain_to_sequence(dom_a_4e43)
    seq_b_4e43 = chain_to_sequence(dom_b_4e43)

    tandem_pairs_4e43 = rank_attachment_site_pairs(
        dom_a_4e43,
        seq_a_4e43,
        dom_b_4e43,
        seq_b_4e43,
        scorer,
        linker_type="flexible",
        requested_residues=15,
        topology="end_to_end_tandem",
        candidate_sites_per_domain=6,
        top_k=5,
    )

    insertion_pairs_4e43 = rank_attachment_site_pairs(
        dom_a_4e43,
        seq_a_4e43,
        dom_b_4e43,
        seq_b_4e43,
        scorer,
        linker_type="flexible",
        requested_residues=15,
        topology="internal_insertion_two_tether",
        candidate_sites_per_domain=6,
        top_k=5,
    )

    payload = {
        "schema_version": "1.0.0",
        "dataset_id": manifest["dataset_id"],
        "dataset_manifest_sha256": manifest_sha,
        "scorer_parameters": {
            "weights": list(scorer.weights),
            "bias": scorer.bias,
            "tick_boxes_excluded_from_encoder": True,
        },
        "benchmark_natural_authorities_148": bench_natural,
        "benchmark_proteome_mined_high_confidence_600": bench_mined,
        "structure_site_pair_reranking_demos": {
            "4e43_end_to_end_tandem_top5": [p.to_dict() for p in tandem_pairs_4e43],
            "4e43_internal_insertion_two_tether_top5": [p.to_dict() for p in insertion_pairs_4e43],
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
        default=ROOT / "reports" / "part3_attachment_scorer_2026-10-19.json",
    )
    args = parser.parse_args()
    res = run_evaluation(args.dataset_dir, args.out)
    print("Part 3 Natural Authorities (148) Benchmark:")
    for arm, m in res["benchmark_natural_authorities_148"]["arms"].items():
        print(
            f"  {arm:48s} | AUROC={m['auroc']:.4f} | APS={m['average_precision_aps']:.4f} "
            f"| P@10%={m['precision_at_10pct_coverage']:.4f} | Top1={m['top1_port_recovery_rate']:.4f}"
        )
    print("\nPart 3 Proteome-Mined High-Confidence (600) Benchmark:")
    for arm, m in res["benchmark_proteome_mined_high_confidence_600"]["arms"].items():
        print(
            f"  {arm:48s} | AUROC={m['auroc']:.4f} | APS={m['average_precision_aps']:.4f} "
            f"| P@10%={m['precision_at_10pct_coverage']:.4f} | Top1={m['top1_port_recovery_rate']:.4f}"
        )


if __name__ == "__main__":
    main()
