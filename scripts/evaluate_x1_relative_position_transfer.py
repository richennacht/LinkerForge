"""Experiment X1: Do relative-position 3D structural features transfer across held-out domain pairs?

Implements `CHECKLIST.md` item X1 using the 89 coordinate-complete CATH v4.4.0
natural two-domain proteins (`cath_full_pilot.jsonl` + cached RCSB PDB coordinate
files in `data/raw/rcsb/*.pdb`).

Splits examples strictly by held-out domain-pair architecture (`grouped_split`
over `homology_cluster_id` and `domain_a`/`domain_b` CATH families so zero
domain-pair architectures cross splits). Compares on identical held-out domain
pairs across 5 seeds (`7, 17, 31, 43, 97`):
  1. `sequence_only_ranker`: 6 candidate linker sequence features only.
  2. `sequence_plus_local_port_ranker`: 12 features (sequence + within-domain
     local port contact density and accessibility proxy).
  3. `sequence_plus_relative_position_3d_ranker`: 18 features adding Li et al.-style
     inter-domain 3D relative-position C-alpha geometry (inter-port CA distance,
     physical 3.8 A/res span slack, inter-centroid distance, relative port-to-partner
     centroid ratios, and terminal exit/entry vector cosine alignment).
  4. Non-learned baselines: `fixed_gs`, `length_heuristic`, `nearest_training`,
     and `geometry_span_match`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from linkerforge.baselines import score_candidates
from linkerforge.candidates import candidate_pool
from linkerforge.models.context import CHARGED, HYDROPHOBIC, LocalPort
from linkerforge.splits import audit_split, grouped_split

if __package__:
    from scripts.build_rcsb_pilot import read_ca, sha256_bytes
else:
    from build_rcsb_pilot import read_ca, sha256_bytes


SEQ_FEATURES = (
    "length_fraction",
    "glycine_fraction",
    "serine_fraction",
    "hydrophobic_fraction",
    "charged_fraction",
    "sequence_entropy",
)

PORT_FEATURES = SEQ_FEATURES + (
    "port_a_accessibility",
    "port_b_accessibility",
    "port_a_contact_density",
    "port_b_contact_density",
    "length_by_mean_accessibility",
    "hydrophobic_by_mean_accessibility",
)

REL_POS_3D_FEATURES = PORT_FEATURES + (
    "length_minus_inter_port_span",
    "abs_length_minus_inter_port_span",
    "physical_slack_3p8A_fraction",
    "length_by_inter_centroid_ratio",
    "length_by_relative_port_asymmetry",
    "gly_ser_by_terminal_alignment",
)


def _vec_sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _norm(v: tuple[float, float, float]) -> float:
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def _unit_vec(v: tuple[float, float, float]) -> tuple[float, float, float]:
    n = _norm(v)
    return (0.0, 0.0, 0.0) if n < 1e-8 else (v[0] / n, v[1] / n, v[2] / n)


def _dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _centroid(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    n = len(points)
    return (
        sum(p[0] for p in points) / n,
        sum(p[1] for p in points) / n,
        sum(p[2] for p in points) / n,
    )


def _radius_of_gyration(points: list[tuple[float, float, float]], center: tuple[float, float, float]) -> float:
    return math.sqrt(sum(math.dist(p, center) ** 2 for p in points) / max(1, len(points)))


def extract_relative_position_geometry(record: dict, rcsb_cache: Path) -> dict[str, float]:
    """Compute 3D inter-domain relative-position features from the verified PDB file."""
    pdb_id = record["source"]["pdb_id"].upper()
    chain_id = record["source"]["chain_id"]
    pdb_path = rcsb_cache / f"{pdb_id}.pdb"
    payload = pdb_path.read_bytes()
    if sha256_bytes(payload) != record["source"]["coordinate_sha256"]:
        raise ValueError(f"PDB coordinate SHA-256 mismatch for {pdb_id}.{chain_id}")
    ca = read_ca(payload, chain_id)
    first, last = record["source"]["cath_author_numbered_range"]
    a_end_seq = record["domain_a"]["end"]
    b_start_seq = record["domain_b"]["start"]
    a_port_num = first + a_end_seq - 1
    b_port_num = first + b_start_seq - 1

    dom_a_pts = [ca[n][1] for n in range(first, a_port_num + 1) if n in ca]
    dom_b_pts = [ca[n][1] for n in range(b_port_num, last + 1) if n in ca]
    if a_port_num not in ca or b_port_num not in ca or len(dom_a_pts) < 5 or len(dom_b_pts) < 5:
        raise ValueError(f"Insufficient resolved domain coordinates for {pdb_id}.{chain_id}")

    port_a_xyz = ca[a_port_num][1]
    port_b_xyz = ca[b_port_num][1]
    inter_port_dist = math.dist(port_a_xyz, port_b_xyz)

    centroid_a = _centroid(dom_a_pts)
    centroid_b = _centroid(dom_b_pts)
    rg_a = max(1.0, _radius_of_gyration(dom_a_pts, centroid_a))
    rg_b = max(1.0, _radius_of_gyration(dom_b_pts, centroid_b))
    inter_centroid_dist = max(1.0, math.dist(centroid_a, centroid_b))

    rel_port_a_to_b = math.dist(port_a_xyz, centroid_b) / inter_centroid_dist
    rel_port_b_to_a = math.dist(port_b_xyz, centroid_a) / inter_centroid_dist
    inter_centroid_ratio = min(2.0, inter_centroid_dist / (rg_a + rg_b)) / 2.0

    u_ab = _unit_vec(_vec_sub(port_b_xyz, port_a_xyz))
    prev_a = ca.get(a_port_num - 1, ca[a_port_num])[1]
    next_b = ca.get(b_port_num + 1, ca[b_port_num])[1]
    exit_a = _unit_vec(_vec_sub(port_a_xyz, prev_a))
    entry_b = _unit_vec(_vec_sub(next_b, port_b_xyz))
    terminal_alignment = ((_dot(exit_a, u_ab) + _dot(entry_b, u_ab)) / 2.0 + 1.0) / 2.0

    ports = record["connection_features"]
    return {
        "port_a_accessibility": float(ports["port_a"]["accessibility_proxy"]),
        "port_b_accessibility": float(ports["port_b"]["accessibility_proxy"]),
        "port_a_contact_density": float(ports["port_a"]["contact_density"]),
        "port_b_contact_density": float(ports["port_b"]["contact_density"]),
        "inter_port_ca_distance_A": inter_port_dist,
        "inter_centroid_distance_A": inter_centroid_dist,
        "inter_centroid_ratio": inter_centroid_ratio,
        "rel_port_asymmetry": min(1.0, abs(rel_port_a_to_b - rel_port_b_to_a)),
        "terminal_alignment": max(0.0, min(1.0, terminal_alignment)),
    }


def compute_feature_vector(geom: dict[str, float], linker_sequence: str, feature_names: tuple[str, ...]) -> tuple[float, ...]:
    n = len(linker_sequence)
    counts = {aa: linker_sequence.count(aa) for aa in set(linker_sequence)}
    entropy = -sum((c / n) * math.log2(c / n) for c in counts.values()) / math.log2(20)
    length_fraction = n / 30.0
    g_frac = linker_sequence.count("G") / n
    s_frac = linker_sequence.count("S") / n
    h_frac = sum(aa in HYDROPHOBIC for aa in linker_sequence) / n
    c_frac = sum(aa in CHARGED for aa in linker_sequence) / n
    mean_acc = (geom["port_a_accessibility"] + geom["port_b_accessibility"]) / 2.0

    # Physical 3.8 A/residue bound (SPEC section 2.4) and empirical 2.5 A/residue flexible expectation
    inter_port_A = geom["inter_port_ca_distance_A"]
    expected_span_residues = min(30.0, inter_port_A / 2.5) / 30.0
    span_diff = length_fraction - expected_span_residues
    slack_3p8 = max(-1.0, min(1.0, (3.8 * n - inter_port_A) / (30.0 * 3.8)))

    all_feats = {
        "length_fraction": length_fraction,
        "glycine_fraction": g_frac,
        "serine_fraction": s_frac,
        "hydrophobic_fraction": h_frac,
        "charged_fraction": c_frac,
        "sequence_entropy": entropy,
        "port_a_accessibility": geom["port_a_accessibility"],
        "port_b_accessibility": geom["port_b_accessibility"],
        "port_a_contact_density": geom["port_a_contact_density"],
        "port_b_contact_density": geom["port_b_contact_density"],
        "length_by_mean_accessibility": length_fraction * mean_acc,
        "hydrophobic_by_mean_accessibility": h_frac * mean_acc,
        "length_minus_inter_port_span": span_diff,
        "abs_length_minus_inter_port_span": abs(span_diff),
        "physical_slack_3p8A_fraction": slack_3p8,
        "length_by_inter_centroid_ratio": length_fraction * geom["inter_centroid_ratio"],
        "length_by_relative_port_asymmetry": length_fraction * geom["rel_port_asymmetry"],
        "gly_ser_by_terminal_alignment": (g_frac + s_frac) * geom["terminal_alignment"],
    }
    return tuple(all_feats[name] for name in feature_names)


def fit_pairwise_linear(
    train_records: list[dict],
    geom_by_id: dict[str, dict[str, float]],
    feature_names: tuple[str, ...],
    *,
    epochs: int = 50,
    lr: float = 0.05,
    l2: float = 0.001,
) -> tuple[float, ...]:
    diffs: list[tuple[float, ...]] = []
    for rec in train_records:
        native = rec["linker"]["sequence"]
        alts = tuple(sorted({
            other["linker"]["sequence"]
            for other in train_records
            if other["example_id"] != rec["example_id"] and other["linker"]["sequence"] != native
        }))[:5]
        if not alts:
            continue
        obs_vec = compute_feature_vector(geom_by_id[rec["example_id"]], native, feature_names)
        for alt in alts:
            neg_vec = compute_feature_vector(geom_by_id[rec["example_id"]], alt, feature_names)
            diff = tuple(o - n for o, n in zip(obs_vec, neg_vec))
            if any(diff):
                diffs.append(diff)
    weights = [0.0] * len(feature_names)
    for _ in range(epochs):
        for diff in diffs:
            margin = sum(w * d for w, d in zip(weights, diff))
            weights = [w * (1.0 - lr * l2) for w in weights]
            if margin < 1.0:
                weights = [w + lr * d for w, d in zip(weights, diff)]
    return tuple(weights)


def rank_candidates_linear(
    rec: dict,
    geom: dict[str, float],
    benchmark_pool: tuple[dict, ...],
    weights: tuple[float, ...],
    feature_names: tuple[str, ...],
    ref_id: str,
) -> int:
    scored = sorted(
        (
            {
                "candidate_id": c["candidate_id"],
                "score": sum(
                    w * x
                    for w, x in zip(
                        weights,
                        compute_feature_vector(geom, c["linker_sequence"], feature_names),
                    )
                ),
            }
            for c in benchmark_pool
        ),
        key=lambda item: (-item["score"], item["candidate_id"]),
    )
    return next(i for i, item in enumerate(scored, 1) if item["candidate_id"] == ref_id)


def rank_geometry_span_match(geom: dict[str, float], benchmark_pool: tuple[dict, ...], ref_id: str) -> int:
    dist_A = geom["inter_port_ca_distance_A"]
    scored = []
    for c in benchmark_pool:
        n = len(c["linker_sequence"])
        feasible = 1.0 if 3.8 * n >= dist_A else 0.0
        target_n = dist_A / 2.5
        score = feasible * 100.0 - abs(n - target_n)
        scored.append({"candidate_id": c["candidate_id"], "score": score})
    scored.sort(key=lambda item: (-item["score"], item["candidate_id"]))
    return next(i for i, item in enumerate(scored, 1) if item["candidate_id"] == ref_id)


def cluster_bootstrap_diff(
    rows: list[dict],
    key_proposal: str,
    key_baseline: str,
    *,
    draws: int = 500,
    seed: int = 17,
) -> dict[str, float]:
    by_group: dict[str, list[dict]] = {}
    for row in rows:
        by_group.setdefault(row["group"], []).append(row)
    groups = sorted(by_group.keys())
    rng = random.Random(seed)
    diffs = []
    for _ in range(draws):
        sampled = [ rng.choice(groups) for _ in range(len(groups)) ]
        sampled_rows = [r for g in sampled for r in by_group[g]]
        mrr_prop = sum(1.0 / r[key_proposal] for r in sampled_rows) / len(sampled_rows)
        mrr_base = sum(1.0 / r[key_baseline] for r in sampled_rows) / len(sampled_rows)
        diffs.append(mrr_prop - mrr_base)
    diffs.sort()
    observed = (
        sum(1.0 / r[key_proposal] for r in rows) / len(rows)
        - sum(1.0 / r[key_baseline] for r in rows) / len(rows)
    )
    return {
        "observed_delta_mrr": observed,
        "ci95_low": diffs[int(0.025 * draws)],
        "ci95_high": diffs[int(0.975 * draws)],
    }


def evaluate_x1(records_path: Path, rcsb_cache: Path, seeds: tuple[int, ...] = (7, 17, 31, 43, 97)) -> dict:
    payload = records_path.read_bytes()
    records = [json.loads(line) for line in payload.splitlines() if line.strip()]
    usable = [
        r for r in records
        if r.get("connection_features", {}).get("port_a")
        and r.get("connection_features", {}).get("port_b")
    ]
    geom_by_id = {r["example_id"]: extract_relative_position_geometry(r, rcsb_cache) for r in usable}

    all_pooled_rows = []
    per_seed_summaries = []
    methods = (
        "sequence_only",
        "sequence_plus_local_port",
        "sequence_plus_relative_position_3d",
        "geometry_span_match",
        "nearest_training",
        "length_heuristic",
        "fixed_gs",
    )

    for seed in seeds:
        assignment = grouped_split(usable, seed=seed)
        counts = audit_split(usable, assignment)
        train = [dict(r, split="train") for r in usable if assignment[r["example_id"]] == "train"]
        held_out = [r for r in usable if assignment[r["example_id"]] in ("validation", "blind")]

        w_seq = fit_pairwise_linear(train, geom_by_id, SEQ_FEATURES)
        w_port = fit_pairwise_linear(train, geom_by_id, PORT_FEATURES)
        w_rel3d = fit_pairwise_linear(train, geom_by_id, REL_POS_3D_FEATURES)

        seed_rows = []
        for rec in held_out:
            ref_seq = rec["linker"]["sequence"]
            benchmark = candidate_pool(rec["example_id"], train, reference_sequence=ref_seq)
            ref_id = next(c["candidate_id"] for c in benchmark if c["linker_sequence"] == ref_seq)
            query = {k: v for k, v in rec.items() if k != "linker"}
            geom = geom_by_id[rec["example_id"]]

            row = {
                "seed": seed,
                "example_id": rec["example_id"],
                "group": rec["homology_cluster_id"],
                "candidate_count": len(benchmark),
                "reference_length": len(ref_seq),
                "inter_port_ca_distance_A": round(geom["inter_port_ca_distance_A"], 3),
                "sequence_only": rank_candidates_linear(rec, geom, benchmark, w_seq, SEQ_FEATURES, ref_id),
                "sequence_plus_local_port": rank_candidates_linear(rec, geom, benchmark, w_port, PORT_FEATURES, ref_id),
                "sequence_plus_relative_position_3d": rank_candidates_linear(rec, geom, benchmark, w_rel3d, REL_POS_3D_FEATURES, ref_id),
                "geometry_span_match": rank_geometry_span_match(geom, benchmark, ref_id),
            }
            for base_method in ("nearest_training", "length_heuristic", "fixed_gs"):
                scored = score_candidates(query, benchmark, train, base_method)
                row[base_method] = next(i for i, item in enumerate(scored, 1) if item["candidate_id"] == ref_id)
            seed_rows.append(row)
            all_pooled_rows.append(row)

        seed_metrics = {}
        for m in methods:
            ranks = [r[m] for r in seed_rows]
            seed_metrics[m] = {
                "n": len(ranks),
                "mrr": sum(1.0 / rk for rk in ranks) / len(ranks),
                "top1_count": sum(rk == 1 for rk in ranks),
                "top5_count": sum(rk <= 5 for rk in ranks),
                "top1_rate": sum(rk == 1 for rk in ranks) / len(ranks),
                "top5_rate": sum(rk <= 5 for rk in ranks) / len(ranks),
            }
        per_seed_summaries.append({
            "seed": seed,
            "split_counts": counts,
            "held_out_examples": len(seed_rows),
            "held_out_domain_pair_groups": len({r["group"] for r in seed_rows}),
            "metrics": seed_metrics,
        })

    aggregate = {}
    for m in methods:
        mrrs = [s["metrics"][m]["mrr"] for s in per_seed_summaries]
        top1_counts = sum(s["metrics"][m]["top1_count"] for s in per_seed_summaries)
        top5_counts = sum(s["metrics"][m]["top5_count"] for s in per_seed_summaries)
        total_evals = sum(s["held_out_examples"] for s in per_seed_summaries)
        aggregate[m] = {
            "total_held_out_evaluations": total_evals,
            "mean_mrr_across_seeds": sum(mrrs) / len(mrrs),
            "min_mrr": min(mrrs),
            "max_mrr": max(mrrs),
            "pooled_top1_raw_count": f"{top1_counts}/{total_evals}",
            "pooled_top5_raw_count": f"{top5_counts}/{total_evals}",
            "pooled_top1_rate": top1_counts / total_evals,
            "pooled_top5_rate": top5_counts / total_evals,
        }

    paired_comparisons = {
        "relative_position_3d_vs_sequence_only": cluster_bootstrap_diff(
            all_pooled_rows, "sequence_plus_relative_position_3d", "sequence_only"
        ),
        "local_port_vs_sequence_only": cluster_bootstrap_diff(
            all_pooled_rows, "sequence_plus_local_port", "sequence_only"
        ),
        "relative_position_3d_vs_local_port": cluster_bootstrap_diff(
            all_pooled_rows, "sequence_plus_relative_position_3d", "sequence_plus_local_port"
        ),
    }

    return {
        "experiment": "X1_held_out_domain_pair_relative_position_transfer",
        "status": "completed_held_out_domain_pair_ablation",
        "source_records_sha256": sha256_bytes(payload),
        "total_records": len(records),
        "coordinate_complete_records": len(usable),
        "unique_domain_pair_groups": len({r["homology_cluster_id"] for r in usable}),
        "seeds": list(seeds),
        "aggregate_metrics": aggregate,
        "paired_cluster_bootstrap_mrr_differences": paired_comparisons,
        "per_seed": per_seed_summaries,
        "limitations": [
            f"Evaluated on {len(usable)} coordinate-complete CATH v4.4.0 two-domain natural proteins across {len({r['homology_cluster_id'] for r in usable})} CATH domain-pair groups.",
            "No Li et al. 2026 figures were used; features were computed directly from verified PDB coordinates.",
            "Pilot records have unreviewed linker interiors and test structural-context transfer for natural-linker ranking, not wet-lab expression.",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--rcsb-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = evaluate_x1(args.records, args.rcsb_cache)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in result.items() if k != "per_seed"}, indent=2, sort_keys=True))
