"""Non-circular evaluation of the Part 4 pipeline.

Everything here uses only quantities recorded in the dataset:

* `scorecard_on_real_linkers` -- how often the deterministic rules reject *real* linkers
  (sensitivity check).  Geometry rules run only for the coordinate-complete CATH tier.
* `corruption_check` -- rejection / ranking of the six constructed corruptions.  These are built
  to violate the same properties the rules and several features measure, so this is a
  by-construction sanity check and is labelled as such.
* `decoy_ranking_benchmark` -- can a method pick the held-out real linker out of length- and
  rigidity-matched real linkers from other clusters?  Chance is 1/(K+1).
* `generator_comparison_on_cath` -- six design strategies compared on the 91 coordinate-complete
  natural linkers using the real inter-port distance, with no oracle access to the true length.

Confidence intervals are cluster bootstraps over `homology_cluster_id`.
"""

from __future__ import annotations

import math
import random
from collections import Counter, defaultdict
from typing import Any, Callable, Mapping, Sequence

from linkerforge.design.engine import flanking_context_similarity
from linkerforge.geometry.constants import DEFAULT_GS_COMPOSITION
from linkerforge.geometry.span import comfortable_residues, minimum_residues, required_length
from linkerforge.retrieval.baseline import assign_cluster_disjoint_folds

from .ranker import (
    SYNTHETIC_NEGATIVE_TYPES,
    AuthoritySequenceIndex,
    evaluate_hard_rejection_scorecard,
    fit_pairwise_compatibility_ranker,
    generate_synthetic_negatives,
    record_query_context,
)

AA = "ACDEFGHIKLMNPQRSTVWY"
RULES = (
    "physical_span_closure_pass",
    "domain_clash_clearance_pass",
    "anchor_accessibility_pass",
    "hydrophobic_stretch_pass",
    "prohibited_cleavage_pass",
    "composite_rigid_cleavable_pass",
)
GENERATORS = (
    "direct_0aa_fusion",
    "fixed_ggggs3",
    "random_length_matched_gs",
    "geometry_only",
    "nearest_neighbour_retrieval",
    "linkerforge_ranker",
)


def _comp(seq: str) -> list[float]:
    n = len(seq) or 1
    return [seq.count(a) / n for a in AA]


def _cos(a: Sequence[float], b: Sequence[float]) -> float:
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(x * x for x in b))
    return 0.0 if na == 0 or nb == 0 else sum(x * y for x, y in zip(a, b)) / (na * nb)


def _cluster_id(rec: Mapping[str, Any]) -> str:
    return str(rec.get("homology_cluster_id") or rec.get("record_id"))


def cluster_bootstrap_mean(
    values: Sequence[float], clusters: Sequence[str], *, reps: int = 500, seed: int = 0
) -> list[float]:
    sums: dict[str, float] = defaultdict(float)
    counts: dict[str, int] = defaultdict(int)
    for v, c in zip(values, clusters):
        sums[c] += v
        counts[c] += 1
    ids = sorted(sums)
    if not ids:
        return [float("nan"), float("nan")]
    rng = random.Random(seed)
    stats = []
    for _ in range(reps):
        picked = [rng.choice(ids) for _ in ids]
        stats.append(sum(sums[c] for c in picked) / sum(counts[c] for c in picked))
    stats.sort()
    return [round(stats[int(0.025 * reps)], 4), round(stats[int(0.975 * reps) - 1], 4)]


def _mean_ci(values: Sequence[float], clusters: Sequence[str], seed: int = 0) -> dict[str, Any]:
    if not values:
        return {"n": 0, "mean": None, "ci95": None}
    return {
        "n": len(values),
        "mean": round(sum(values) / len(values), 4),
        "ci95": cluster_bootstrap_mean(values, clusters, seed=seed),
    }


# --------------------------------------------------------------------------------------
# 1. rule sensitivity on real linkers
# --------------------------------------------------------------------------------------
def scorecard_on_real_linkers(records_by_tier: Mapping[str, Sequence[Mapping[str, Any]]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for tier, recs in records_by_tier.items():
        fails: Counter[str] = Counter()
        n = passes = geom_n = dataset_flag_ok = dataset_flag_n = 0
        for rec in recs:
            ctx = record_query_context(rec)
            if not ctx["sequence"]:
                continue
            card = evaluate_hard_rejection_scorecard(
                ctx["sequence"],
                required_path_length_A=ctx["path_A"],
                target_length=ctx["target_length"],
                port_a_sasa=ctx["sasa_a"],
                port_b_sasa=ctx["sasa_b"],
                requested_rigidity=ctx["rigidity"],
                requested_cleavable=False,
            )
            n += 1
            passes += int(card.hard_rejection_pass)
            for rule in RULES:
                if not getattr(card, rule):
                    fails[rule] += 1
            if ctx["path_A"] is not None:
                geom_n += 1
            flag = (rec.get("geometry") or {}).get("physical_3p8A_feasible")
            if flag is not None:
                dataset_flag_n += 1
                dataset_flag_ok += int(bool(flag))
        out[tier] = {
            "n_records": n,
            "n_with_real_geometry": geom_n,
            "all_rules_pass": passes,
            "all_rules_pass_rate": round(passes / n, 4) if n else None,
            "rule_failures": {r: fails.get(r, 0) for r in RULES},
            "dataset_own_physical_3p8A_flag_true": dataset_flag_ok if dataset_flag_n else None,
            "dataset_own_physical_3p8A_flag_n": dataset_flag_n or None,
        }
    return out


# --------------------------------------------------------------------------------------
# 2. constructed corruptions (by construction)
# --------------------------------------------------------------------------------------
def corruption_check(
    records: Sequence[dict[str, Any]], *, n_folds: int = 5, seed: int = 23
) -> dict[str, Any]:
    rejected: Counter[str] = Counter()
    total = 0
    for rec in records:
        ctx = record_query_context(rec)
        if len(ctx["sequence"]) < 4:
            continue
        total += 1
        negs = generate_synthetic_negatives(
            ctx["sequence"],
            target_length=ctx["target_length"],
            requested_rigidity=ctx["rigidity"],
            left_flank=ctx["left_flank"],
            right_flank=ctx["right_flank"],
        )
        for name, seq in negs.items():
            card = evaluate_hard_rejection_scorecard(
                seq,
                required_path_length_A=ctx["path_A"],
                target_length=ctx["target_length"],
                port_a_sasa=ctx["sasa_a"],
                port_b_sasa=ctx["sasa_b"],
                requested_rigidity=ctx["rigidity"],
            )
            rejected[name] += int(not card.hard_rejection_pass)

    wins: dict[str, list[float]] = {c: [] for c in SYNTHETIC_NEGATIVE_TYPES}
    clusters: dict[str, list[str]] = {c: [] for c in SYNTHETIC_NEGATIVE_TYPES}
    for train_rows, test_rows in assign_cluster_disjoint_folds(records, n_folds=n_folds, seed=seed):
        ranker = fit_pairwise_compatibility_ranker(train_rows, seed=seed)
        for rec in test_rows:
            ctx = record_query_context(rec)
            if len(ctx["sequence"]) < 4:
                continue
            kw = dict(
                target_length=ctx["target_length"],
                required_path_length_A=ctx["path_A"],
                requested_rigidity=ctx["rigidity"],
                left_flank=ctx["left_flank"],
                right_flank=ctx["right_flank"],
                port_a_sasa=ctx["sasa_a"],
                port_b_sasa=ctx["sasa_b"],
            )
            pos, _, _ = ranker.score_candidate(ctx["sequence"], **kw)
            negs = generate_synthetic_negatives(
                ctx["sequence"],
                target_length=ctx["target_length"],
                requested_rigidity=ctx["rigidity"],
                left_flank=ctx["left_flank"],
                right_flank=ctx["right_flank"],
            )
            for name, seq in negs.items():
                s, _, _ = ranker.score_candidate(seq, **kw)
                wins[name].append(1.0 if pos > s else 0.0)
                clusters[name].append(_cluster_id(rec))
    return {
        "status": "by_construction_sanity_check_not_evidence_of_fusion_quality",
        "n_records": total,
        "scorecard_rejection_rate_per_corruption": {
            c: round(rejected[c] / total, 4) if total else None for c in SYNTHETIC_NEGATIVE_TYPES
        },
        "ranker_heldout_pairwise_accuracy_per_corruption": {
            c: _mean_ci(wins[c], clusters[c]) for c in SYNTHETIC_NEGATIVE_TYPES
        },
    }


# --------------------------------------------------------------------------------------
# 3. decoy ranking
# --------------------------------------------------------------------------------------
def _mid_rank(true_score: float, other_scores: Sequence[float]) -> float:
    greater = sum(1 for s in other_scores if s > true_score)
    equal = sum(1 for s in other_scores if s == true_score)
    return 1.0 + greater + 0.5 * equal


def decoy_ranking_benchmark(
    records: Sequence[dict[str, Any]],
    decoy_pool: Sequence[dict[str, Any]],
    *,
    n_folds: int = 5,
    seed: int = 23,
    n_decoys: int = 19,
    length_tolerance: int = 2,
    min_decoys: int = 9,
    flank_neighbours: int = 10,
) -> dict[str, Any]:
    arms = ("random", "training_composition_prior", "flank_conditioned_nn", "linkerforge_ranker")
    pool_idx: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for r in decoy_pool:
        c = record_query_context(r)
        if len(c["sequence"]) >= 4:
            pool_idx[(c["rigidity"], len(c["sequence"]))].append(r)

    rng = random.Random(seed)
    ranks: dict[str, list[float]] = {a: [] for a in arms}
    sizes: list[int] = []
    clusters: list[str] = []
    skipped = 0

    for train_rows, test_rows in assign_cluster_disjoint_folds(records, n_folds=n_folds, seed=seed):
        ranker = fit_pairwise_compatibility_ranker(train_rows, seed=seed)
        train_ctx = []
        by_rig: dict[str, list[list[float]]] = defaultdict(list)
        for r in train_rows:
            c = record_query_context(r)
            if len(c["sequence"]) < 4:
                continue
            c["comp"] = _comp(c["sequence"])
            train_ctx.append(c)
            by_rig[c["rigidity"]].append(c["comp"])
        mean_comp = {
            k: [sum(col) / len(v) for col in zip(*v)] for k, v in by_rig.items() if v
        }
        for rec in test_rows:
            ctx = record_query_context(rec)
            true_seq = ctx["sequence"]
            if len(true_seq) < 4:
                continue
            cid = _cluster_id(rec)
            cands = [
                d
                for dl in range(-length_tolerance, length_tolerance + 1)
                for d in pool_idx.get((ctx["rigidity"], len(true_seq) + dl), [])
                if _cluster_id(d) != cid and d["linker_sequence"] != true_seq
            ]
            if len(cands) < min_decoys:
                skipped += 1
                continue
            decoys = rng.sample(cands, min(n_decoys, len(cands)))
            seqs = [true_seq] + [str(d["linker_sequence"]).upper() for d in decoys]

            sims = sorted(
                (
                    (
                        flanking_context_similarity(
                            ctx["left_flank"], ctx["right_flank"], None, None,
                            {"left_flank_sequence": t["left_flank"], "right_flank_sequence": t["right_flank"]},
                        ),
                        t["comp"],
                    )
                    for t in train_ctx
                ),
                key=lambda x: -x[0],
            )[:flank_neighbours]
            nn_comp = [sum(col) / len(sims) for col in zip(*[c for _, c in sims])] if sims else None
            prior = mean_comp.get(ctx["rigidity"]) or next(iter(mean_comp.values()), None)

            scores: dict[str, list[float]] = {a: [] for a in arms}
            for s in seqs:
                comp = _comp(s)
                scores["random"].append(rng.random())
                scores["training_composition_prior"].append(_cos(comp, prior) if prior else 0.0)
                scores["flank_conditioned_nn"].append(_cos(comp, nn_comp) if nn_comp else 0.0)
                sc, _, _ = ranker.score_candidate(
                    s,
                    target_length=ctx["target_length"],
                    required_path_length_A=ctx["path_A"],
                    requested_rigidity=ctx["rigidity"],
                    left_flank=ctx["left_flank"],
                    right_flank=ctx["right_flank"],
                    port_a_sasa=ctx["sasa_a"],
                    port_b_sasa=ctx["sasa_b"],
                )
                scores["linkerforge_ranker"].append(sc)
            for a in arms:
                ranks[a].append(_mid_rank(scores[a][0], scores[a][1:]))
            sizes.append(len(seqs))
            clusters.append(cid)

    def chance(fn: Callable[[int], float]) -> float:
        return round(sum(fn(n) for n in sizes) / len(sizes), 4) if sizes else float("nan")

    out_arms = {}
    for a in arms:
        rr = ranks[a]
        out_arms[a] = {
            "hit_at_1": _mean_ci([1.0 if r <= 1.0 else 0.0 for r in rr], clusters, seed),
            "hit_at_5": _mean_ci([1.0 if r <= 5.0 else 0.0 for r in rr], clusters, seed),
            "mrr": _mean_ci([1.0 / r for r in rr], clusters, seed),
        }
    return {
        "task": "true_linker_vs_length_and_rigidity_matched_real_linkers_from_other_clusters",
        "n_queries": len(sizes),
        "n_skipped_too_few_decoys": skipped,
        "mean_candidates_per_query": round(sum(sizes) / len(sizes), 2) if sizes else None,
        "chance": {
            "hit_at_1": chance(lambda n: 1.0 / n),
            "hit_at_5": chance(lambda n: min(5, n) / n),
            "mrr": chance(lambda n: sum(1.0 / i for i in range(1, n + 1)) / n),
        },
        "arms": out_arms,
    }


# --------------------------------------------------------------------------------------
# 4. design-generator comparison on coordinate-complete natural linkers
# --------------------------------------------------------------------------------------
LENGTH_POLICIES = ("unstrained_reach", "design_minimum_plus_2")


def geometry_length(distance_A: float, policy: str = "unstrained_reach") -> int:
    """Linker length proposed from geometry alone (no access to the true length).

    `unstrained_reach`: fewest residues whose Gaussian-chain RMS end-to-end distance reaches the
    span (G4S2 statistics).  `design_minimum_plus_2`: the 3.5 A/residue design minimum plus two.
    Both policies were fixed before looking at any result and both are reported.
    """
    if policy == "unstrained_reach":
        comfy = comfortable_residues(distance_A, DEFAULT_GS_COMPOSITION)
        base = comfy if comfy is not None else minimum_residues(distance_A) + 2
    elif policy == "design_minimum_plus_2":
        base = int(required_length(distance_A, linker_type="flexible")["minimum_residues_design"]) + 2
    else:
        raise ValueError(f"policy must be one of {LENGTH_POLICIES}")
    return max(5, int(base))


def _repeat_to_length(unit: str, n: int) -> str:
    return (unit * (n // len(unit) + 1))[:n] if n > 0 else ""


def generator_comparison_on_cath(
    natural_records: Sequence[dict[str, Any]],
    *,
    n_folds: int = 5,
    seed: int = 23,
    length_policy: str = "unstrained_reach",
) -> dict[str, Any]:
    metrics: dict[str, dict[str, list[float]]] = {
        g: defaultdict(list) for g in GENERATORS
    }
    clusters: list[str] = []
    rng = random.Random(seed)
    for train_rows, test_rows in assign_cluster_disjoint_folds(natural_records, n_folds=n_folds, seed=seed):
        index = AuthoritySequenceIndex(train_rows)
        ranker = fit_pairwise_compatibility_ranker(train_rows, seed=seed)
        for rec in test_rows:
            ctx = record_query_context(rec)
            if ctx["path_A"] is None or len(ctx["sequence"]) < 4:
                continue
            d = float(ctx["path_A"])
            true_seq = ctx["sequence"]
            l_geom = geometry_length(d, length_policy)
            clusters.append(_cluster_id(rec))
            geo_seq = _repeat_to_length("GGGGS", l_geom)
            hits = index.retrieve(
                target_length=l_geom,
                rigidity=ctx["rigidity"],
                left_flank_sequence=ctx["left_flank"],
                right_flank_sequence=ctx["right_flank"],
                top_k=4,
            )
            nn_seq = hits[0].sequence if hits else geo_seq
            pool = [geo_seq] + [h.sequence for h in hits]

            def _rank_key(seq: str) -> tuple[int, float]:
                sc, card, _ = ranker.score_candidate(
                    seq,
                    target_length=l_geom,
                    required_path_length_A=d,
                    requested_rigidity=ctx["rigidity"],
                    left_flank=ctx["left_flank"],
                    right_flank=ctx["right_flank"],
                    port_a_sasa=ctx["sasa_a"],
                    port_b_sasa=ctx["sasa_b"],
                )
                return (int(card.hard_rejection_pass), sc)

            chosen = {
                "direct_0aa_fusion": "",
                "fixed_ggggs3": "GGGGSGGGGSGGGGS",
                "random_length_matched_gs": "".join(rng.choice("GS") for _ in range(l_geom)),
                "geometry_only": geo_seq,
                "nearest_neighbour_retrieval": nn_seq,
                "linkerforge_ranker": max(pool, key=_rank_key),
            }
            true_comp = _comp(true_seq)
            for g, seq in chosen.items():
                card = evaluate_hard_rejection_scorecard(
                    seq,
                    required_path_length_A=d,
                    target_length=l_geom,
                    port_a_sasa=ctx["sasa_a"],
                    port_b_sasa=ctx["sasa_b"],
                    requested_rigidity=ctx["rigidity"],
                )
                m = metrics[g]
                m["span_closure_pass"].append(float(card.physical_span_closure_pass))
                m["all_rules_pass"].append(float(card.hard_rejection_pass))
                m["abs_length_error_vs_natural"].append(float(abs(len(seq) - len(true_seq))))
                m["length_within_2aa_of_natural"].append(float(abs(len(seq) - len(true_seq)) <= 2))
                m["composition_cosine_to_natural"].append(_cos(_comp(seq), true_comp) if seq else 0.0)
    return {
        "task": "design_strategies_on_cath_records_with_real_inter_port_distance",
        "n_queries": len(clusters),
        "length_policy": length_policy,
        "methods": {
            g: {name: _mean_ci(vals, clusters, seed) for name, vals in m.items()}
            for g, m in metrics.items()
        },
    }


