"""Domain-Pair and Family-Disjoint Retrieval Baseline (`X2`) Implementation.

Evaluates linker reconstruction on held-out `homology_cluster_id` folds using
flanking domain context only (and 3D port distance on coordinate-resolved CATH
records), comparing:
  1. `fixed_ggggs3_control`: Fixed 15-aa `(GGGGS)3` control (`flex-ggggs3`).
  2. `training_median_length_gs_control`: Non-retrieval baseline using a `(GGGGS)`
     repeat truncated to the training-fold median linker length.
  3. `flanking_context_template_retrieval`: Predicts target length and rigidity
     from k=3 nearest training-fold flanking contexts and retrieves the closest
     template from `literature_engineered_template` (`linker_property_library.v1`).
  4. `flanking_context_nn_natural_authority`: Nearest-neighbour retrieval from
     the training fold of verified natural authorities (`natural_cath_pdb_verified`
     / `natural_disprot_curated_dfl`) using flanking sequence and domain context only.
  5. `flanking_context_nn_natural_plus_mined`: Nearest-neighbour retrieval from
     the training fold augmented with benchmark/cluster-disjoint mined inter-domain
     candidates (`mined_high_confidence_dfl_call` + `mined_provisional_dfl_candidate`).
  6. `geometry_only_control` (3D cohort): Uses `inter_port_ca_distance_A`
     (`comfortable_residues` / 2.5 A/res span rule) without flanking sequence identity.
  7. `geometry_plus_flanking_hybrid_retrieval` (3D cohort): Combines physical 3.8 A
     bridgeability and span proximity with flanking-context similarity.
  8. `length_matched_template_control` (Oracle-length control): Closest-length
     template from `literature_engineered_template` given the target length `L_true`.
  9. `length_matched_frequent_natural_control` (Oracle-length non-retrieval control):
     Most frequent amino-acid repeat from the training fold at length `L_true`
     (`CHECKLIST.md` item `X2`).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
import hashlib
from math import sqrt
import random
from typing import Any, Sequence

from linkerforge.geometry.constants import CA_VIRTUAL_BOND_ANGSTROM
from linkerforge.geometry.span import comfortable_residues, minimum_residues
from linkerforge.properties import featurise

try:
    import parasail  # type: ignore
except ImportError:  # pragma: no cover
    parasail = None

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def _digest_key(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


@lru_cache(maxsize=65536)
def sequence_nw_identity(seq_a: str, seq_b: str) -> float:
    """Global Needleman-Wunsch sequence identity normalized by max(len(a), len(b))."""
    a = seq_a.strip().upper()
    b = seq_b.strip().upper()
    if not a or not b:
        return 0.0
    denom = float(max(len(a), len(b)))
    if parasail is not None:
        res = parasail.nw_stats_striped_16(a, b, 10, 1, parasail.blosum62)
        return float(res.matches) / denom
    # Pure-Python NW fallback
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        ca = a[i - 1]
        row_prev = dp[i - 1]
        row_cur = dp[i]
        for j in range(1, m + 1):
            match = row_prev[j - 1] + (1 if ca == b[j - 1] else 0)
            delete = row_prev[j]
            insert = row_cur[j - 1]
            row_cur[j] = max(match, delete, insert)
    return float(dp[n][m]) / denom


def positional_aa_recovery(pred_seq: str, true_seq: str) -> float:
    """Exact position-wise amino-acid recovery divided by len(true_seq)."""
    if not true_seq:
        return 0.0
    matches = sum(1 for a, b in zip(pred_seq, true_seq) if a == b)
    return float(matches) / float(len(true_seq))


def _cosine(vec_a: Sequence[float], vec_b: Sequence[float]) -> float:
    dot = sum(a * b for a, b in zip(vec_a, vec_b))
    na = sqrt(sum(a * a for a in vec_a))
    nb = sqrt(sum(b * b for b in vec_b))
    if na <= 0.0 or nb <= 0.0:
        return 0.0
    return max(0.0, min(1.0, dot / (na * nb)))


def composition_cosine_similarity(seq_a: str, seq_b: str) -> float:
    """Cosine similarity on the 6-feature Chen et al. (2013) composition vector."""
    return _cosine(featurise(seq_a), featurise(seq_b))


@lru_cache(maxsize=65536)
def _flank_feature_vector(seq: str) -> tuple[float, ...]:
    clean = [ch for ch in seq.upper() if ch in ALPHABET]
    if not clean:
        return tuple(0.0 for _ in range(24))
    n = float(len(clean))
    aa_freqs = [clean.count(aa) / n for aa in ALPHABET]
    charged = sum(clean.count(aa) for aa in "DEKR") / n
    polar = sum(clean.count(aa) for aa in "STNQYH") / n
    hydrophobic = sum(clean.count(aa) for aa in "AVILMFW") / n
    gly_pro = sum(clean.count(aa) for aa in "GP") / n
    return tuple(aa_freqs + [charged, polar, hydrophobic, gly_pro])


def _domain_prefix_similarity(dom_a: str | None, dom_b: str | None) -> float:
    if not dom_a or not dom_b:
        return 0.0
    if dom_a == dom_b:
        return 1.0
    # For CATH superfamilies (e.g., 3.30.40.10), compare class/architecture/topology prefixes
    if "." in dom_a and "." in dom_b:
        parts_a = dom_a.split(".")
        parts_b = dom_b.split(".")
        shared = 0
        for pa, pb in zip(parts_a, parts_b):
            if pa == pb:
                shared += 1
            else:
                break
        return shared / 4.0
    return 0.0


def flanking_pair_similarity(flank_a: dict[str, Any], flank_b: dict[str, Any]) -> float:
    """Compute flanking-context-only similarity between two records.

    Never reads or touches the linker sequence or linker length. Uses:
      - Boundary 25-aa N-flank and C-flank amino-acid & biophysical composition cosine
      - Boundary 20-aa N-flank and C-flank sequence alignment identity
      - Left and right domain family / CATH hierarchy similarity
    """
    la = str(flank_a.get("left_flank_sequence") or "")
    ra = str(flank_a.get("right_flank_sequence") or "")
    lb = str(flank_b.get("left_flank_sequence") or "")
    rb = str(flank_b.get("right_flank_sequence") or "")

    left_comp = _cosine(_flank_feature_vector(la[-25:]), _flank_feature_vector(lb[-25:])) if (la and lb) else 0.0
    right_comp = _cosine(_flank_feature_vector(ra[:25]), _flank_feature_vector(rb[:25])) if (ra and rb) else 0.0
    left_id = sequence_nw_identity(la[-20:], lb[-20:]) if (la and lb) else 0.0
    right_id = sequence_nw_identity(ra[:20], rb[:20]) if (ra and rb) else 0.0

    seq_comp = 0.5 * (left_comp + right_comp)
    seq_id = 0.5 * (left_id + right_id)

    src_a = flank_a.get("left_domain_source") or flank_a.get("right_domain_source")
    src_b = flank_b.get("left_domain_source") or flank_b.get("right_domain_source")
    if src_a and src_b and src_a == src_b:
        dom_left = _domain_prefix_similarity(flank_a.get("left_domain_id"), flank_b.get("left_domain_id"))
        dom_right = _domain_prefix_similarity(flank_a.get("right_domain_id"), flank_b.get("right_domain_id"))
        swap_dom = 0.5 * (
            _domain_prefix_similarity(flank_a.get("left_domain_id"), flank_b.get("right_domain_id"))
            + _domain_prefix_similarity(flank_a.get("right_domain_id"), flank_b.get("left_domain_id"))
        )
        dom_score = max(0.5 * (dom_left + dom_right), 0.7 * swap_dom)
        return 0.45 * seq_comp + 0.35 * seq_id + 0.20 * dom_score

    return 0.5625 * seq_comp + 0.4375 * seq_id


def _fast_flank_comp_and_domain(flank_a: dict[str, Any], flank_b: dict[str, Any]) -> float:
    la = str(flank_a.get("left_flank_sequence") or "")
    ra = str(flank_a.get("right_flank_sequence") or "")
    lb = str(flank_b.get("left_flank_sequence") or "")
    rb = str(flank_b.get("right_flank_sequence") or "")
    left_comp = _cosine(_flank_feature_vector(la[-25:]), _flank_feature_vector(lb[-25:])) if (la and lb) else 0.0
    right_comp = _cosine(_flank_feature_vector(ra[:25]), _flank_feature_vector(rb[:25])) if (ra and rb) else 0.0
    seq_comp = 0.5 * (left_comp + right_comp)
    src_a = flank_a.get("left_domain_source") or flank_a.get("right_domain_source")
    src_b = flank_b.get("left_domain_source") or flank_b.get("right_domain_source")
    if src_a and src_b and src_a == src_b:
        dom_left = _domain_prefix_similarity(flank_a.get("left_domain_id"), flank_b.get("left_domain_id"))
        dom_right = _domain_prefix_similarity(flank_a.get("right_domain_id"), flank_b.get("right_domain_id"))
        swap_dom = 0.5 * (
            _domain_prefix_similarity(flank_a.get("left_domain_id"), flank_b.get("right_domain_id"))
            + _domain_prefix_similarity(flank_a.get("right_domain_id"), flank_b.get("left_domain_id"))
        )
        dom_score = max(0.5 * (dom_left + dom_right), 0.7 * swap_dom)
        return 0.70 * seq_comp + 0.30 * dom_score
    return 0.70 * seq_comp


def _shortlist_by_flank_composition(
    query_flank: dict[str, Any],
    pool: Sequence[dict[str, Any]],
    *,
    limit: int = 250,
) -> Sequence[dict[str, Any]]:
    if len(pool) <= limit:
        return pool
    scored = sorted(
        pool,
        key=lambda cand: (
            -_fast_flank_comp_and_domain(query_flank, cand.get("flanking_context") or {}),
            cand["record_id"],
        ),
    )
    return scored[:limit]


def assign_cluster_disjoint_folds(
    records: Sequence[dict[str, Any]],
    *,
    n_folds: int = 5,
    seed: int = 17,
) -> list[tuple[list[dict[str, Any]], list[dict[str, Any]]]]:
    """Partition records into `n_folds` strictly disjoint `homology_cluster_id` folds."""
    by_cluster: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in records:
        by_cluster[str(r["homology_cluster_id"])].append(r)

    ordered_clusters = sorted(
        by_cluster.keys(),
        key=lambda cid: (-len(by_cluster[cid]), _digest_key(f"{seed}:{cid}")),
    )
    fold_bins: list[list[dict[str, Any]]] = [[] for _ in range(n_folds)]
    for cid in ordered_clusters:
        # Greedily place cluster in the currently smallest fold bin
        target_idx = min(range(n_folds), key=lambda idx: (len(fold_bins[idx]), idx))
        fold_bins[target_idx].extend(by_cluster[cid])

    folds: list[tuple[list[dict[str, Any]], list[dict[str, Any]]]] = []
    for test_idx in range(n_folds):
        test_rows = fold_bins[test_idx]
        train_rows = [r for idx, b in enumerate(fold_bins) if idx != test_idx for r in b]
        train_clusters = {r["homology_cluster_id"] for r in train_rows}
        test_clusters = {r["homology_cluster_id"] for r in test_rows}
        if train_clusters & test_clusters:
            raise AssertionError("Cluster leakage detected across train/test folds")
        folds.append((train_rows, test_rows))
    return folds


def _repeat_to_length(unit: str, length: int) -> str:
    reps = (length // len(unit)) + 2
    return (unit * reps)[:length]


def _most_frequent_dipeptide_repeat(train_rows: Sequence[dict[str, Any]], length: int) -> str:
    counts: Counter[str] = Counter()
    for r in train_rows:
        seq = str(r["linker_sequence"])
        for i in range(len(seq) - 1):
            counts[seq[i : i + 2]] += 1
    best_pair = counts.most_common(1)[0][0] if counts else "GS"
    return _repeat_to_length(best_pair, length)


def _select_template_by_length_and_rigidity(
    templates: Sequence[dict[str, Any]],
    target_len: float,
    target_rigidity: str,
) -> dict[str, Any]:
    matching = [
        t
        for t in templates
        if t["biophysical_properties"].get("rigidity_predicted_label") == target_rigidity
    ]
    pool = matching if matching else list(templates)
    return min(
        pool,
        key=lambda t: (abs(int(t["linker_length"]) - target_len), int(t["linker_length"]), t["record_id"]),
    )


def _score_prediction(
    pred_seq: str,
    pred_rigidity: str,
    true_row: dict[str, Any],
) -> dict[str, Any]:
    true_seq = str(true_row["linker_sequence"])
    true_len = int(true_row["linker_length"])
    pred_len = len(pred_seq)
    len_err = abs(pred_len - true_len)
    true_rig = str(true_row["biophysical_properties"]["rigidity_predicted_label"])

    geom = true_row.get("geometry")
    if geom and geom.get("inter_port_ca_distance_A") is not None:
        ca_dist = float(geom["inter_port_ca_distance_A"])
        feasible_3p8 = ca_dist <= (pred_len + 1) * CA_VIRTUAL_BOND_ANGSTROM
    else:
        feasible_3p8 = None

    return {
        "record_id": true_row["record_id"],
        "homology_cluster_id": true_row["homology_cluster_id"],
        "true_length": true_len,
        "pred_length": pred_len,
        "length_abs_error": float(len_err),
        "length_within_2aa": int(len_err <= 2),
        "length_within_5aa": int(len_err <= 5),
        "sequence_nw_identity": sequence_nw_identity(pred_seq, true_seq),
        "positional_aa_recovery": positional_aa_recovery(pred_seq, true_seq),
        "composition_cosine_similarity": composition_cosine_similarity(pred_seq, true_seq),
        "rigidity_match": int(pred_rigidity == true_rig),
        "physical_3p8A_feasible": None if feasible_3p8 is None else int(feasible_3p8),
    }


def _aggregate_arm(
    per_example: Sequence[dict[str, Any]],
    *,
    bootstrap_iters: int = 1000,
    seed: int = 17,
) -> dict[str, Any]:
    n = len(per_example)
    if n == 0:
        raise ValueError("Empty evaluation list")
    mean_mae = sum(x["length_abs_error"] for x in per_example) / n
    within_2 = sum(x["length_within_2aa"] for x in per_example)
    within_5 = sum(x["length_within_5aa"] for x in per_example)
    mean_nw = sum(x["sequence_nw_identity"] for x in per_example) / n
    mean_pos = sum(x["positional_aa_recovery"] for x in per_example) / n
    mean_comp = sum(x["composition_cosine_similarity"] for x in per_example) / n
    rig_matches = sum(x["rigidity_match"] for x in per_example)

    geom_items = [x["physical_3p8A_feasible"] for x in per_example if x["physical_3p8A_feasible"] is not None]
    geom_feasible_count = sum(int(v) for v in geom_items) if geom_items else None
    geom_feasible_rate = (geom_feasible_count / len(geom_items)) if geom_items else None

    # Cluster-level bootstrap 95% CI on length_mae and composition_cosine_similarity
    by_cluster: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for x in per_example:
        by_cluster[str(x["homology_cluster_id"])].append(x)
    clusters = sorted(by_cluster.keys())
    rng = random.Random(seed)
    boot_mae: list[float] = []
    boot_nw: list[float] = []
    boot_comp: list[float] = []
    for _ in range(bootstrap_iters):
        sampled = [ rng.choice(clusters) for _ in range(len(clusters)) ]
        rows = [item for cid in sampled for item in by_cluster[cid]]
        boot_mae.append(sum(r["length_abs_error"] for r in rows) / len(rows))
        boot_nw.append(sum(r["sequence_nw_identity"] for r in rows) / len(rows))
        boot_comp.append(sum(r["composition_cosine_similarity"] for r in rows) / len(rows))
    boot_mae.sort()
    boot_nw.sort()
    boot_comp.sort()
    lo_idx = int(0.025 * bootstrap_iters)
    hi_idx = int(0.975 * bootstrap_iters)

    return {
        "evaluated_records": n,
        "distinct_homology_clusters": len(clusters),
        "length_mae_residues": round(mean_mae, 4),
        "length_mae_cluster_boot_95ci": [round(boot_mae[lo_idx], 4), round(boot_mae[hi_idx], 4)],
        "length_within_2aa_count": within_2,
        "length_within_2aa_rate": round(within_2 / n, 4),
        "length_within_5aa_count": within_5,
        "length_within_5aa_rate": round(within_5 / n, 4),
        "sequence_nw_identity_mean": round(mean_nw, 4),
        "sequence_nw_identity_cluster_boot_95ci": [round(boot_nw[lo_idx], 4), round(boot_nw[hi_idx], 4)],
        "positional_aa_recovery_mean": round(mean_pos, 4),
        "composition_cosine_similarity_mean": round(mean_comp, 4),
        "composition_cosine_cluster_boot_95ci": [round(boot_comp[lo_idx], 4), round(boot_comp[hi_idx], 4)],
        "rigidity_axis_match_count": rig_matches,
        "rigidity_axis_match_rate": round(rig_matches / n, 4),
        "physical_3p8A_feasible_count": geom_feasible_count,
        "physical_3p8A_feasible_rate": None if geom_feasible_rate is None else round(geom_feasible_rate, 4),
    }


def evaluate_x2_retrieval_benchmark(
    natural_records: Sequence[dict[str, Any]],
    mined_records: Sequence[dict[str, Any]],
    template_records: Sequence[dict[str, Any]],
    *,
    has_3d_geometry: bool,
    n_folds: int = 5,
    seed: int = 17,
) -> dict[str, Any]:
    """Run the complete X2 held-out cluster-disjoint retrieval benchmark."""
    folds = assign_cluster_disjoint_folds(natural_records, n_folds=n_folds, seed=seed)
    arm_examples: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for train_rows, test_rows in folds:
        test_clusters = {r["homology_cluster_id"] for r in test_rows}
        # Disjoint mined pool: exclude any mined record sharing a homology_cluster_id with test fold
        disjoint_mined = [
            m for m in mined_records if m["homology_cluster_id"] not in test_clusters
        ]
        combined_pool = list(train_rows) + disjoint_mined

        # Training median length for blind non-retrieval length control
        sorted_lens = sorted(int(r["linker_length"]) for r in train_rows)
        median_len = sorted_lens[len(sorted_lens) // 2]
        median_gs_seq = _repeat_to_length("GGGGS", median_len)

        for test_row in test_rows:
            query_flank = test_row.get("flanking_context") or {}
            true_len = int(test_row["linker_length"])

            # 1. Fixed (GGGGS)3 control (15 aa)
            arm_examples["fixed_ggggs3_control"].append(
                _score_prediction("GGGGSGGGGSGGGGS", "flexible", test_row)
            )

            # 2. Training median-length GS control (non-retrieval)
            arm_examples["training_median_length_gs_control"].append(
                _score_prediction(median_gs_seq, "flexible", test_row)
            )

            # Rank training natural authorities by flanking-context similarity ONLY
            ranked_train = sorted(
                train_rows,
                key=lambda cand: (
                    -flanking_pair_similarity(query_flank, cand.get("flanking_context") or {}),
                    cand["record_id"],
                ),
            )
            best_nat = ranked_train[0]
            top3_nat = ranked_train[:3]

            # 3. Flanking-context Template Retrieval (PR #2 template library from flanking context only)
            pred_len_from_flank = sum(int(r["linker_length"]) for r in top3_nat) / len(top3_nat)
            rig_votes = [
                str(r["biophysical_properties"]["rigidity_predicted_label"]) for r in top3_nat
            ]
            pred_rig_from_flank = Counter(rig_votes).most_common(1)[0][0]
            chosen_tmpl = _select_template_by_length_and_rigidity(
                template_records, pred_len_from_flank, pred_rig_from_flank
            )
            arm_examples["flanking_context_template_retrieval"].append(
                _score_prediction(
                    str(chosen_tmpl["linker_sequence"]),
                    str(chosen_tmpl["biophysical_properties"]["rigidity_predicted_label"]),
                    test_row,
                )
            )

            # 4. Flanking-context Nearest-Neighbour from Natural Authorities
            arm_examples["flanking_context_nn_natural_authority"].append(
                _score_prediction(
                    str(best_nat["linker_sequence"]),
                    str(best_nat["biophysical_properties"]["rigidity_predicted_label"]),
                    test_row,
                )
            )

            # 5. Flanking-context Nearest-Neighbour from Natural + Mined Authorities
            short_combined = _shortlist_by_flank_composition(query_flank, combined_pool, limit=250)
            ranked_combined = sorted(
                short_combined,
                key=lambda cand: (
                    -flanking_pair_similarity(query_flank, cand.get("flanking_context") or {}),
                    cand["record_id"],
                ),
            )
            best_comb = ranked_combined[0]
            arm_examples["flanking_context_nn_natural_plus_mined"].append(
                _score_prediction(
                    str(best_comb["linker_sequence"]),
                    str(best_comb["biophysical_properties"]["rigidity_predicted_label"]),
                    test_row,
                )
            )

            # 6 & 7. 3D Geometry Controls (when 3D attachment coordinates are present)
            if has_3d_geometry and test_row.get("geometry"):
                ca_dist = float(test_row["geometry"]["inter_port_ca_distance_A"])
                min_req = minimum_residues(ca_dist)
                # Physical target length from X1 (2.5 A/res span match, bounded by comfortable_residues)
                comf_req = comfortable_residues(ca_dist) or max(5, round(ca_dist / 2.5))
                geom_target_len = max(min_req, min(30, round(0.5 * (comf_req + ca_dist / 2.5))))

                # 6. Geometry-only control (closest training linker to geometric target length, no flanking seq)
                best_geom_only = min(
                    train_rows,
                    key=lambda cand: (
                        int(cand["linker_length"]) < min_req,
                        abs(int(cand["linker_length"]) - geom_target_len),
                        cand["record_id"],
                    ),
                )
                arm_examples["geometry_only_control"].append(
                    _score_prediction(
                        str(best_geom_only["linker_sequence"]),
                        str(best_geom_only["biophysical_properties"]["rigidity_predicted_label"]),
                        test_row,
                    )
                )

                # 7. Geometry + Flanking Hybrid Retrieval (bridgeable window + flanking similarity)
                bridgeable_pool = [
                    cand
                    for cand in combined_pool
                    if (int(cand["linker_length"]) + 1) * CA_VIRTUAL_BOND_ANGSTROM >= ca_dist
                ]
                search_pool = bridgeable_pool if bridgeable_pool else combined_pool
                short_search = _shortlist_by_flank_composition(query_flank, search_pool, limit=250)
                best_hybrid = max(
                    short_search,
                    key=lambda cand: (
                        flanking_pair_similarity(query_flank, cand.get("flanking_context") or {})
                        - 0.03 * abs(int(cand["linker_length"]) - geom_target_len),
                        cand["record_id"],
                    ),
                )
                arm_examples["geometry_plus_flanking_hybrid_retrieval"].append(
                    _score_prediction(
                        str(best_hybrid["linker_sequence"]),
                        str(best_hybrid["biophysical_properties"]["rigidity_predicted_label"]),
                        test_row,
                    )
                )

            # 8. Oracle Length-Matched Template Control (`CURRENT_DIRECTION_2026-10-18.md` item 4)
            oracle_tmpl = _select_template_by_length_and_rigidity(
                template_records, float(true_len), pred_rig_from_flank
            )
            arm_examples["length_matched_template_control"].append(
                _score_prediction(
                    str(oracle_tmpl["linker_sequence"]),
                    str(oracle_tmpl["biophysical_properties"]["rigidity_predicted_label"]),
                    test_row,
                )
            )

            # 9. Oracle Length-Matched Frequent Dipeptide Control (`CHECKLIST.md` item X2)
            freq_seq = _most_frequent_dipeptide_repeat(train_rows, true_len)
            arm_examples["length_matched_frequent_natural_control"].append(
                _score_prediction(freq_seq, "flexible", test_row)
            )

    summary_arms = {
        arm_name: _aggregate_arm(examples, seed=seed)
        for arm_name, examples in arm_examples.items()
    }
    return {
        "n_folds": n_folds,
        "seed": seed,
        "total_evaluated_natural_linkers": len(natural_records),
        "distinct_homology_clusters": len({r["homology_cluster_id"] for r in natural_records}),
        "mined_candidate_pool_size": len(mined_records),
        "template_pool_size": len(template_records),
        "arms": summary_arms,
    }
