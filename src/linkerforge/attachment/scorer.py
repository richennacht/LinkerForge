"""Position-controlled attachment-site scorer and closed-form geometric site-pair re-ranker.

Implements Component 4 (`ARCHITECTURE_HANDOFF.md` section 6):
  1. `TwoTowerAttachmentScorer` -- a regularized logistic model over ten
     *position-free* local sequence features (loop/turn propensity, hydrophilicity,
     backbone-flexibility scales, Gly/Pro/Ser/Asn/Asp, bulky-hydrophobic and charged
     fractions, loop contrast).  The two positional features of the 12-slot feature
     vector (`term_prox`, `rel_interior`) are retained in the vector for API
     compatibility but their weights are fixed to zero, because training or
     evaluating with them lets a model recover the label from residue position alone.
     `rigidity` / `cleavability` tick-boxes never enter the encoder.
  2. `score_residue_with_structure` -- a hand-set (NOT learned) heuristic that mixes
     the sequence score with Shrake-Rupley relative accessibility and a coarse C-alpha
     curvature class when coordinates are supplied.  Its weights (0.52 / 0.38 / +-0.2)
     are design constants and have no held-out validation.
  3. `rank_attachment_site_pairs` -- combines the per-site scores with the
     deterministic junction geometry from `linkerforge.geometry`.
  4. `evaluate_attachment_scorer_benchmark` -- a position-matched control: for every
     linker-flanking sequence window, residues adjacent to the linker (positives) are
     compared with residues at the same offset from the opposite truncation edge of the
     same window (negatives), so residue position carries no information.  Evaluation is
     cluster-disjoint (`homology_cluster_id`) with cluster-bootstrap confidence
     intervals.  `diagnose_original_task_leakage` reproduces the earlier benchmark
     design to document why its near-perfect scores were uninformative.

Scope: positives mean "adjacent to a natural inter-domain linker".  This is NOT an
experimentally measured insertion/attachment tolerance label.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass
from math import acos, degrees, dist, exp, pi, sqrt
import random
from typing import Any, Sequence

from linkerforge.geometry.junction import evaluate_junction
from linkerforge.geometry.sasa import relative_accessibility, sphere_points
from linkerforge.geometry.structure import Chain, Residue
from linkerforge.geometry.tether import implied_chain_topology

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# Chou-Fasman / Levitt normalized loop/turn & flexibility propensities
LOOP_TURN_PROPENSITY: dict[str, float] = {
    "G": 1.56,
    "P": 1.52,
    "S": 1.43,
    "D": 1.46,
    "N": 1.56,
    "T": 0.96,
    "K": 1.01,
    "E": 0.74,
    "Q": 0.98,
    "R": 0.95,
    "A": 0.66,
    "H": 0.95,
    "Y": 1.14,
    "C": 1.19,
    "W": 0.96,
    "F": 0.60,
    "M": 0.60,
    "L": 0.59,
    "V": 0.50,
    "I": 0.47,
}

# Kyte-Doolittle hydropathy (normalized so hydrophilic surface residues are positive)
HYDROPHILICITY: dict[str, float] = {
    "R": 1.00,
    "K": 0.87,
    "D": 0.78,
    "E": 0.78,
    "N": 0.78,
    "Q": 0.78,
    "H": 0.71,
    "P": 0.36,
    "Y": 0.29,
    "W": 0.20,
    "S": 0.18,
    "T": 0.16,
    "G": 0.09,
    "A": -0.40,
    "M": -0.42,
    "C": -0.56,
    "F": -0.62,
    "L": -0.84,
    "V": -0.93,
    "I": -1.00,
}

# Vihinen et al. B-factor backbone flexibility scale
BACKBONE_FLEXIBILITY: dict[str, float] = {
    "G": 1.04,
    "S": 1.05,
    "D": 1.03,
    "K": 1.09,
    "E": 1.10,
    "P": 1.05,
    "N": 1.04,
    "Q": 1.04,
    "R": 1.03,
    "T": 1.00,
    "A": 0.98,
    "H": 0.95,
    "Y": 0.93,
    "M": 0.95,
    "V": 0.94,
    "L": 0.94,
    "I": 0.93,
    "F": 0.91,
    "C": 0.92,
    "W": 0.90,
}


@dataclass(frozen=True)
class AttachmentSiteCandidate:
    residue_index_1based: int
    author_residue_number: int
    amino_acid: str
    role: str  # "c_terminal_port" | "n_terminal_port" | "internal_loop_site"
    two_tower_tolerance_score: float
    sequence_flexibility_score: float
    relative_sasa: float | None
    local_ca_curvature_degrees: float | None
    secondary_structure_proxy: str
    admissible_surface_site: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class AttachmentSitePairRanking:
    rank: int
    domain_a_site: AttachmentSiteCandidate
    domain_b_site: AttachmentSiteCandidate
    topology: str
    inter_port_ca_distance_A: float
    required_path_length_A: float
    path_detour_ratio: float
    minimum_physical_residues: int
    minimum_design_residues: int
    maximum_design_residues: int | None
    geometric_verdict: str
    two_tower_pair_tolerance: float
    closed_form_geometric_score: float
    combined_site_pair_score: float
    flags: list[str]

    def to_dict(self) -> dict[str, Any]:
        out = asdict(self)
        out["domain_a_site"] = self.domain_a_site.to_dict()
        out["domain_b_site"] = self.domain_b_site.to_dict()
        return out


def _sigmoid(x: float) -> float:
    x_clamped = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + exp(-x_clamped))


def _local_ca_curvature_and_ss(residues: Sequence[Residue], idx: int) -> tuple[float, str]:
    """Estimate local C-alpha bend angle and secondary structure proxy (`loop`, `helix`, `sheet`)."""
    n = len(residues)
    if n < 3 or idx <= 0 or idx >= n - 1:
        return 65.0, "terminal_loop"
    ca_prev = residues[idx - 1].ca
    ca_curr = residues[idx].ca
    ca_next = residues[idx + 1].ca
    if ca_prev is None or ca_curr is None or ca_next is None:
        return 60.0, "loop"

    v1 = (ca_curr[0] - ca_prev[0], ca_curr[1] - ca_prev[1], ca_curr[2] - ca_prev[2])
    v2 = (ca_next[0] - ca_curr[0], ca_next[1] - ca_curr[1], ca_next[2] - ca_curr[2])
    n1 = sqrt(sum(c * c for c in v1))
    n2 = sqrt(sum(c * c for c in v2))
    if n1 <= 1e-6 or n2 <= 1e-6:
        return 60.0, "loop"
    cos_theta = max(-1.0, min(1.0, sum(a * b for a, b in zip(v1, v2)) / (n1 * n2)))
    bend_deg = degrees(acos(cos_theta))

    # Check i-1 to i+2 C-alpha distance to distinguish alpha-helix (~5.1-5.8 A) vs extended sheet (~9.5-11.0 A)
    if 1 <= idx <= n - 3:
        c0 = residues[idx - 1].ca
        c3 = residues[idx + 2].ca
        if c0 is not None and c3 is not None:
            d14 = dist(c0, c3)
            if 4.8 <= d14 <= 6.0 and 75.0 <= bend_deg <= 100.0:
                return bend_deg, "helix_core"
            if d14 >= 9.2 and bend_deg <= 48.0:
                return bend_deg, "beta_strand"
    return bend_deg, "loop"


def extract_residue_tolerance_features(
    sequence: str,
    pos_0based: int,
    *,
    attachment_mode: str = "any",
) -> list[float]:
    """Extract 12 single-domain sequence features around `pos_0based`.

    Tick-boxes (`rigidity`, `cleavability`) are intentionally excluded (`ARCHITECTURE_HANDOFF.md` §6).
    """
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        return [0.0] * 12

    win_3 = seq[max(0, pos_0based - 1) : min(n, pos_0based + 2)]
    win_7 = seq[max(0, pos_0based - 3) : min(n, pos_0based + 4)]
    win_11 = seq[max(0, pos_0based - 5) : min(n, pos_0based + 6)]

    def _mean_prop(sub: str, table: dict[str, float]) -> float:
        vals = [table.get(ch, 0.0) for ch in sub]
        return sum(vals) / len(vals) if vals else 0.0

    loop_3 = _mean_prop(win_3, LOOP_TURN_PROPENSITY)
    loop_7 = _mean_prop(win_7, LOOP_TURN_PROPENSITY)
    hydro_3 = _mean_prop(win_3, HYDROPHILICITY)
    hydro_7 = _mean_prop(win_7, HYDROPHILICITY)
    flex_3 = _mean_prop(win_3, BACKBONE_FLEXIBILITY)
    flex_11 = _mean_prop(win_11, BACKBONE_FLEXIBILITY)

    gly_pro_ser_7 = sum(1 for ch in win_7 if ch in "GPSND") / len(win_7)
    bulky_hydro_7 = sum(1 for ch in win_7 if ch in "FILVWY") / len(win_7)
    charged_7 = sum(1 for ch in win_7 if ch in "DEKR") / len(win_7)

    dist_n = pos_0based / float(max(1, n - 1))
    dist_c = (n - 1 - pos_0based) / float(max(1, n - 1))
    if attachment_mode == "c_terminal":
        term_prox = exp(-float(n - 1 - pos_0based) / 6.0)
    elif attachment_mode == "n_terminal":
        term_prox = exp(-float(pos_0based) / 6.0)
    else:
        term_prox = max(exp(-float(pos_0based) / 6.0), exp(-float(n - 1 - pos_0based) / 6.0))

    loop_contrast = loop_3 - _mean_prop(win_11, LOOP_TURN_PROPENSITY)
    rel_interior = 4.0 * dist_n * dist_c

    return [
        loop_3,
        loop_7,
        hydro_3,
        hydro_7,
        flex_3,
        flex_11,
        gly_pro_ser_7,
        bulky_hydro_7,
        charged_7,
        term_prox,
        loop_contrast,
        rel_interior,
    ]


@dataclass(frozen=True)
class TwoTowerAttachmentScorer:
    """Learned single-domain attachment-site tolerance scorer."""

    weights: tuple[float, ...]
    bias: float
    feature_means: tuple[float, ...]
    feature_stds: tuple[float, ...]

    def score_sequence_position(
        self,
        sequence: str,
        pos_0based: int,
        *,
        attachment_mode: str = "any",
    ) -> float:
        feats = extract_residue_tolerance_features(
            sequence, pos_0based, attachment_mode=attachment_mode
        )
        z = self.bias
        for w, x, mu, sd in zip(self.weights, feats, self.feature_means, self.feature_stds):
            z += w * ((x - mu) / max(1e-6, sd))
        return _sigmoid(z)

    def score_residue_with_structure(
        self,
        sequence: str,
        pos_0based: int,
        *,
        relative_sasa: float | None = None,
        ss_proxy: str = "loop",
        attachment_mode: str = "any",
    ) -> float:
        seq_prob = self.score_sequence_position(
            sequence, pos_0based, attachment_mode=attachment_mode
        )
        if relative_sasa is None:
            return seq_prob
        # Combine sequence tower with 3D surface accessibility & secondary-structure preservation
        sasa_term = min(1.0, max(0.0, relative_sasa / 0.45))
        if ss_proxy in ("loop", "terminal_loop"):
            ss_bonus = 0.15
        elif ss_proxy == "helix_core":
            ss_bonus = -0.22
        else:
            ss_bonus = -0.18
        combined = 0.52 * seq_prob + 0.38 * sasa_term + ss_bonus
        return max(0.01, min(0.99, combined))


POSITION_FREE_FEATURE_INDICES: tuple[int, ...] = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10)
POSITIONAL_FEATURE_INDICES: tuple[int, ...] = (9, 11)
MIN_WINDOW_LENGTH = 25
EDGE_OFFSETS: tuple[int, ...] = (0, 1, 2)

_FALLBACK_WEIGHTS = (0.35, 0.30, 0.25, 0.20, 0.35, 0.20, 0.40, -0.45, 0.15, 0.0, 0.25, 0.0)


def matched_offset_positions(sequence_length: int, *, linker_side: str) -> tuple[list[int], list[int]]:
    """Return (positive, negative) 0-based positions with identical edge offsets.

    `linker_side="c_terminal"`: the linker follows the end of this window (left flank), so
    positives are the last residues and negatives the first residues.  `"n_terminal"` is the
    mirror image (right flank).  Both sets sit at offsets 0, 1, 2 from a truncation edge of
    the same window, so position-dependent window clipping affects them identically.
    """
    n = sequence_length
    near_edge = [n - 1 - k for k in EDGE_OFFSETS]
    far_edge = list(EDGE_OFFSETS)
    if linker_side == "c_terminal":
        return near_edge, far_edge
    if linker_side == "n_terminal":
        return far_edge, near_edge
    raise ValueError(f"linker_side must be 'c_terminal' or 'n_terminal', got {linker_side!r}")


def _flank_windows(row: dict[str, Any]) -> list[tuple[str, str]]:
    flank = row.get("flanking_context") or {}
    out: list[tuple[str, str]] = []
    for key, side in (("left_flank_sequence", "c_terminal"), ("right_flank_sequence", "n_terminal")):
        seq = str(flank.get(key) or "").strip().upper()
        if len(seq) >= MIN_WINDOW_LENGTH:
            out.append((side, seq))
    return out


def fit_two_tower_attachment_scorer(
    training_records: Sequence[dict[str, Any]],
    *,
    epochs: int = 120,
    lr: float = 0.08,
    l2_reg: float = 0.01,
) -> TwoTowerAttachmentScorer:
    """Fit the scorer on position-matched examples using position-free features only.

    Positives: the three residues adjacent to a natural linker.  Negatives: the three
    residues at the same offset from the opposite edge of the same flank window.  The
    positional feature slots (`term_prox`, `rel_interior`) receive weight 0.
    """
    x_rows: list[list[float]] = []
    y_rows: list[int] = []
    for rec in training_records:
        for side, seq in _flank_windows(rec):
            pos_idx, neg_idx = matched_offset_positions(len(seq), linker_side=side)
            for idx, y in [(i, 1) for i in pos_idx] + [(i, 0) for i in neg_idx]:
                feats = extract_residue_tolerance_features(seq, idx, attachment_mode="any")
                x_rows.append([feats[d] for d in POSITION_FREE_FEATURE_INDICES])
                y_rows.append(y)

    if not x_rows:
        return TwoTowerAttachmentScorer(
            weights=_FALLBACK_WEIGHTS,
            bias=0.0,
            feature_means=tuple(0.0 for _ in range(12)),
            feature_stds=tuple(1.0 for _ in range(12)),
        )

    dim = len(x_rows[0])
    n_samples = float(len(x_rows))
    means = [sum(row[d] for row in x_rows) / n_samples for d in range(dim)]
    stds = [
        max(1e-4, sqrt(sum((row[d] - means[d]) ** 2 for row in x_rows) / n_samples))
        for d in range(dim)
    ]
    x_norm = [[(row[d] - means[d]) / stds[d] for d in range(dim)] for row in x_rows]

    w = [0.0] * dim
    b = 0.0
    pos_count = max(1, sum(y_rows))
    neg_count = max(1, len(y_rows) - pos_count)
    pos_weight = len(y_rows) / (2.0 * pos_count)
    neg_weight = len(y_rows) / (2.0 * neg_count)
    for _ in range(epochs):
        grad_w = [0.0] * dim
        grad_b = 0.0
        for feats, y in zip(x_norm, y_rows):
            pred = _sigmoid(b + sum(wj * xj for wj, xj in zip(w, feats)))
            err = (pred - float(y)) * (pos_weight if y == 1 else neg_weight)
            grad_b += err
            for d in range(dim):
                grad_w[d] += err * feats[d]
        b -= lr * (grad_b / n_samples)
        for d in range(dim):
            w[d] -= lr * ((grad_w[d] / n_samples) + l2_reg * w[d])

    full_w = [0.0] * 12
    full_mu = [0.0] * 12
    full_sd = [1.0] * 12
    for local, slot in enumerate(POSITION_FREE_FEATURE_INDICES):
        full_w[slot] = round(w[local], 6)
        full_mu[slot] = round(means[local], 6)
        full_sd[slot] = round(stds[local], 6)
    return TwoTowerAttachmentScorer(
        weights=tuple(full_w),
        bias=round(b, 6),
        feature_means=tuple(full_mu),
        feature_stds=tuple(full_sd),
    )


def score_domain_attachment_sites(
    domain_sequence: str,
    scorer: TwoTowerAttachmentScorer,
    *,
    chain: Chain | None = None,
    attachment_mode: str = "any",
    top_k: int | None = None,
) -> list[AttachmentSiteCandidate]:
    """Score every residue in a single domain for linker attachment / insertion tolerance."""
    seq = domain_sequence.strip().upper()
    sasa_map: dict[int, float] = {}
    residues_list: list[Residue] = []
    if chain is not None and chain.residues:
        residues_list = list(chain.residues)
        ctx_atoms = tuple(atom for res in residues_list for atom in res.atoms)
        pts = sphere_points(96)
        for res in residues_list:
            acc = relative_accessibility(res, ctx_atoms, points=pts)
            rel_val = acc.get("relative_solvent_accessibility")
            sasa_map[res.number] = float(rel_val) if rel_val is not None else 0.35
        if not seq:
            seq = "".join(
                ALPHABET[0] if len(r.name) != 3 else "G" for r in residues_list
            )

    n = len(residues_list) if residues_list else len(seq)
    candidates: list[AttachmentSiteCandidate] = []
    for idx in range(n):
        aa = seq[idx] if idx < len(seq) else "G"
        res_num = residues_list[idx].number if residues_list else (idx + 1)
        rel_sasa = sasa_map.get(res_num) if residues_list else None
        curv_deg, ss_proxy = (
            _local_ca_curvature_and_ss(residues_list, idx)
            if residues_list
            else (60.0, "loop")
        )
        seq_score = scorer.score_sequence_position(
            seq, idx, attachment_mode=attachment_mode
        )
        tol_score = scorer.score_residue_with_structure(
            seq,
            idx,
            relative_sasa=rel_sasa,
            ss_proxy=ss_proxy,
            attachment_mode=attachment_mode,
        )
        if idx >= n - 5:
            role = "c_terminal_port"
        elif idx < 5:
            role = "n_terminal_port"
        else:
            role = "internal_loop_site"
        admissible = True if rel_sasa is None else (rel_sasa >= 0.15 and ss_proxy != "helix_core")
        candidates.append(
            AttachmentSiteCandidate(
                residue_index_1based=idx + 1,
                author_residue_number=res_num,
                amino_acid=aa,
                role=role,
                two_tower_tolerance_score=round(tol_score, 6),
                sequence_flexibility_score=round(seq_score, 6),
                relative_sasa=None if rel_sasa is None else round(rel_sasa, 4),
                local_ca_curvature_degrees=round(curv_deg, 2),
                secondary_structure_proxy=ss_proxy,
                admissible_surface_site=admissible,
            )
        )

    candidates.sort(
        key=lambda c: (-c.two_tower_tolerance_score, c.residue_index_1based)
    )
    return candidates[:top_k] if top_k is not None else candidates


def rank_attachment_site_pairs(
    domain_a: Chain,
    domain_a_sequence: str,
    domain_b: Chain,
    domain_b_sequence: str,
    scorer: TwoTowerAttachmentScorer,
    *,
    linker_type: str = "flexible",
    requested_residues: int | None = None,
    topology: str = "end_to_end_tandem",
    candidate_sites_per_domain: int = 8,
    top_k: int = 10,
) -> list[AttachmentSitePairRanking]:
    """Rank attachment site pairs `(site_a, site_b)` via Two-Tower + Closed-Form Geometry."""
    mode_a = "c_terminal" if topology == "end_to_end_tandem" else "any"
    mode_b = "n_terminal" if topology == "end_to_end_tandem" else "any"

    sites_a = score_domain_attachment_sites(
        domain_a_sequence,
        scorer,
        chain=domain_a,
        attachment_mode=mode_a,
        top_k=candidate_sites_per_domain,
    )
    sites_b = score_domain_attachment_sites(
        domain_b_sequence,
        scorer,
        chain=domain_b,
        attachment_mode=mode_b,
        top_k=candidate_sites_per_domain,
    )

    topo_info = (
        implied_chain_topology(1, peptide_tethers=1)
        if topology == "end_to_end_tandem"
        else implied_chain_topology(2, peptide_tethers=2)
    )

    pair_items: list[AttachmentSitePairRanking] = []
    for sa in sites_a:
        for sb in sites_b:
            junc = evaluate_junction(
                domain_a,
                sa.author_residue_number,
                domain_b,
                sb.author_residue_number,
                linker_type=linker_type,
                requested_residues=requested_residues,
            )
            span_win = junc["required_linker_length"]
            span_info = junc["span"]
            orient_a = junc["orientation"]["domain_a"]
            orient_b = junc["orientation"]["domain_b"]

            ca_dist = float(span_info["straight_ca_distance_angstrom"])
            req_path = float(span_info["path_length_used_angstrom"])
            detour_ratio = req_path / max(1.0, ca_dist)

            min_phys = int(span_win["minimum_residues"])
            min_des = int(span_win.get("minimum_residues_design", min_phys))
            max_des = span_win.get("maximum_residues")

            if requested_residues is not None:
                verdict = str(junc["requested_length_verdict"]["verdict"])
            else:
                verdict = "bridgeable" if min_des <= 45 else "long_span"

            # Closed-form geometric score in [0, 1]
            faces_a = 1.0 if orient_a.get("anchor_faces_partner_domain") else 0.45
            faces_b = 1.0 if orient_b.get("anchor_faces_partner_domain") else 0.45
            detour_score = max(0.0, min(1.0, 1.0 - 0.8 * max(0.0, detour_ratio - 1.0)))
            span_score = max(0.1, min(1.0, exp(-max(0.0, ca_dist - 25.0) / 35.0)))
            verdict_bonus = 1.0 if verdict == "bridgeable" else 0.35

            geom_score = verdict_bonus * (
                0.35 * detour_score + 0.35 * (0.5 * (faces_a + faces_b)) + 0.30 * span_score
            )
            pair_tol = 0.5 * (sa.two_tower_tolerance_score + sb.two_tower_tolerance_score)
            combined_score = 0.55 * pair_tol + 0.45 * geom_score

            flags: list[str] = []
            if sa.admissible_surface_site and sb.admissible_surface_site:
                flags.append("both_anchors_surface_accessible")
            else:
                flags.append("buried_or_structured_anchor_warning")
            if detour_ratio > 1.08:
                flags.append("steric_domain_detour_required")
            if verdict == "bridgeable":
                flags.append("span_window_bridgeable")
            else:
                flags.append(f"span_verdict_{verdict}")
            if not bool(topo_info.get("single_chain_realizable", True)):
                flags.append("multi_chain_or_circular_topology_required")

            pair_items.append(
                AttachmentSitePairRanking(
                    rank=0,
                    domain_a_site=sa,
                    domain_b_site=sb,
                    topology=topology,
                    inter_port_ca_distance_A=round(ca_dist, 4),
                    required_path_length_A=round(req_path, 4),
                    path_detour_ratio=round(detour_ratio, 4),
                    minimum_physical_residues=min_phys,
                    minimum_design_residues=min_des,
                    maximum_design_residues=None if max_des is None else int(max_des),
                    geometric_verdict=verdict,
                    two_tower_pair_tolerance=round(pair_tol, 6),
                    closed_form_geometric_score=round(geom_score, 6),
                    combined_site_pair_score=round(combined_score, 6),
                    flags=flags,
                )
            )

    pair_items.sort(
        key=lambda p: (
            -p.combined_site_pair_score,
            p.inter_port_ca_distance_A,
            p.domain_a_site.author_residue_number,
            p.domain_b_site.author_residue_number,
        )
    )
    ranked: list[AttachmentSitePairRanking] = []
    for idx, item in enumerate(pair_items[:top_k], start=1):
        ranked.append(
            AttachmentSitePairRanking(
                rank=idx,
                domain_a_site=item.domain_a_site,
                domain_b_site=item.domain_b_site,
                topology=item.topology,
                inter_port_ca_distance_A=item.inter_port_ca_distance_A,
                required_path_length_A=item.required_path_length_A,
                path_detour_ratio=item.path_detour_ratio,
                minimum_physical_residues=item.minimum_physical_residues,
                minimum_design_residues=item.minimum_design_residues,
                maximum_design_residues=item.maximum_design_residues,
                geometric_verdict=item.geometric_verdict,
                two_tower_pair_tolerance=item.two_tower_pair_tolerance,
                closed_form_geometric_score=item.closed_form_geometric_score,
                combined_site_pair_score=item.combined_site_pair_score,
                flags=item.flags,
            )
        )
    return ranked


def _compute_auroc_and_aps(labels: Sequence[int], scores: Sequence[float]) -> tuple[float, float]:
    paired = sorted(zip(scores, labels), key=lambda x: -x[0])
    total_pos = sum(labels)
    total_neg = len(labels) - total_pos
    if total_pos == 0 or total_neg == 0:
        return 0.5, 0.0

    tp = 0
    fp = 0
    prev_tpr = 0.0
    prev_fpr = 0.0
    auroc = 0.0
    aps = 0.0

    for sc, y in paired:
        if y == 1:
            tp += 1
            prec = tp / float(tp + fp)
            aps += prec / float(total_pos)
        else:
            fp += 1
        tpr = tp / float(total_pos)
        fpr = fp / float(total_neg)
        auroc += 0.5 * (tpr + prev_tpr) * (fpr - prev_fpr)
        prev_tpr = tpr
        prev_fpr = fpr

    return auroc, aps


def _precision_at_coverage(labels: Sequence[int], scores: Sequence[float], coverage: float) -> float:
    if not labels:
        return 0.0
    k = max(1, round(len(labels) * coverage))
    paired = sorted(zip(scores, labels), key=lambda x: -x[0])[:k]
    return sum(y for _, y in paired) / float(k)


def _auroc_ranks(labels: Sequence[int], scores: Sequence[float]) -> float:
    """Tie-aware rank-sum AUROC (ties count 0.5)."""
    n_pos = sum(1 for y in labels if y == 1)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        avg = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    rank_sum_pos = sum(r for r, y in zip(ranks, labels) if y == 1)
    return (rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def _cluster_bootstrap_auroc(
    cluster_ids: Sequence[str],
    labels: Sequence[int],
    scores: Sequence[float],
    *,
    reps: int = 300,
    seed: int = 0,
) -> list[float]:
    by_cluster: dict[str, list[int]] = defaultdict(list)
    for idx, cid in enumerate(cluster_ids):
        by_cluster[cid].append(idx)
    ids = sorted(by_cluster)
    rng = random.Random(seed)
    stats: list[float] = []
    for _ in range(reps):
        picked = [rng.choice(ids) for _ in ids]
        idxs = [i for cid in picked for i in by_cluster[cid]]
        stats.append(_auroc_ranks([labels[i] for i in idxs], [scores[i] for i in idxs]))
    stats.sort()
    return [round(stats[int(0.025 * reps)], 4), round(stats[int(0.975 * reps) - 1], 4)]


def _flexibility_propensity(residue: str) -> float:
    return 0.55 * BACKBONE_FLEXIBILITY.get(residue, 0.98) + 0.45 * (
        LOOP_TURN_PROPENSITY.get(residue, 0.8) / 1.56
    )


def evaluate_attachment_scorer_benchmark(
    natural_records: Sequence[dict[str, Any]],
    *,
    n_folds: int = 5,
    seed: int = 17,
    bootstrap_reps: int = 300,
) -> dict[str, Any]:
    """Position-matched, cluster-disjoint evaluation of the sequence-context scorer.

    For every flank window the six candidate residues are the three residues adjacent to the
    linker (label 1) and the three residues at the same offsets from the opposite edge (label 0).
    Residue position therefore carries no information; the chance level is AUROC 0.5.
    No label-derived quantity (accessibility, secondary structure) is passed to any arm.
    """
    from linkerforge.retrieval.baseline import assign_cluster_disjoint_folds

    folds = assign_cluster_disjoint_folds(natural_records, n_folds=n_folds, seed=seed)
    rng = random.Random(seed)

    arm_names = (
        "sequence_context_logistic",
        "single_residue_flexibility_propensity",
        "random_control",
    )
    labels: list[int] = []
    clusters: list[str] = []
    scores: dict[str, list[float]] = {a: [] for a in arm_names}
    window_acc: dict[str, list[float]] = {a: [] for a in arm_names}
    n_windows = 0

    for train_rows, test_rows in folds:
        scorer = fit_two_tower_attachment_scorer(train_rows)
        for row in test_rows:
            cid = str(row.get("homology_cluster_id") or row.get("record_id"))
            for side, seq in _flank_windows(row):
                n_windows += 1
                pos_idx, neg_idx = matched_offset_positions(len(seq), linker_side=side)
                idxs = pos_idx + neg_idx
                ys = [1] * len(pos_idx) + [0] * len(neg_idx)
                per_arm: dict[str, list[float]] = {a: [] for a in arm_names}
                for i in idxs:
                    per_arm["sequence_context_logistic"].append(
                        scorer.score_sequence_position(seq, i, attachment_mode="any")
                    )
                    per_arm["single_residue_flexibility_propensity"].append(
                        _flexibility_propensity(seq[i])
                    )
                    per_arm["random_control"].append(rng.random())
                labels.extend(ys)
                clusters.extend([cid] * len(ys))
                for a in arm_names:
                    scores[a].extend(per_arm[a])
                    window_acc[a].append(_auroc_ranks(ys, per_arm[a]))

    arms: dict[str, Any] = {}
    for a in arm_names:
        arms[a] = {
            "auroc_pooled": round(_auroc_ranks(labels, scores[a]), 4),
            "auroc_cluster_bootstrap_95ci": _cluster_bootstrap_auroc(
                clusters, labels, scores[a], reps=bootstrap_reps, seed=seed
            ),
            "within_window_auroc_mean": round(sum(window_acc[a]) / max(1, len(window_acc[a])), 4),
        }
    return {
        "task": "position_matched_adjacent_to_linker_vs_opposite_edge_same_window",
        "chance_auroc": 0.5,
        "n_folds": n_folds,
        "seed": seed,
        "n_records": len(natural_records),
        "n_windows": n_windows,
        "n_candidate_residues": len(labels),
        "n_homology_clusters": len(set(clusters)),
        "label_meaning": "adjacent to a natural inter-domain linker (not an experimental tolerance label)",
        "arms": arms,
    }


def diagnose_original_task_leakage(natural_records: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Reproduce the earlier benchmark design and show position alone solves it.

    Original design: positives = last three (left) / first three (right) residues of the window,
    negatives = interior residues.  Any score that is a function of residue position alone
    separates these perfectly, so the earlier near-1.0 AUROC (and its "no terminal prior"
    ablation, which left the second positional feature in place) carried no information.
    The earlier structure arm additionally set accessibility and secondary-structure inputs
    from the label itself; that arm is not reproduced here.
    """
    labels: list[int] = []
    pos_only: list[float] = []
    rel_interior_only: list[float] = []
    for row in natural_records:
        for side, seq in _flank_windows(row):
            n = len(seq)
            if side == "c_terminal":
                pos_idx = list(range(n - 3, n))
                neg_idx = list(range(4, n - 6, 3))
            else:
                pos_idx = list(range(0, 3))
                neg_idx = list(range(6, n - 4, 3))
            for idx, y in [(i, 1) for i in pos_idx] + [(i, 0) for i in neg_idx]:
                labels.append(y)
                dist_n = idx / float(max(1, n - 1))
                dist_c = (n - 1 - idx) / float(max(1, n - 1))
                if side == "c_terminal":
                    pos_only.append(exp(-float(n - 1 - idx) / 6.0))
                else:
                    pos_only.append(exp(-float(idx) / 6.0))
                rel_interior_only.append(1.0 - 4.0 * dist_n * dist_c)
    return {
        "original_design": "last/first 3 residues vs interior residues of the same window",
        "n_candidate_residues": len(labels),
        "n_positive_residues": sum(labels),
        "auroc_position_only_term_prox": round(_auroc_ranks(labels, pos_only), 4),
        "auroc_position_only_rel_interior": round(_auroc_ranks(labels, rel_interior_only), 4),
        "interpretation": (
            "A score that depends only on distance to the window edge separates the original "
            "labels almost perfectly; the original benchmark therefore did not test sequence "
            "or structure."
        ),
    }
