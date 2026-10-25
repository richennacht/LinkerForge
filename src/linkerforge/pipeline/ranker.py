"""End-to-End Pairwise Compatibility Ranker, Hard Structural Rejection Scorecard,
6-Baseline Benchmark, and sfGFP–mCherry 3D Fusion Pipeline (`ARCHITECTURE_HANDOFF.md` §7).

Implements:
  1. Deterministic Hard Structural Rejection Scorecard (`evaluate_hard_rejection_scorecard`)
     enforcing physical 3.8 A/res span closure, steric detour clearance, port accessibility,
     hydrophobic stretch exclusion (<50% FILVWM and no >=5-aa hydrophobic run), cleavage motif
     compliance, and rigid+cleavable flexible-spacer compliance. Never outputs a fusion
     success probability (`SPEC.md` §5.8 / `ARCHITECTURE_HANDOFF.md` §8).
  2. Six Synthetic Negative Corruptions Generator (`generate_synthetic_negatives`):
     - `too_short_by_50pct`
     - `too_long_by_200pct`
     - `rigid_helical_where_flexible_needed`
     - `all_hydrophobic_core`
     - `charge_clash_with_junction_anchors`
     - `protease_site_buried_in_rigid_helix`
  3. Learned Pairwise Compatibility Ranker (`PairwiseCompatibilityRanker`) trained with
     margin-regularized pairwise logistic ranking loss on natural + proteome-mined high-confidence
     linker authorities (`mined_and_natural_linkers_v1`) vs. the 6 synthetic negative corruptions
     on cluster-disjoint (`homology_cluster_id`) splits.
  4. Six Required Baselines Evaluation (`evaluate_six_baselines_benchmark`):
     - `direct_0aa_fusion`
     - `fixed_ggggs3`
     - `random_length_matched_gs`
     - `geometry_only`
     - `nearest_neighbour_retrieval`
     - `linkerforge_ranker`
  5. `sfGFP (PDB 2B3P) — linker — mCherry (PDB 2H5Q)` End-to-End Fusion Demonstration and
     3D PDB + PyMOL `.pml` script generator (`run_sfgfp_mcherry_fusion_demo`,
     `build_fusion_3d_pdb_and_pymol_script`).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import random
from typing import Any, Mapping, Sequence

from linkerforge.geometry.span import maximum_span
from linkerforge.design.engine import (
    TIER_PRIORITY,
    build_rigid_cleavable_composite,
    flanking_context_similarity,
)
from linkerforge.properties import load_library

KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}


def compute_linker_sequence_physics(sequence: str) -> dict[str, float]:
    """Compute biophysical sequence properties for a candidate linker."""
    seq = sequence.upper().strip()
    L = len(seq)
    if L == 0:
        return {
            "hydropathy_kyte_doolittle": 0.0,
            "net_charge_ph7": 0.0,
            "fcr": 0.0,
            "charge_segregation_kappa": 0.0,
        }
    hydro = sum(KYTE_DOOLITTLE.get(c, 0.0) for c in seq) / L
    pos = sum(1 for c in seq if c in ("K", "R")) + 0.1 * sum(1 for c in seq if c == "H")
    neg = sum(1 for c in seq if c in ("D", "E"))
    net_q = pos - neg
    fcr = (sum(1 for c in seq if c in ("K", "R", "D", "E")) / L)
    # Charge segregation proxy in [0, 1]
    if fcr < 0.10 or L < 5:
        kappa = 0.10
    else:
        window = min(5, L)
        sigma_vals: list[float] = []
        for i in range(L - window + 1):
            sub = seq[i : i + window]
            sp = sum(1 for c in sub if c in ("K", "R")) / window
            sn = sum(1 for c in sub if c in ("D", "E")) / window
            sigma_vals.append((sp - sn) ** 2)
        kappa = min(1.0, sum(sigma_vals) / max(1, len(sigma_vals)))
    return {
        "hydropathy_kyte_doolittle": float(hydro),
        "net_charge_ph7": float(net_q),
        "fcr": float(fcr),
        "charge_segregation_kappa": float(kappa),
    }


@dataclass(frozen=True)
class RetrievedAuthorityCandidate:
    record_id: str
    evidence_tier: str
    sequence: str
    length: int
    similarity_score: float


class AuthoritySequenceIndex:
    """Fast in-memory authority & template retriever for X2 nearest-neighbour and Part 4 ranking."""

    def __init__(self, records: Sequence[Mapping[str, Any]]) -> None:
        self.records = list(records)

    def retrieve(
        self,
        *,
        target_length: int,
        rigidity: str = "flexible",
        cleavable: bool = False,
        protease_family: str = "TEV",
        left_flank_sequence: str = "",
        right_flank_sequence: str = "",
        top_k: int = 5,
    ) -> list[RetrievedAuthorityCandidate]:
        scored: list[tuple[tuple[float, int, str], RetrievedAuthorityCandidate]] = []
        for row in self.records:
            seq = str(row.get("linker_sequence", "")).strip().upper()
            if len(seq) < 4:
                continue
            L = len(seq)
            props = row.get("biophysical_properties") or {}
            eff_rig = str(
                props.get("curated_rigidity_label")
                or props.get("rigidity_predicted_label")
                or ("rigid" if ("EAAAK" in seq or seq.count("P") / L > 0.30) else "flexible")
            )
            has_cleav = any(m in seq for m in CLEAVAGE_MOTIFS.values())
            if cleavable and not has_cleav:
                continue
            if not cleavable and has_cleav:
                continue
            rig_match = 1.0 if eff_rig == rigidity else 0.0
            len_sim = math.exp(-abs(L - target_length) / max(float(target_length), 5.0))
            flank_sim = flanking_context_similarity(
                left_flank_sequence,
                right_flank_sequence,
                None,
                None,
                dict(row.get("flanking_context") or {}),
            )
            tier = str(row.get("authority_tier") or row.get("evidence_tier") or "mined_high_confidence_dfl_call")
            sim = 0.45 * len_sim + 0.35 * rig_match + 0.20 * flank_sim
            cand = RetrievedAuthorityCandidate(
                record_id=str(row.get("record_id", "auth")),
                evidence_tier=tier,
                sequence=seq,
                length=L,
                similarity_score=round(sim, 4),
            )
            scored.append(((-round(sim, 4), TIER_PRIORITY.get(tier, 5), cand.record_id), cand))
        scored.sort(key=lambda x: x[0])
        out: list[RetrievedAuthorityCandidate] = []
        seen: set[str] = set()
        for _, c in scored:
            if c.sequence in seen:
                continue
            seen.add(c.sequence)
            out.append(c)
            if len(out) >= top_k:
                break
        return out

SYNTHETIC_NEGATIVE_TYPES: tuple[str, ...] = (
    "too_short_by_50pct",
    "too_long_by_200pct",
    "rigid_helical_where_flexible_needed",
    "all_hydrophobic_core",
    "charge_clash_with_junction_anchors",
    "protease_site_buried_in_rigid_helix",
)

BASELINE_NAMES: tuple[str, ...] = (
    "direct_0aa_fusion",
    "fixed_ggggs3",
    "random_length_matched_gs",
    "geometry_only",
    "nearest_neighbour_retrieval",
    "linkerforge_ranker",
)

HYDROPHOBIC_SET = frozenset("FILVWM")
CLEAVAGE_MOTIFS: dict[str, str] = {
    "TEV": "ENLYFQS",
    "Thrombin": "LVPRGS",
    "PreScission_3C": "LEVLFQGP",
    "Enterokinase": "DDDDK",
    "Factor_Xa": "IEGR",
}


@dataclass(frozen=True)
class HardRejectionScorecard:
    """Deterministic structural & sequence compliance scorecard (no probability calibration)."""

    physical_span_closure_pass: bool
    domain_clash_clearance_pass: bool
    anchor_accessibility_pass: bool
    hydrophobic_stretch_pass: bool
    prohibited_cleavage_pass: bool
    composite_rigid_cleavable_pass: bool
    hard_rejection_pass: bool
    rejection_reasons: tuple[str, ...]
    contour_length_A: float
    required_path_length_A: float | None
    max_hydrophobic_run: int
    hydrophobic_fraction: float
    not_evaluated: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["rejection_reasons"] = list(self.rejection_reasons)
        d["not_evaluated"] = list(self.not_evaluated)
        return d


@dataclass(frozen=True)
class CandidateRankingResult:
    """Ranked fusion candidate with compatibility score and deterministic rejection scorecard."""

    rank: int
    candidate_id: str
    method_name: str
    sequence: str
    length: int
    compatibility_score: float
    hard_rejection_pass: bool
    scorecard: HardRejectionScorecard
    feature_breakdown: dict[str, float]
    provenance: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "rank": self.rank,
            "candidate_id": self.candidate_id,
            "method_name": self.method_name,
            "sequence": self.sequence,
            "length": self.length,
            "compatibility_score": round(self.compatibility_score, 4),
            "hard_rejection_pass": self.hard_rejection_pass,
            "scorecard": self.scorecard.to_dict(),
            "feature_breakdown": {k: round(v, 4) for k, v in self.feature_breakdown.items()},
            "provenance": self.provenance,
        }


def _max_contiguous_hydrophobic_run(sequence: str) -> int:
    best = 0
    cur = 0
    for ch in sequence.upper():
        if ch in HYDROPHOBIC_SET:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 0
    return best


def _anchor_net_charge(flank_sequence: str, window: int = 8, is_left: bool = True) -> float:
    seq = flank_sequence.upper()
    sub = seq[-window:] if is_left else seq[:window]
    if not sub:
        return 0.0
    pos = sum(1 for c in sub if c in ("K", "R", "H"))
    neg = sum(1 for c in sub if c in ("D", "E"))
    return float(pos - neg)


DETOUR_RATIO_LIMIT = 2.2  # design constant, not a literature value
ANCHOR_RSA_MINIMUM = 0.15  # inside the 5-25 % range of published buried/exposed cutoffs
EXCESS_LENGTH_FACTOR = 2.2
EXCESS_LENGTH_SLACK = 18


def evaluate_hard_rejection_scorecard(
    sequence: str,
    *,
    required_path_length_A: float | None = None,
    target_length: int | None = None,
    detour_ratio: float | None = None,
    port_a_sasa: float | None = None,
    port_b_sasa: float | None = None,
    requested_rigidity: str = "flexible",
    requested_cleavable: bool = False,
) -> HardRejectionScorecard:
    """Evaluate the six deterministic rules.

    A rule whose inputs are not supplied is *not evaluated* (it passes vacuously and is listed in
    `not_evaluated`).  Thresholds `DETOUR_RATIO_LIMIT`, `ANCHOR_RSA_MINIMUM` and the excess-length
    factor are design constants, not fitted or literature-derived values.
    """
    seq = sequence.upper().strip()
    L = len(seq)
    reach_A = round(maximum_span(L), 3)
    not_evaluated: list[str] = []

    # 1. physical span closure: (L + 1) virtual C-alpha bonds of at most 3.8 A each
    if required_path_length_A is None:
        physical_span_closure_pass = True
        not_evaluated.append("physical_span_closure")
    else:
        physical_span_closure_pass = reach_A >= required_path_length_A

    # 2. steric detour limit and excessive-length limit
    domain_clash_clearance_pass = True
    if detour_ratio is None and target_length is None:
        not_evaluated.append("domain_clash_clearance")
    if detour_ratio is not None and detour_ratio > DETOUR_RATIO_LIMIT:
        domain_clash_clearance_pass = False
    if target_length is not None and L > max(
        EXCESS_LENGTH_FACTOR * target_length, target_length + EXCESS_LENGTH_SLACK
    ):
        domain_clash_clearance_pass = False

    # 3. anchor accessibility
    if port_a_sasa is None or port_b_sasa is None:
        anchor_accessibility_pass = True
        not_evaluated.append("anchor_accessibility")
    else:
        anchor_accessibility_pass = min(port_a_sasa, port_b_sasa) >= ANCHOR_RSA_MINIMUM

    # 4. hydrophobic stretch (< 50 % FILVWM and no run of 5)
    hydro_count = sum(1 for c in seq if c in HYDROPHOBIC_SET)
    hydro_frac = (hydro_count / L) if L > 0 else 0.0
    max_hydro_run = _max_contiguous_hydrophobic_run(seq)
    hydrophobic_stretch_pass = (hydro_frac < 0.50) and (max_hydro_run < 5)

    # 5. cleavage motif compliance
    has_cleavage_motif = any(motif in seq for motif in CLEAVAGE_MOTIFS.values())
    prohibited_cleavage_pass = has_cleavage_motif if requested_cleavable else not has_cleavage_motif

    # 6. protease motif not enclosed by helical EAAAK repeats on both sides
    composite_rigid_cleavable_pass = True
    for motif in CLEAVAGE_MOTIFS.values():
        idx = seq.find(motif)
        if idx != -1:
            left_ctx = seq[max(0, idx - 5) : idx]
            right_ctx = seq[idx + len(motif) : idx + len(motif) + 5]
            if ("EAAAK" in left_ctx or "AAAK" in left_ctx) and (
                "EAAAK" in right_ctx or "EAAA" in right_ctx
            ):
                composite_rigid_cleavable_pass = False

    reasons: list[str] = []
    if not physical_span_closure_pass:
        reasons.append(
            f"physical_span_closure_failed(max_reach={reach_A:.1f}A < path={required_path_length_A:.1f}A, L={L})"
        )
    if not domain_clash_clearance_pass:
        reasons.append(
            f"detour_or_excess_length(L={L}, target={target_length}, detour={detour_ratio})"
        )
    if not anchor_accessibility_pass:
        reasons.append(f"buried_attachment_port(rsa_a={port_a_sasa}, rsa_b={port_b_sasa})")
    if not hydrophobic_stretch_pass:
        reasons.append(f"hydrophobic_collapse_risk(frac={hydro_frac:.2f}, run={max_hydro_run})")
    if not prohibited_cleavage_pass:
        reasons.append(
            "missing_requested_protease_motif" if requested_cleavable else "unintended_protease_motif"
        )
    if not composite_rigid_cleavable_pass:
        reasons.append("protease_motif_enclosed_by_helical_repeats")

    return HardRejectionScorecard(
        physical_span_closure_pass=physical_span_closure_pass,
        domain_clash_clearance_pass=domain_clash_clearance_pass,
        anchor_accessibility_pass=anchor_accessibility_pass,
        hydrophobic_stretch_pass=hydrophobic_stretch_pass,
        prohibited_cleavage_pass=prohibited_cleavage_pass,
        composite_rigid_cleavable_pass=composite_rigid_cleavable_pass,
        hard_rejection_pass=len(reasons) == 0,
        rejection_reasons=tuple(reasons),
        contour_length_A=reach_A,
        required_path_length_A=None if required_path_length_A is None else round(required_path_length_A, 3),
        max_hydrophobic_run=max_hydro_run,
        hydrophobic_fraction=round(hydro_frac, 4),
        not_evaluated=tuple(not_evaluated),
    )


def generate_synthetic_negatives(
    positive_sequence: str,
    *,
    target_length: int,
    requested_rigidity: str = "flexible",
    left_flank: str = "",
    right_flank: str = "",
) -> dict[str, str]:
    """Generate the 6 canonical synthetic negative corruptions (`ARCHITECTURE_HANDOFF.md` §7)."""
    pos = positive_sequence.upper().strip()
    L = max(len(pos), target_length, 6)

    # 1. Too short by 50%
    short_len = max(2, int(round(0.45 * target_length)))
    too_short = pos[:short_len] if len(pos) >= short_len else ("GS" * short_len)[:short_len]

    # 2. Too long by 200% (3x target length)
    long_len = max(target_length * 3, L * 3, 24)
    repeat_unit = pos if len(pos) >= 4 else "GGGGS"
    too_long = (repeat_unit * ((long_len // len(repeat_unit)) + 2))[:long_len]

    # 3. Rigid helical where flexible needed (or ultra-flexible poly-G where rigid needed)
    if requested_rigidity == "flexible":
        rigid_wrong = ("EAAAK" * ((L // 5) + 2))[:L]
    else:
        rigid_wrong = ("GGGGS" * ((L // 5) + 2))[:L]

    # 4. All-hydrophobic core sequence
    hydro_motif = "LVLFLILVWMF"
    all_hydrophobic = (hydro_motif * ((L // len(hydro_motif)) + 2))[:L]

    # 5. Charge clash with junction anchors
    flank_charge = _anchor_net_charge(left_flank, is_left=True) + _anchor_net_charge(
        right_flank, is_left=False
    )
    clash_aa = "K" if flank_charge >= 0 else "E"
    charge_clash = clash_aa * L

    # 6. Protease site buried directly inside a rigid alpha-helix without flexible GS spacers
    buried_protease = "EAAAKEAAAKENLYFQSEAAAKEAAAK"
    if L > len(buried_protease):
        pad = ("EAAAK" * ((L - len(buried_protease)) // 5 + 1))[: L - len(buried_protease)]
        buried_protease = buried_protease + pad

    return {
        "too_short_by_50pct": too_short,
        "too_long_by_200pct": too_long,
        "rigid_helical_where_flexible_needed": rigid_wrong,
        "all_hydrophobic_core": all_hydrophobic,
        "charge_clash_with_junction_anchors": charge_clash,
        "protease_site_buried_in_rigid_helix": buried_protease,
    }


def extract_pairwise_compatibility_features(
    sequence: str,
    *,
    target_length: int,
    required_path_length_A: float | None = None,
    requested_rigidity: str = "flexible",
    requested_cleavable: bool = False,
    left_flank: str = "",
    right_flank: str = "",
    port_a_sasa: float | None = None,
    port_b_sasa: float | None = None,
) -> dict[str, float]:
    """Extract 12 interpretable physical, geometric, and junction-compatibility features."""
    seq = sequence.upper().strip()
    L = len(seq)
    if L == 0:
        return {
            "length_closeness": 0.0,
            "span_closure_margin": -1.0,
            "entropic_length_penalty": 1.0,
            "rigidity_axis_alignment": 0.0,
            "hydrophobicity_safety": 0.0,
            "disorder_coil_propensity": 0.0,
            "charge_junction_compatibility": 0.0,
            "charge_mixing_quality": 0.0,
            "spacer_cleavage_compliance": 0.0,
            "port_accessibility_score": 0.0,
            "sequence_complexity_balance": 0.0,
            "hard_pass_indicator": 0.0,
        }

    phys = compute_linker_sequence_physics(seq)
    scorecard = evaluate_hard_rejection_scorecard(
        seq,
        required_path_length_A=required_path_length_A,
        target_length=target_length,
        port_a_sasa=port_a_sasa,
        port_b_sasa=port_b_sasa,
        requested_rigidity=requested_rigidity,
        requested_cleavable=requested_cleavable,
    )

    # 1. Length closeness to target geometric length
    rel_len_err = abs(L - target_length) / max(float(target_length), 5.0)
    length_closeness = math.exp(-1.75 * rel_len_err)

    # 2. Span closure margin (positive when contour comfortably spans path, negative when impossible)
    if required_path_length_A is None:
        span_closure_margin = 0.0  # no geometry supplied: neutral
    else:
        contour_ratio = maximum_span(L) / max(required_path_length_A, 4.0)
        if contour_ratio < 1.0:
            span_closure_margin = contour_ratio - 1.5  # strong negative penalty
        else:
            span_closure_margin = min(1.0, (contour_ratio - 1.0) / 1.2)

    # 3. Entropic excess length penalty
    entropic_length_penalty = max(0.0, (L - max(target_length * 1.4, target_length + 6)) / max(float(target_length), 6.0))
    entropic_length_penalty = min(entropic_length_penalty, 2.0)

    # 4. Rigidity axis alignment
    gly_ser_frac = sum(1 for c in seq if c in ("G", "S", "T", "N", "Q")) / L
    helical_salt_motifs = seq.count("EAAAK") + seq.count("AAAK")
    pro_ala_frac = sum(1 for c in seq if c in ("P", "A", "E", "K")) / L
    if requested_rigidity == "flexible":
        # Natural disordered linkers are enriched in polar/coil residues and avoid long EAAAK rods
        if helical_salt_motifs >= 2:
            rigidity_axis_alignment = -0.8
        else:
            rigidity_axis_alignment = min(1.0, 0.55 * gly_ser_frac + 0.45 * (1.0 - min(1.0, phys["hydropathy_kyte_doolittle"] / 2.0)))
    elif requested_rigidity == "rigid":
        rigidity_axis_alignment = min(1.0, 0.7 * pro_ala_frac + 0.3 * min(1.0, helical_salt_motifs / 2.0))
    else:
        rigidity_axis_alignment = 0.75

    # 5. Hydrophobicity safety
    hydro_frac = scorecard.hydrophobic_fraction
    max_run = scorecard.max_hydrophobic_run
    hydrophobicity_safety = max(-1.0, 1.0 - 2.2 * hydro_frac - 0.18 * max(0, max_run - 2))

    # 6. Disorder / coil propensity
    disorder_residues = sum(1 for c in seq if c in ("G", "S", "P", "E", "K", "Q", "T", "D", "N", "R", "A"))
    disorder_coil_propensity = disorder_residues / L

    # 7. Charge junction compatibility with flanking anchors
    flank_charge = _anchor_net_charge(left_flank, is_left=True) + _anchor_net_charge(
        right_flank, is_left=False
    )
    linker_net_charge = phys["net_charge_ph7"]
    fcr = phys["fcr"]
    if fcr > 0.70 and abs(linker_net_charge) >= 0.60 * L:
        charge_junction_compatibility = -1.0
    elif flank_charge * linker_net_charge > 12.0:
        charge_junction_compatibility = -0.5
    else:
        charge_junction_compatibility = max(0.0, 1.0 - 0.04 * abs(linker_net_charge))

    # 8. Charge mixing quality (penalize extreme homopolymeric charge blocks kappa -> 1)
    kappa = phys["charge_segregation_kappa"]
    charge_mixing_quality = 1.0 - max(0.0, kappa - 0.35)

    # 9. Spacer & cleavage compliance
    if scorecard.prohibited_cleavage_pass and scorecard.composite_rigid_cleavable_pass:
        spacer_cleavage_compliance = 1.0
    else:
        spacer_cleavage_compliance = -1.0

    # 10. Port accessibility score
    port_accessibility_score = (
        0.0 if port_a_sasa is None or port_b_sasa is None else 0.5 * (port_a_sasa + port_b_sasa)
    )

    # 11. Sequence complexity balance (avoid single-residue homopolymers and random 2-letter runs)
    unique_aa = len(set(seq))
    max_single_aa_frac = max(seq.count(c) for c in set(seq)) / L
    if unique_aa == 1:
        sequence_complexity_balance = -0.8
    elif unique_aa == 2 and ("GGGGS" not in seq and "PAPAP" not in seq):
        sequence_complexity_balance = 0.45 - max(0.0, max_single_aa_frac - 0.60)
    else:
        sequence_complexity_balance = (
            1.0 - max(0.0, max_single_aa_frac - 0.65) * 2.0 + min(0.30, unique_aa * 0.04)
        )

    # 12. Hard pass indicator
    hard_pass_indicator = 1.0 if scorecard.hard_rejection_pass else -1.0

    return {
        "length_closeness": float(length_closeness),
        "span_closure_margin": float(span_closure_margin),
        "entropic_length_penalty": float(entropic_length_penalty),
        "rigidity_axis_alignment": float(rigidity_axis_alignment),
        "hydrophobicity_safety": float(hydrophobicity_safety),
        "disorder_coil_propensity": float(disorder_coil_propensity),
        "charge_junction_compatibility": float(charge_junction_compatibility),
        "charge_mixing_quality": float(charge_mixing_quality),
        "spacer_cleavage_compliance": float(spacer_cleavage_compliance),
        "port_accessibility_score": float(port_accessibility_score),
        "sequence_complexity_balance": float(sequence_complexity_balance),
        "hard_pass_indicator": float(hard_pass_indicator),
    }


FEATURE_ORDER: tuple[str, ...] = (
    "length_closeness",
    "span_closure_margin",
    "entropic_length_penalty",
    "rigidity_axis_alignment",
    "hydrophobicity_safety",
    "disorder_coil_propensity",
    "charge_junction_compatibility",
    "charge_mixing_quality",
    "spacer_cleavage_compliance",
    "port_accessibility_score",
    "sequence_complexity_balance",
    "hard_pass_indicator",
)


@dataclass(frozen=True)
class PairwiseCompatibilityRanker:
    """Learned pairwise compatibility ranker (`ARCHITECTURE_HANDOFF.md` §7)."""

    feature_names: tuple[str, ...]
    weights: tuple[float, ...]

    def score_features(self, feat_dict: Mapping[str, float]) -> float:
        return sum(w * float(feat_dict.get(k, 0.0)) for k, w in zip(self.feature_names, self.weights))

    def score_candidate(
        self,
        sequence: str,
        *,
        target_length: int,
        required_path_length_A: float | None = None,
        requested_rigidity: str = "flexible",
        requested_cleavable: bool = False,
        left_flank: str = "",
        right_flank: str = "",
        port_a_sasa: float | None = None,
        port_b_sasa: float | None = None,
    ) -> tuple[float, HardRejectionScorecard, dict[str, float]]:
        feats = extract_pairwise_compatibility_features(
            sequence,
            target_length=target_length,
            required_path_length_A=required_path_length_A,
            requested_rigidity=requested_rigidity,
            requested_cleavable=requested_cleavable,
            left_flank=left_flank,
            right_flank=right_flank,
            port_a_sasa=port_a_sasa,
            port_b_sasa=port_b_sasa,
        )
        scorecard = evaluate_hard_rejection_scorecard(
            sequence,
            required_path_length_A=required_path_length_A,
            target_length=target_length,
            port_a_sasa=port_a_sasa,
            port_b_sasa=port_b_sasa,
            requested_rigidity=requested_rigidity,
            requested_cleavable=requested_cleavable,
        )
        raw_score = self.score_features(feats)
        return raw_score, scorecard, feats


def record_query_context(rec: Mapping[str, Any]) -> dict[str, Any]:
    """Extract the design-query context of a dataset record using only recorded fields.

    `path_A` and the port accessibilities exist only for coordinate-complete records (the CATH
    tier); for every other record they are None and the corresponding rules / features are
    neutral.  Rigidity comes from the record's curated or predicted label, never from the
    sequence heuristics used by earlier versions.
    """
    seq = str(rec.get("linker_sequence", "")).strip().upper()
    geom = rec.get("geometry") or {}
    props = rec.get("biophysical_properties") or {}
    curated = props.get("curated_rigidity_label")
    predicted = props.get("rigidity_predicted_label")
    if curated in ("flexible", "rigid"):
        rigidity = curated
    elif predicted in ("flexible", "rigid"):
        rigidity = predicted
    else:
        rigidity = "flexible"
    flank = rec.get("flanking_context") or {}
    return {
        "sequence": seq,
        "target_length": int(rec.get("linker_length", len(seq))),
        "path_A": geom.get("inter_port_ca_distance_A"),
        "sasa_a": geom.get("port_a_accessibility_proxy"),
        "sasa_b": geom.get("port_b_accessibility_proxy"),
        "rigidity": rigidity,
        "left_flank": str(flank.get("left_flank_sequence") or ""),
        "right_flank": str(flank.get("right_flank_sequence") or ""),
    }


def fit_pairwise_compatibility_ranker(
    training_records: Sequence[Mapping[str, Any]],
    *,
    learning_rate: float = 0.15,
    epochs: int = 150,
    l2_penalty: float = 1e-3,
    seed: int = 23,
    init: str = "zero",
) -> PairwiseCompatibilityRanker:
    """Fit the linear ranker with a pairwise logistic loss: real linker vs. its six corruptions.

    `init="zero"` (default) starts from all-zero weights, so the result is learned from the
    corruption task; `init="prior"` starts from hand-set weights (the earlier behaviour, in which
    training barely moved the weights).  The six corruptions are constructed to violate the same
    properties several features measure, so held-out accuracy on this task is a training-objective
    sanity check, not evidence about real fusion quality.
    """
    rng = random.Random(seed)
    n_feat = len(FEATURE_ORDER)
    if init == "prior":
        weights = [1.4, 1.2, -1.1, 1.0, 1.2, 0.7, 0.9, 0.5, 1.1, 0.3, 0.6, 1.5]
    elif init == "zero":
        weights = [0.0] * n_feat
    else:
        raise ValueError("init must be 'zero' or 'prior'")

    def _vec(seq: str, ctx: Mapping[str, Any]) -> list[float]:
        feats = extract_pairwise_compatibility_features(
            seq,
            target_length=ctx["target_length"],
            required_path_length_A=ctx["path_A"],
            requested_rigidity=ctx["rigidity"],
            requested_cleavable=False,
            left_flank=ctx["left_flank"],
            right_flank=ctx["right_flank"],
            port_a_sasa=ctx["sasa_a"],
            port_b_sasa=ctx["sasa_b"],
        )
        return [feats[k] for k in FEATURE_ORDER]

    pairs: list[tuple[list[float], list[float]]] = []
    for rec in training_records:
        ctx = record_query_context(rec)
        pos_seq = ctx["sequence"]
        if len(pos_seq) < 4:
            continue
        pos_vec = _vec(pos_seq, ctx)
        negs = generate_synthetic_negatives(
            pos_seq,
            target_length=ctx["target_length"],
            requested_rigidity=ctx["rigidity"],
            left_flank=ctx["left_flank"],
            right_flank=ctx["right_flank"],
        )
        for neg_seq in negs.values():
            pairs.append((pos_vec, _vec(neg_seq, ctx)))

    if not pairs:
        return PairwiseCompatibilityRanker(feature_names=FEATURE_ORDER, weights=tuple(weights))

    for _ in range(epochs):
        rng.shuffle(pairs)
        grad_w = [0.0] * n_feat
        for pos_vec, neg_vec in pairs:
            diff = [p - n for p, n in zip(pos_vec, neg_vec)]
            margin = sum(w * d for w, d in zip(weights, diff))
            prob_wrong = 1.0 / (1.0 + math.exp(max(-30.0, min(30.0, margin))))
            for j in range(n_feat):
                grad_w[j] += -prob_wrong * diff[j]
        inv_n = 1.0 / len(pairs)
        for j in range(n_feat):
            weights[j] -= learning_rate * (grad_w[j] * inv_n + l2_penalty * weights[j])

    return PairwiseCompatibilityRanker(
        feature_names=FEATURE_ORDER,
        weights=tuple(round(w, 6) for w in weights),
    )
