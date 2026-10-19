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
    required_path_length_A: float
    max_hydrophobic_run: int
    hydrophobic_fraction: float

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["rejection_reasons"] = list(self.rejection_reasons)
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


def evaluate_hard_rejection_scorecard(
    sequence: str,
    *,
    required_path_length_A: float,
    target_length: int,
    detour_ratio: float = 1.15,
    port_a_sasa: float = 0.45,
    port_b_sasa: float = 0.45,
    requested_rigidity: str = "flexible",
    requested_cleavable: bool = False,
) -> HardRejectionScorecard:
    """Evaluate the 6 deterministic hard-rejection rules on a candidate linker sequence."""
    seq = sequence.upper().strip()
    L = len(seq)
    contour_A = round(L * 3.8, 3)

    # 1. Physical span closure (3.8 A per residue maximum extended C-alpha contour)
    # Require at least contour_A >= required_path_length_A and L >= 1 when path > 3.8 A
    min_contour_needed = max(required_path_length_A, min(target_length * 2.1, required_path_length_A))
    physical_span_closure_pass = (contour_A >= min_contour_needed) and (
        L >= max(2, int(math.floor(target_length * 0.60))) if required_path_length_A > 4.0 else True
    )

    # 2. Domain clash & excessive loop detour clearance
    # Direct 0-aa fusion or >2.5x target length causes domain clash or excessive entropic loop
    domain_clash_clearance_pass = (
        (detour_ratio <= 2.2)
        and (L > 0 or required_path_length_A <= 3.8)
        and (L <= max(target_length * 2.2, target_length + 18))
    )

    # 3. Anchor accessibility
    anchor_accessibility_pass = (port_a_sasa >= 0.15) and (port_b_sasa >= 0.15)

    # 4. Hydrophobic stretch check (< 50% FILVWM and max contiguous run < 5)
    hydro_count = sum(1 for c in seq if c in HYDROPHOBIC_SET)
    hydro_frac = (hydro_count / L) if L > 0 else 0.0
    max_hydro_run = _max_contiguous_hydrophobic_run(seq)
    hydrophobic_stretch_pass = (hydro_frac < 0.50) and (max_hydro_run < 5)

    # 5. Cleavage site compliance
    has_cleavage_motif = any(motif in seq for motif in CLEAVAGE_MOTIFS.values())
    if requested_cleavable:
        prohibited_cleavage_pass = has_cleavage_motif
    else:
        prohibited_cleavage_pass = not has_cleavage_motif

    # 6. Composite rigid + cleavable spacer rule
    # Protease recognition motif must not be directly flanked on both sides by rigid EAAAK/AAA helices
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
            f"physical_span_closure_failed(contour={contour_A:.1f}A < required={min_contour_needed:.1f}A, L={L})"
        )
    if not domain_clash_clearance_pass:
        reasons.append(
            f"domain_clash_or_entropic_excess(L={L}, target={target_length}, detour={detour_ratio:.2f})"
        )
    if not anchor_accessibility_pass:
        reasons.append(
            f"buried_attachment_port(port_a_sasa={port_a_sasa:.2f}, port_b_sasa={port_b_sasa:.2f})"
        )
    if not hydrophobic_stretch_pass:
        reasons.append(
            f"hydrophobic_collapse_risk(hydro_frac={hydro_frac:.2f}, max_run={max_hydro_run})"
        )
    if not prohibited_cleavage_pass:
        if requested_cleavable:
            reasons.append("missing_requested_protease_cleavage_motif")
        else:
            reasons.append("unintended_protease_cleavage_motif_present")
    if not composite_rigid_cleavable_pass:
        reasons.append("protease_site_buried_inside_rigid_alpha_helix_without_gs_spacers")

    hard_pass = len(reasons) == 0
    return HardRejectionScorecard(
        physical_span_closure_pass=physical_span_closure_pass,
        domain_clash_clearance_pass=domain_clash_clearance_pass,
        anchor_accessibility_pass=anchor_accessibility_pass,
        hydrophobic_stretch_pass=hydrophobic_stretch_pass,
        prohibited_cleavage_pass=prohibited_cleavage_pass,
        composite_rigid_cleavable_pass=composite_rigid_cleavable_pass,
        hard_rejection_pass=hard_pass,
        rejection_reasons=tuple(reasons),
        contour_length_A=contour_A,
        required_path_length_A=round(required_path_length_A, 3),
        max_hydrophobic_run=max_hydro_run,
        hydrophobic_fraction=round(hydro_frac, 4),
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
    required_path_length_A: float,
    requested_rigidity: str = "flexible",
    requested_cleavable: bool = False,
    left_flank: str = "",
    right_flank: str = "",
    port_a_sasa: float = 0.45,
    port_b_sasa: float = 0.45,
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
            "port_accessibility_score": min(port_a_sasa, port_b_sasa),
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
    contour_ratio = (L * 3.8) / max(required_path_length_A, 4.0)
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
    port_accessibility_score = 0.5 * (port_a_sasa + port_b_sasa)

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
        required_path_length_A: float,
        requested_rigidity: str = "flexible",
        requested_cleavable: bool = False,
        left_flank: str = "",
        right_flank: str = "",
        port_a_sasa: float = 0.45,
        port_b_sasa: float = 0.45,
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


def fit_pairwise_compatibility_ranker(
    training_records: Sequence[Mapping[str, Any]],
    *,
    learning_rate: float = 0.08,
    epochs: int = 80,
    l2_penalty: float = 1e-3,
    seed: int = 23,
) -> PairwiseCompatibilityRanker:
    """Train `PairwiseCompatibilityRanker` using pairwise logistic ranking loss over
    natural/mined positive linkers vs. the 6 synthetic negative corruptions.
    """
    rng = random.Random(seed)
    n_feat = len(FEATURE_ORDER)
    weights = [
        1.4,   # length_closeness
        1.2,   # span_closure_margin
        -1.1,  # entropic_length_penalty
        1.0,   # rigidity_axis_alignment
        1.2,   # hydrophobicity_safety
        0.7,   # disorder_coil_propensity
        0.9,   # charge_junction_compatibility
        0.5,   # charge_mixing_quality
        1.1,   # spacer_cleavage_compliance
        0.3,   # port_accessibility_score
        0.6,   # sequence_complexity_balance
        1.5,   # hard_pass_indicator
    ]

    pairs: list[tuple[list[float], list[float]]] = []
    for rec in training_records:
        pos_seq = str(rec.get("linker_sequence", "")).strip().upper()
        if len(pos_seq) < 4:
            continue
        target_len = int(rec.get("linker_length", len(pos_seq)))
        geom = rec.get("geometry") or {}
        req_path = float(
            geom.get("estimated_end_to_end_ca_distance_A")
            or round(1.95 * target_len, 2)
        )
        flank = rec.get("flanking_context") or {}
        left_flank = str(flank.get("left_flank_sequence", ""))
        right_flank = str(flank.get("right_flank_sequence", ""))
        pa_sasa = float(geom.get("port_a_accessibility_proxy") or 0.45)
        pb_sasa = float(geom.get("port_b_accessibility_proxy") or 0.45)
        rigidity = "rigid" if ("EAAAK" in pos_seq or pos_seq.count("P") / len(pos_seq) > 0.30) else "flexible"

        pos_feats_dict = extract_pairwise_compatibility_features(
            pos_seq,
            target_length=target_len,
            required_path_length_A=req_path,
            requested_rigidity=rigidity,
            requested_cleavable=False,
            left_flank=left_flank,
            right_flank=right_flank,
            port_a_sasa=pa_sasa,
            port_b_sasa=pb_sasa,
        )
        pos_vec = [pos_feats_dict[k] for k in FEATURE_ORDER]

        negs = generate_synthetic_negatives(
            pos_seq,
            target_length=target_len,
            requested_rigidity=rigidity,
            left_flank=left_flank,
            right_flank=right_flank,
        )
        for neg_seq in negs.values():
            neg_feats_dict = extract_pairwise_compatibility_features(
                neg_seq,
                target_length=target_len,
                required_path_length_A=req_path,
                requested_rigidity=rigidity,
                requested_cleavable=False,
                left_flank=left_flank,
                right_flank=right_flank,
                port_a_sasa=pa_sasa,
                port_b_sasa=pb_sasa,
            )
            neg_vec = [neg_feats_dict[k] for k in FEATURE_ORDER]
            pairs.append((pos_vec, neg_vec))

    if not pairs:
        return PairwiseCompatibilityRanker(feature_names=FEATURE_ORDER, weights=tuple(weights))

    for _ in range(epochs):
        rng.shuffle(pairs)
        grad_w = [0.0] * n_feat
        for pos_vec, neg_vec in pairs:
            diff = [p - n for p, n in zip(pos_vec, neg_vec)]
            margin = sum(w * d for w, d in zip(weights, diff))
            # Pairwise logistic loss: log(1 + exp(-margin))
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


def evaluate_six_baselines_benchmark(
    records: Sequence[Mapping[str, Any]],
    *,
    n_folds: int = 5,
    seed: int = 23,
) -> dict[str, Any]:
    """Evaluate all 6 required baselines (`ARCHITECTURE_HANDOFF.md` §7) + pairwise synthetic
    negative discrimination on cluster-disjoint (`homology_cluster_id`) splits.
    """
    clusters = sorted({str(r.get("homology_cluster_id", r.get("record_id", ""))) for r in records})
    rng = random.Random(seed)
    rng.shuffle(clusters)
    cluster_to_fold = {cid: idx % n_folds for idx, cid in enumerate(clusters)}

    method_records: dict[str, list[dict[str, float]]] = {m: [] for m in BASELINE_NAMES}
    corruption_wins: dict[str, list[float]] = {c: [] for c in SYNTHETIC_NEGATIVE_TYPES}
    natural_vs_corruptions_top1: list[float] = []
    natural_vs_corruptions_mrr: list[float] = []

    for fold in range(n_folds):
        train_recs = [
            r
            for r in records
            if cluster_to_fold[str(r.get("homology_cluster_id", r.get("record_id", "")))] != fold
        ]
        test_recs = [
            r
            for r in records
            if cluster_to_fold[str(r.get("homology_cluster_id", r.get("record_id", "")))] == fold
        ]
        if not train_recs or not test_recs:
            continue

        ranker = fit_pairwise_compatibility_ranker(train_recs, seed=seed + fold)
        auth_index = AuthoritySequenceIndex(train_recs)

        for rec_idx, rec in enumerate(test_recs):
            true_seq = str(rec.get("linker_sequence", "")).strip().upper()
            if len(true_seq) < 4:
                continue
            target_len = int(rec.get("linker_length", len(true_seq)))
            geom = rec.get("geometry") or {}
            req_path = float(
                geom.get("estimated_end_to_end_ca_distance_A")
                or round(1.95 * target_len, 2)
            )
            flank = rec.get("flanking_context") or {}
            left_flank = str(flank.get("left_flank_sequence", ""))
            right_flank = str(flank.get("right_flank_sequence", ""))
            pa_sasa = float(geom.get("port_a_accessibility_proxy") or 0.45)
            pb_sasa = float(geom.get("port_b_accessibility_proxy") or 0.45)
            rigidity = (
                "rigid"
                if ("EAAAK" in true_seq or true_seq.count("P") / len(true_seq) > 0.30)
                else "flexible"
            )

            # Score natural linker vs 6 synthetic corruptions
            pos_score, _, _ = ranker.score_candidate(
                true_seq,
                target_length=target_len,
                required_path_length_A=req_path,
                requested_rigidity=rigidity,
                left_flank=left_flank,
                right_flank=right_flank,
                port_a_sasa=pa_sasa,
                port_b_sasa=pb_sasa,
            )
            negs = generate_synthetic_negatives(
                true_seq,
                target_length=target_len,
                requested_rigidity=rigidity,
                left_flank=left_flank,
                right_flank=right_flank,
            )
            neg_scores: list[float] = []
            for c_type, neg_seq in negs.items():
                n_score, _, _ = ranker.score_candidate(
                    neg_seq,
                    target_length=target_len,
                    required_path_length_A=req_path,
                    requested_rigidity=rigidity,
                    left_flank=left_flank,
                    right_flank=right_flank,
                    port_a_sasa=pa_sasa,
                    port_b_sasa=pb_sasa,
                )
                neg_scores.append(n_score)
                corruption_wins[c_type].append(1.0 if pos_score > n_score else 0.0)

            better_negs = sum(1 for ns in neg_scores if ns >= pos_score)
            rank_pos = 1 + better_negs
            natural_vs_corruptions_top1.append(1.0 if rank_pos == 1 else 0.0)
            natural_vs_corruptions_mrr.append(1.0 / rank_pos)

            # Build the 6 baseline candidate sequences for this query
            local_rng = random.Random(seed * 1000 + fold * 100 + rec_idx)
            # 1. Direct 0-aa fusion
            seq_direct = ""
            # 2. Fixed (GGGGS)3
            seq_ggggs3 = "GGGGSGGGGSGGGGS"
            # 3. Random length-matched GS
            seq_rand_gs = "".join(local_rng.choice(("G", "S")) for _ in range(target_len))
            # 4. Geometry-only (length-matched homopolymer/unfiltered template ignoring flank/rigidity)
            seq_geom_only = ("AAAGGG" * ((target_len // 6) + 2))[:target_len]
            # 5. Nearest-neighbour retrieval (X2)
            x2_hits = auth_index.retrieve(
                target_length=target_len,
                rigidity=rigidity,
                cleavable=False,
                left_flank_sequence=left_flank,
                right_flank_sequence=right_flank,
                top_k=4,
            )
            seq_x2 = x2_hits[0].sequence if x2_hits else seq_ggggs3

            # 6. LinkerForge ranker over Part 2 template candidates + X2 authority retrieval pool
            if rigidity == "rigid":
                tmpl_seqs = [
                    ("EAAAK" * ((target_len // 5) + 2))[:target_len],
                    ("PAPAP" * ((target_len // 5) + 2))[:target_len],
                    "AEAAAKEAAAKEAAAKA"[: max(10, target_len)],
                ]
            else:
                tmpl_seqs = [
                    ("GGGGS" * ((target_len // 5) + 2))[:target_len],
                    ("GSAGSAAGSGEF" * ((target_len // 12) + 2))[:target_len],
                    ("KESGSVSSEQLAQFRSLD" * ((target_len // 18) + 2))[:target_len],
                ]
            pool_seqs = tmpl_seqs + [h.sequence for h in x2_hits]
            if not pool_seqs:
                pool_seqs = [seq_ggggs3]
            best_lf_seq = pool_seqs[0]
            best_lf_tuple = (-1e9, -1e9)
            for cand_seq in pool_seqs:
                c_score, c_card, _ = ranker.score_candidate(
                    cand_seq,
                    target_length=target_len,
                    required_path_length_A=req_path,
                    requested_rigidity=rigidity,
                    left_flank=left_flank,
                    right_flank=right_flank,
                    port_a_sasa=pa_sasa,
                    port_b_sasa=pb_sasa,
                )
                cand_tuple = (1.0 if c_card.hard_rejection_pass else 0.0, c_score)
                if cand_tuple > best_lf_tuple:
                    best_lf_tuple = cand_tuple
                    best_lf_seq = cand_seq

            baseline_seqs = {
                "direct_0aa_fusion": seq_direct,
                "fixed_ggggs3": seq_ggggs3,
                "random_length_matched_gs": seq_rand_gs,
                "geometry_only": seq_geom_only,
                "nearest_neighbour_retrieval": seq_x2,
                "linkerforge_ranker": best_lf_seq,
            }

            for m_name, b_seq in baseline_seqs.items():
                b_score, b_card, b_feats = ranker.score_candidate(
                    b_seq,
                    target_length=target_len,
                    required_path_length_A=req_path,
                    requested_rigidity=rigidity,
                    left_flank=left_flank,
                    right_flank=right_flank,
                    port_a_sasa=pa_sasa,
                    port_b_sasa=pb_sasa,
                )
                win_vs_negs = (
                    sum(1.0 if b_score > ns else 0.0 for ns in neg_scores) / len(neg_scores)
                    if neg_scores
                    else 0.0
                )
                method_records[m_name].append(
                    {
                        "hard_rejection_pass": 1.0 if b_card.hard_rejection_pass else 0.0,
                        "compatibility_score": b_score,
                        "length_within_2aa": 1.0 if abs(len(b_seq) - target_len) <= 2 else 0.0,
                        "rigidity_axis_alignment": b_feats["rigidity_axis_alignment"],
                        "pairwise_win_rate_vs_synthetic_negatives": win_vs_negs,
                    }
                )

    def _agg(rows: list[dict[str, float]]) -> dict[str, float]:
        if not rows:
            return {}
        keys = list(rows[0].keys())
        return {
            f"mean_{k}": round(sum(r[k] for r in rows) / len(rows), 4)
            for k in keys
        }

    return {
        "n_records_evaluated": len(natural_vs_corruptions_top1),
        "n_folds": n_folds,
        "split_field": "homology_cluster_id",
        "natural_vs_6_synthetic_corruptions": {
            "top1_natural_rank_accuracy": round(
                sum(natural_vs_corruptions_top1) / max(1, len(natural_vs_corruptions_top1)), 4
            ),
            "mrr_natural_rank": round(
                sum(natural_vs_corruptions_mrr) / max(1, len(natural_vs_corruptions_mrr)), 4
            ),
            "per_corruption_pairwise_win_rate": {
                k: round(sum(v) / max(1, len(v)), 4) for k, v in corruption_wins.items()
            },
        },
        "six_baselines_comparison": {m: _agg(rows) for m, rows in method_records.items()},
    }


ONE_TO_THREE: dict[str, str] = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def build_fusion_3d_pdb_and_pymol_script(
    *,
    domain_a_name: str,
    domain_a_sequence: str,
    linker_sequence: str,
    domain_b_name: str,
    domain_b_sequence: str,
    inter_port_distance_A: float = 24.8,
) -> dict[str, Any]:
    """Generate a complete 3D backbone PDB file (`N, CA, C, O`) and ready-to-run PyMOL `.pml`
    script for a two-domain fusion construct (`Domain A — linker — Domain B`), suitable for
    both Mol* (`molstar`) and PyMOL / 3Dmol WebGL visualization.

    Models Domain A (`sfGFP`, 11-stranded beta-barrel cylinder on the left), the 3D collision-free
    curved spline trajectory of the designed linker spanning `inter_port_distance_A`, and Domain B
    (`mCherry`, 11-stranded beta-barrel cylinder on the right).
    """
    seq_a = domain_a_sequence.upper().strip()
    seq_l = linker_sequence.upper().strip()
    seq_b = domain_b_sequence.upper().strip()

    pdb_lines: list[str] = [
        f"HEADER    FUSION PROTEIN                          19-OCT-26   LFG1",
        f"TITLE     LINKERFORGE END-TO-END FUSION: {domain_a_name} - LINKER ({len(seq_l)} AA) - {domain_b_name}",
        f"REMARK   1 DOMAIN_A: {domain_a_name} (RESIDUES 1-{len(seq_a)})",
        f"REMARK   2 LINKER:   {seq_l or 'DIRECT_0AA'} (RESIDUES {len(seq_a)+1}-{len(seq_a)+len(seq_l)})",
        f"REMARK   3 DOMAIN_B: {domain_b_name} (RESIDUES {len(seq_a)+len(seq_l)+1}-{len(seq_a)+len(seq_l)+len(seq_b)})",
    ]

    atom_serial = 1
    res_serial = 1

    def _emit_residue(res_char: str, chain_id: str, res_num: int, ca_xyz: tuple[float, float, float], b_factor: float) -> None:
        nonlocal atom_serial
        res3 = ONE_TO_THREE.get(res_char, "GLY")
        cx, cy, cz = ca_xyz
        atoms = [
            ("N", (cx - 1.22, cy - 0.45, cz - 0.30)),
            ("CA", (cx, cy, cz)),
            ("C", (cx + 1.25, cy + 0.42, cz + 0.28)),
            ("O", (cx + 1.45, cy + 1.58, cz + 0.52)),
        ]
        for atom_name, (ax, ay, az) in atoms:
            element = atom_name[0]
            pdb_lines.append(
                f"ATOM  {atom_serial:5d}  {atom_name:<3s} {res3:>3s} {chain_id:1s}{res_num:4d}    "
                f"{ax:8.3f}{ay:8.3f}{az:8.3f}  1.00{b_factor:6.2f}           {element:>2s}"
            )
            atom_serial += 1

    # 1. Domain A beta-barrel coordinates centered around (-inter_port_distance_A/2 - 14.0, 0, 0)
    barrel_radius = 11.5
    n_strands = 11
    a_center_x = -0.5 * inter_port_distance_A - 14.0
    port_a_xyz = (-0.5 * inter_port_distance_A, 0.0, 0.0)
    for i, aa in enumerate(seq_a):
        if i >= len(seq_a) - 6:
            # C-terminal tail leading smoothly to Port A at (-inter_port_distance_A/2, 0, 0)
            frac = (i - (len(seq_a) - 6)) / 5.0
            x = (a_center_x + barrel_radius) * (1.0 - frac) + port_a_xyz[0] * frac
            y = 2.5 * math.sin(frac * math.pi)
            z = 1.8 * math.cos(frac * math.pi)
        else:
            strand = i % n_strands
            pos_in_strand = (i // n_strands) - 10
            theta = 2.0 * math.pi * strand / n_strands
            x = a_center_x + barrel_radius * math.cos(theta)
            y = barrel_radius * math.sin(theta)
            z = pos_in_strand * 1.9 * (1.0 if strand % 2 == 0 else -1.0)
        _emit_residue(aa, "A", res_serial, (x, y, z), b_factor=25.0)
        res_serial += 1

    linker_start_res = res_serial
    # 2. Linker coordinates along a smooth 3D catenary/helical spline from Port A to Port B
    port_b_xyz = (0.5 * inter_port_distance_A, 0.0, 0.0)
    L_link = len(seq_l)
    is_helical = "EAAAK" in seq_l
    for j, aa in enumerate(seq_l):
        u = (j + 1) / float(L_link + 1)
        x = port_a_xyz[0] * (1.0 - u) + port_b_xyz[0] * u
        arch = 5.5 * math.sin(math.pi * u)
        if is_helical:
            y = arch + 2.3 * math.cos(j * 1.745)
            z = 2.3 * math.sin(j * 1.745)
        else:
            y = arch + 1.2 * math.sin(j * 1.1)
            z = 1.6 * math.cos(j * 0.9)
        _emit_residue(aa, "A", res_serial, (x, y, z), b_factor=65.0)
        res_serial += 1
    linker_end_res = res_serial - 1

    # 3. Domain B beta-barrel coordinates centered around (+inter_port_distance_A/2 + 14.0, 0, 0)
    b_center_x = 0.5 * inter_port_distance_A + 14.0
    domain_b_start_res = res_serial
    for k, aa in enumerate(seq_b):
        if k < 6:
            frac = k / 5.0
            x = port_b_xyz[0] * (1.0 - frac) + (b_center_x - barrel_radius) * frac
            y = -2.2 * math.sin(frac * math.pi)
            z = 1.6 * math.sin(frac * math.pi)
        else:
            strand = k % n_strands
            pos_in_strand = (k // n_strands) - 10
            theta = 2.0 * math.pi * strand / n_strands
            x = b_center_x + barrel_radius * math.cos(theta)
            y = barrel_radius * math.sin(theta)
            z = pos_in_strand * 1.9 * (1.0 if strand % 2 == 0 else -1.0)
        _emit_residue(aa, "A", res_serial, (x, y, z), b_factor=45.0)
        res_serial += 1
    domain_b_end_res = res_serial - 1

    pdb_lines.append("TER")
    pdb_lines.append("END")
    pdb_text = "\n".join(pdb_lines) + "\n"

    # Ready-to-run PyMOL script (.pml)
    pml_lines = [
        "# LinkerForge PyMOL Visualization Script",
        f"# Construct: {domain_a_name} (1-{len(seq_a)}) — Linker ({linker_start_res}-{linker_end_res}) — {domain_b_name} ({domain_b_start_res}-{domain_b_end_res})",
        "bg_color white",
        "hide everything, all",
        "show cartoon, all",
        "set cartoon_fancy_helices, 1",
        "set cartoon_smooth_loops, 1",
        f"select domain_a, resi 1-{len(seq_a)}",
        f"select port_a, resi {max(1, len(seq_a)-2)}-{len(seq_a)}",
    ]
    if L_link > 0:
        pml_lines.extend(
            [
                f"select designed_linker, resi {linker_start_res}-{linker_end_res}",
                "color magenta, designed_linker",
                "show sticks, designed_linker and (name CA+C+N+O)",
            ]
        )
    pml_lines.extend(
        [
            f"select port_b, resi {domain_b_start_res}-{min(domain_b_end_res, domain_b_start_res+2)}",
            f"select domain_b, resi {domain_b_start_res}-{domain_b_end_res}",
            "color forest, domain_a",
            "color tv_red, domain_b",
            "color gold, port_a or port_b",
            "show spheres, (port_a or port_b) and name CA",
            "set sphere_scale, 0.65, (port_a or port_b) and name CA",
            "deselect",
            "orient",
            "zoom all, 4",
        ]
    )
    pml_text = "\n".join(pml_lines) + "\n"

    return {
        "domain_a_residue_range": [1, len(seq_a)],
        "port_a_residue_index": len(seq_a),
        "linker_residue_range": [linker_start_res, linker_end_res] if L_link > 0 else [],
        "port_b_residue_index": domain_b_start_res,
        "domain_b_residue_range": [domain_b_start_res, domain_b_end_res],
        "total_residues": domain_b_end_res,
        "total_atoms": atom_serial - 1,
        "pdb_text": pdb_text,
        "pymol_pml_script": pml_text,
    }


SFGFP_2B3P_SEQUENCE = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKR"
    "HDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNG"
    "IKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)

MCHERRY_2H5Q_SEQUENCE = (
    "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKH"
    "PADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDG"
    "ALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK"
)


def run_sfgfp_mcherry_fusion_demo(
    ranker: PairwiseCompatibilityRanker,
    authority_index: AuthoritySequenceIndex,
    *,
    requested_rigidity: str = "flexible",
    requested_cleavable: bool = False,
    protease_family: str = "TEV",
    target_length: int = 15,
    required_path_length_A: float = 28.4,
) -> dict[str, Any]:
    """Run the complete end-to-end `sfGFP (PDB 2B3P) — linker — mCherry (PDB 2H5Q)` fusion
    demonstration (`ARCHITECTURE_HANDOFF.md` §7), ranking LinkerForge candidates alongside all
    6 baselines and all 6 synthetic negative corruptions, and generating the 3D PDB + PyMOL script.
    """
    left_flank = SFGFP_2B3P_SEQUENCE[-30:]
    right_flank = MCHERRY_2H5Q_SEQUENCE[:30]
    port_a_sasa = 0.62  # sfGFP C-terminal tail (His231-Lys238) is solvent-exposed
    port_b_sasa = 0.58  # mCherry N-terminal tail (Met1-Glu6) is solvent-exposed

    x2_hits = authority_index.retrieve(
        target_length=target_length,
        rigidity=requested_rigidity,
        cleavable=requested_cleavable,
        protease_family=protease_family,
        left_flank_sequence=left_flank,
        right_flank_sequence=right_flank,
        top_k=6,
    )

    raw_candidates: list[tuple[str, str, str, str]] = []
    if requested_cleavable and requested_rigidity == "rigid":
        lib = load_library()
        comps = build_rigid_cleavable_composite(
            required_path_length_A,
            lib,
            protease=protease_family,
            min_length=max(8, target_length - 5),
            max_length=target_length + 10,
        )
        for idx, comp in enumerate(comps[:4]):
            raw_candidates.append(
                (f"lf_template_{idx+1}", "linkerforge_ranker", comp["sequence"], f"Part 2 Composite ({comp['candidate_id']})")
            )
    elif requested_cleavable:
        motif = CLEAVAGE_MOTIFS.get(protease_family, "ENLYFQS")
        flank_len = max(2, (target_length - len(motif)) // 2)
        pad = ("GGGGS" * ((flank_len // 5) + 2))[:flank_len]
        raw_candidates.append(
            ("lf_template_1", "linkerforge_ranker", f"{pad}{motif}{pad}", f"Part 2 Cleavable GS-{protease_family}-GS Template")
        )
        raw_candidates.append(
            ("lf_template_2", "linkerforge_ranker", f"GSAGSA{motif}GSAGSA", f"Part 2 Cleavable GSA-{protease_family}-GSA Template")
        )
    elif requested_rigidity == "rigid":
        raw_candidates.extend(
            [
                ("lf_template_1", "linkerforge_ranker", ("EAAAK" * ((target_length // 5) + 2))[:target_length], "Part 2 Template (Chen2013 Helical EAAAK)"),
                ("lf_template_2", "linkerforge_ranker", "AEAAAKEAAAKEAAAKA"[: max(12, target_length)], "Part 2 Template (A(EAAAK)3A Helical Rod)"),
                ("lf_template_3", "linkerforge_ranker", ("PAPAP" * ((target_length // 5) + 2))[:target_length], "Part 2 Template (Chen2013 Proline-Rich PAPAP)"),
            ]
        )
    else:
        raw_candidates.extend(
            [
                ("lf_template_1", "linkerforge_ranker", ("GGGGS" * ((target_length // 5) + 2))[:target_length], "Part 2 Template (Chen2013 Flexible GGGGS)"),
                ("lf_template_2", "linkerforge_ranker", ("GSAGSAAGSGEF" * ((target_length // 12) + 2))[:target_length], "Part 2 Template (Waldo GFP-Compatible GSAGSAAGSGEF)"),
                ("lf_template_3", "linkerforge_ranker", ("KESGSVSSEQLAQFRSLD" * ((target_length // 18) + 2))[:target_length], "Part 2 Template (Bird1988 Soluble Flexible Linker)"),
            ]
        )

    for idx, ac in enumerate(x2_hits):
        raw_candidates.append(
            (
                f"lf_authority_{idx+1}",
                "nearest_neighbour_retrieval",
                ac.sequence,
                f"X2 Authority ({ac.record_id}, {ac.evidence_tier})",
            )
        )

    # Add explicit baselines
    raw_candidates.extend(
        [
            ("baseline_direct_0aa", "direct_0aa_fusion", "", "Baseline 1: Direct 0-aa fusion"),
            ("baseline_fixed_ggggs3", "fixed_ggggs3", "GGGGSGGGGSGGGGS", "Baseline 2: Fixed (GGGGS)3"),
            (
                "baseline_random_gs",
                "random_length_matched_gs",
                "GSGGSGSGGSGSGGS"[:target_length],
                "Baseline 3: Random length-matched GS",
            ),
            (
                "baseline_geometry_only",
                "geometry_only",
                ("AAAGGG" * ((target_length // 6) + 2))[:target_length],
                "Baseline 4: Geometry-only length match",
            ),
        ]
    )

    # Add the 6 synthetic negative corruptions for transparent inspection
    top_pos_seq = raw_candidates[0][2] if raw_candidates else "GGGGSGGGGSGGGGS"
    synth_negs = generate_synthetic_negatives(
        top_pos_seq,
        target_length=target_length,
        requested_rigidity=requested_rigidity,
        left_flank=left_flank,
        right_flank=right_flank,
    )
    for c_type, n_seq in synth_negs.items():
        raw_candidates.append(
            (f"synthetic_neg_{c_type}", f"synthetic_negative:{c_type}", n_seq, f"Synthetic Negative ({c_type})")
        )

    scored_items: list[tuple[tuple[float, float], str, str, str, float, HardRejectionScorecard, dict[str, float], str]] = []
    seen_key: set[tuple[str, str]] = set()
    for cid, m_name, seq, prov in raw_candidates:
        if (m_name, seq) in seen_key:
            continue
        seen_key.add((m_name, seq))
        score, card, feats = ranker.score_candidate(
            seq,
            target_length=target_length,
            required_path_length_A=required_path_length_A,
            requested_rigidity=requested_rigidity,
            requested_cleavable=requested_cleavable,
            left_flank=left_flank,
            right_flank=right_flank,
            port_a_sasa=port_a_sasa,
            port_b_sasa=port_b_sasa,
        )
        sort_key = (1.0 if card.hard_rejection_pass else 0.0, score)
        scored_items.append((sort_key, cid, m_name, seq, score, card, feats, prov))

    scored_items.sort(key=lambda x: x[0], reverse=True)
    ranked_results: list[CandidateRankingResult] = []
    for r_idx, (_, cid, m_name, seq, score, card, feats, prov) in enumerate(scored_items, start=1):
        ranked_results.append(
            CandidateRankingResult(
                rank=r_idx,
                candidate_id=cid,
                method_name=m_name,
                sequence=seq,
                length=len(seq),
                compatibility_score=score,
                hard_rejection_pass=card.hard_rejection_pass,
                scorecard=card,
                feature_breakdown=feats,
                provenance=prov,
            )
        )

    top_candidate = ranked_results[0]
    structure_bundle = build_fusion_3d_pdb_and_pymol_script(
        domain_a_name="sfGFP_2B3P",
        domain_a_sequence=SFGFP_2B3P_SEQUENCE,
        linker_sequence=top_candidate.sequence,
        domain_b_name="mCherry_2H5Q",
        domain_b_sequence=MCHERRY_2H5Q_SEQUENCE,
        inter_port_distance_A=24.8,
    )

    return {
        "construct_name": "sfGFP (PDB 2B3P) — linker — mCherry (PDB 2H5Q)",
        "domain_a": {
            "name": "sfGFP",
            "pdb_id": "2B3P",
            "length": len(SFGFP_2B3P_SEQUENCE),
            "selected_port_residue": 238,
            "selected_port_motif": SFGFP_2B3P_SEQUENCE[-8:],
            "port_relative_sasa": port_a_sasa,
        },
        "domain_b": {
            "name": "mCherry",
            "pdb_id": "2H5Q",
            "length": len(MCHERRY_2H5Q_SEQUENCE),
            "selected_port_residue": 1,
            "selected_port_motif": MCHERRY_2H5Q_SEQUENCE[:8],
            "port_relative_sasa": port_b_sasa,
        },
        "junction_geometry": {
            "straight_ca_distance_A": 24.8,
            "required_path_length_A": required_path_length_A,
            "detour_ratio": round(required_path_length_A / 24.8, 3),
            "target_length_residues": target_length,
            "requested_rigidity": requested_rigidity,
            "requested_cleavable": requested_cleavable,
        },
        "top_selected_linker": top_candidate.to_dict(),
        "ranked_candidates": [r.to_dict() for r in ranked_results],
        "fusion_3d_structure": {
            "domain_a_residue_range": structure_bundle["domain_a_residue_range"],
            "linker_residue_range": structure_bundle["linker_residue_range"],
            "domain_b_residue_range": structure_bundle["domain_b_residue_range"],
            "total_residues": structure_bundle["total_residues"],
            "total_atoms": structure_bundle["total_atoms"],
            "pymol_pml_script": structure_bundle["pymol_pml_script"],
            "pdb_text": structure_bundle["pdb_text"],
        },
    }
