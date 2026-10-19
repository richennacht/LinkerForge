"""Deterministic Geometry-Plus-Template and Authority Linker Design Engine.

Implements Step 3 of `CURRENT_DIRECTION_2026-10-18.md` and `SPEC.md` section 7:
  - Q1 (Feasibility): Evaluates required residue windows (`flexible_window` or
    `helical_window`), requested-length bridgeability, Gaussian/worm-like chain
    strain, and chain topology (`end_to_end_tandem` vs `internal_insertion_two_tether`).
  - Q2 (Attachment & Flanking Context): Evaluates 3D anchor accessibility,
    orientation, and domain occlusion when `Chain` structures are supplied, and
    computes flanking sequence/domain context similarity against `mined_and_natural_linkers_v1`.
  - Q3 (Linker Sequence Selection): Retrieves and ranks candidates across cited
    literature templates (`linker_property_library.v1`), verified natural linkers
    (`natural_cath_pdb_verified`, `natural_disprot_curated_dfl`), and domain-gated
    high-confidence mined candidates (`mined_high_confidence_dfl_call`).
  - Enforces the two-part composite rule for `(rigidity="rigid", cleavability="cleavable")`:
    a scissile protease site requires local backbone unfolding (Hubbard et al. 1994;
    Chen et al. 2013) and cannot sit inside an unbroken alpha-helix, so it is
    constructed as a rigid helical/Pro-rich spacer plus a flexible GS-flanked
    protease cassette.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from math import sqrt
from pathlib import Path
from typing import Any, Sequence

from linkerforge.dataset_v1 import DEFAULT_DATASET_DIR, load_linker_records
from linkerforge.geometry.constants import (
    CA_VIRTUAL_BOND_ANGSTROM,
    DEFAULT_GS_COMPOSITION,
    HELIX_RISE_DEFAULT_ANGSTROM,
)
from linkerforge.geometry.junction import evaluate_junction, evaluate_tether_set
from linkerforge.geometry.span import (
    comfortable_residues,
    flexible_strain,
    helical_window,
    minimum_residues,
    required_length,
    verdict_for_requested_length,
)
from linkerforge.geometry.structure import Chain
from linkerforge.geometry.tether import implied_chain_topology
from linkerforge.properties import (
    call_cleavability,
    composition_features,
    fit_rigidity_model,
    load_library,
    modelling_subset,
)
from linkerforge.properties.library import Library

TOPOLOGIES = ("end_to_end_tandem", "internal_insertion_two_tether")
RIGIDITY_OPTIONS = ("flexible", "rigid")
CLEAVABILITY_OPTIONS = ("non_cleavable", "cleavable")

TIER_PRIORITY = {
    "literature_engineered_template": 0,
    "natural_cath_pdb_verified": 1,
    "natural_disprot_curated_dfl": 2,
    "mined_high_confidence_dfl_call": 3,
    "mined_provisional_dfl_candidate": 4,
}


@dataclass(frozen=True)
class DesignRequest:
    distance_angstrom: float
    rigidity: str = "flexible"
    cleavability: str = "non_cleavable"
    protease: str | None = None
    topology: str = "end_to_end_tandem"
    requested_residues: int | None = None
    return_distance_angstrom: float | None = None
    left_flank_sequence: str = ""
    right_flank_sequence: str = ""
    left_domain_id: str | None = None
    right_domain_id: str | None = None
    include_mined_candidates: bool = True
    min_length: int = 5
    max_length: int = 60


@dataclass(frozen=True)
class DesignCandidate:
    rank: int
    candidate_id: str
    authority_tier: str
    sequence: str
    length: int
    construction_mode: str
    geometric_verdict: str
    target_length_deviation_residues: float
    flanking_context_similarity: float
    rigidity_predicted_label: str
    rigidity_p_rigid: float
    cleavability_label: str
    matched_proteases: list[str]
    flags: list[str]
    provenance_note: str


@dataclass(frozen=True)
class DesignReport:
    schema_version: str
    q1_feasibility: dict[str, Any]
    q2_attachment_and_context: dict[str, Any]
    q3_candidates: list[dict[str, Any]]
    scope_and_integrity_notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _aa_composition_vector(seq: str) -> tuple[float, ...]:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    clean = [ch for ch in seq.upper() if ch in alphabet]
    if not clean:
        return tuple(0.0 for _ in alphabet)
    n = float(len(clean))
    return tuple(clean.count(aa) / n for aa in alphabet)


def _cosine_similarity(vec_a: Sequence[float], vec_b: Sequence[float]) -> float:
    dot = sum(a * b for a, b in zip(vec_a, vec_b))
    na = sqrt(sum(a * a for a in vec_a))
    nb = sqrt(sum(b * b for b in vec_b))
    if na <= 0.0 or nb <= 0.0:
        return 0.0
    return max(0.0, min(1.0, dot / (na * nb)))


def flanking_context_similarity(
    query_left_seq: str,
    query_right_seq: str,
    query_left_domain: str | None,
    query_right_domain: str | None,
    candidate_flank: dict[str, Any],
) -> float:
    """Compute [0, 1] flanking domain + boundary sequence similarity using flanks only."""
    cand_left_seq = str(candidate_flank.get("left_flank_sequence") or "")
    cand_right_seq = str(candidate_flank.get("right_flank_sequence") or "")
    cand_left_dom = candidate_flank.get("left_domain_id")
    cand_right_dom = candidate_flank.get("right_domain_id")

    scores: list[float] = []
    # Local 20-aa boundary window adjacent to the linker attachment ports
    if query_left_seq and cand_left_seq:
        scores.append(
            _cosine_similarity(
                _aa_composition_vector(query_left_seq[-25:]),
                _aa_composition_vector(cand_left_seq[-25:]),
            )
        )
    if query_right_seq and cand_right_seq:
        scores.append(
            _cosine_similarity(
                _aa_composition_vector(query_right_seq[:25]),
                _aa_composition_vector(cand_right_seq[:25]),
            )
        )
    if query_left_domain and cand_left_dom:
        scores.append(1.0 if query_left_domain == cand_left_dom else 0.0)
    if query_right_domain and cand_right_dom:
        scores.append(1.0 if query_right_domain == cand_right_dom else 0.0)
    if not scores:
        return 0.0
    return sum(scores) / len(scores)


def build_rigid_cleavable_composite(
    distance_angstrom: float,
    library: Library,
    *,
    protease: str | None = None,
    min_length: int = 12,
    max_length: int = 60,
) -> list[dict[str, Any]]:
    """Construct two-part (rigid spacer + flexible protease cassette) candidates.

    Follows `composite_templates[0]` (`rigid-cleavable-composite`) in
    `data/curated/linker_property_library_sources.v1.json` and `SPEC.md` section 7:
    protease cleavage requires local backbone unfolding (Hubbard et al. 1994), so
    a scissile site cannot sit inside an unbroken alpha-helix.
    """
    cleavable_entries = [
        e
        for e in library.entries
        if e.cleavability_label == "cleavable"
        and e.sequence is not None
        and len(e.sequence) >= 4
        and (
            protease is None
            or (e.protease and protease.lower() in e.protease.lower())
        )
    ]
    rigid_spacers = [
        ("rigid-eaaak1", "AEAAAKA"),
        ("rigid-eaaak2", "AEAAAKEAAAKA"),
        ("rigid-eaaak3", "AEAAAKEAAAKEAAAKA"),
        ("rigid-eaaak4", "AEAAAKEAAAKEAAAKEAAAKA"),
        ("rigid-papap", "PAPAP"),
        ("rigid-ap5", "APAPAPAPAP"),
    ]
    composites: list[dict[str, Any]] = []
    for entry in cleavable_entries:
        motif_seq = entry.sequence or ""
        # Ensure the scissile motif is flanked by flexible GS residues if not already GS-padded
        if motif_seq.startswith("GGGGS") and motif_seq.endswith("GGGGS"):
            cassette = motif_seq
        else:
            cassette = f"GGGGS{motif_seq}GGGGS"
        cassette_reach = len(cassette) * 2.5
        remaining_dist = max(0.0, distance_angstrom - cassette_reach)
        for spacer_id, spacer_seq in rigid_spacers:
            combined_seq = spacer_seq + cassette
            if not (min_length <= len(combined_seq) <= max_length):
                continue
            # Effective reach: helical spacer at 1.5 A/res + flexible cassette max span at 3.8 A/res
            min_helical_span = max(0.0, len(spacer_seq) * HELIX_RISE_DEFAULT_ANGSTROM - 7.6)
            max_combined_span = (
                len(spacer_seq) * HELIX_RISE_DEFAULT_ANGSTROM
                + (len(cassette) + 1) * CA_VIRTUAL_BOND_ANGSTROM
            )
            bridgeable = min_helical_span <= distance_angstrom <= max_combined_span
            target_spacer_len = max(5, round(remaining_dist / HELIX_RISE_DEFAULT_ANGSTROM))
            dev = abs(len(spacer_seq) - target_spacer_len)
            composites.append(
                {
                    "candidate_id": f"composite:{spacer_id}+{entry.entry_id}",
                    "authority_tier": "literature_engineered_template",
                    "sequence": combined_seq,
                    "length": len(combined_seq),
                    "construction_mode": "two_part_rigid_spacer_plus_flexible_protease_cassette",
                    "bridgeable": bridgeable,
                    "deviation": float(dev),
                    "protease": entry.protease or "unspecified",
                    "provenance_note": (
                        f"Two-part composite (`rigid-cleavable-composite`): rigid spacer "
                        f"{spacer_id} ({spacer_seq}) + flexible protease cassette from "
                        f"{entry.entry_id} ({library.citation_for(entry)}). Local unfolding "
                        "required for proteolysis (Hubbard et al. 1994)."
                    ),
                }
            )
    composites.sort(key=lambda c: (not c["bridgeable"], c["deviation"], c["length"]))
    return composites


def design_linker_candidates(
    request: DesignRequest,
    *,
    dataset_dir: Path = DEFAULT_DATASET_DIR,
    library: Library | None = None,
    domain_a: Chain | None = None,
    anchor_a: int | None = None,
    domain_b: Chain | None = None,
    anchor_b: int | None = None,
    top_k: int = 5,
) -> DesignReport:
    """Run Q1 (feasibility), Q2 (attachment/context), and Q3 (candidate selection)."""
    if request.rigidity not in RIGIDITY_OPTIONS:
        raise ValueError(f"rigidity must be one of {RIGIDITY_OPTIONS}")
    if request.cleavability not in CLEAVABILITY_OPTIONS:
        raise ValueError(f"cleavability must be one of {CLEAVABILITY_OPTIONS}")
    if request.topology not in TOPOLOGIES:
        raise ValueError(f"topology must be one of {TOPOLOGIES}")
    if request.distance_angstrom < 0:
        raise ValueError("distance_angstrom must be non-negative")

    lib = library if library is not None else load_library()
    rigidity_entries = modelling_subset(lib, "rigidity")
    rigidity_model = fit_rigidity_model(
        [e.sequence for e in rigidity_entries if e.sequence is not None],
        [e.rigidity_label for e in rigidity_entries if e.sequence is not None],
    )

    # -------------------------------------------------------------------------
    # Q1: Geometric Feasibility & Topology
    # -------------------------------------------------------------------------
    geom_linker_type = "helical" if request.rigidity == "rigid" else "flexible"
    window = required_length(request.distance_angstrom, linker_type=geom_linker_type)
    if request.rigidity == "flexible":
        target_length = comfortable_residues(request.distance_angstrom, DEFAULT_GS_COMPOSITION) or max(
            request.min_length, minimum_residues(request.distance_angstrom) + 2
        )
    else:
        min_r = int(window["minimum_residues"])
        max_r = window["maximum_residues"]
        target_length = (
            round(0.5 * (min_r + int(max_r)))
            if max_r is not None
            else max(request.min_length, min_r)
        )
    target_length = max(request.min_length, min(request.max_length, int(target_length)))
    if request.requested_residues is not None:
        target_length = int(request.requested_residues)

    req_verdict = verdict_for_requested_length(window, target_length)
    strain_info = (
        flexible_strain(request.distance_angstrom, target_length)
        if request.rigidity == "flexible"
        else None
    )

    if request.topology == "end_to_end_tandem":
        topo_info = implied_chain_topology(1, peptide_tethers=1)
        return_window = None
    else:
        topo_info = implied_chain_topology(2, peptide_tethers=2)
        ret_dist = (
            request.return_distance_angstrom
            if request.return_distance_angstrom is not None
            else request.distance_angstrom
        )
        return_window = required_length(ret_dist, linker_type=geom_linker_type)

    q1_feasibility = {
        "topology": request.topology,
        "chain_topology": topo_info,
        "requested_rigidity": request.rigidity,
        "requested_cleavability": request.cleavability,
        "primary_span_window": window,
        "target_residues_used": target_length,
        "target_length_verdict": req_verdict,
        "flexible_strain": strain_info,
        "return_tether_span_window": return_window,
    }

    # -------------------------------------------------------------------------
    # Q2: Attachment Site & Flanking Context Assessment
    # -------------------------------------------------------------------------
    junction_3d = None
    if (
        domain_a is not None
        and anchor_a is not None
        and domain_b is not None
        and anchor_b is not None
    ):
        junction_3d = evaluate_junction(
            domain_a,
            anchor_a,
            domain_b,
            anchor_b,
            linker_type=geom_linker_type,
            requested_residues=target_length,
        )

    q2_attachment = {
        "structure_junction_report": junction_3d,
        "left_domain_id": request.left_domain_id,
        "right_domain_id": request.right_domain_id,
        "left_flank_length": len(request.left_flank_sequence),
        "right_flank_length": len(request.right_flank_sequence),
        "attachment_scope_note": (
            "End-to-end tandem fusion requires 1 polypeptide linker; internal insertion "
            "requires 2 polypeptide tethers leaving and re-entering the host domain."
        ),
    }

    # -------------------------------------------------------------------------
    # Q3: Candidate Linker Selection (Templates + Mined/Natural Authorities)
    # -------------------------------------------------------------------------
    scored_candidates: list[DesignCandidate] = []

    if request.rigidity == "rigid" and request.cleavability == "cleavable":
        comp_list = build_rigid_cleavable_composite(
            request.distance_angstrom,
            lib,
            protease=request.protease,
            min_length=request.min_length,
            max_length=request.max_length,
        )
        for idx, item in enumerate(comp_list[:top_k], start=1):
            seq = item["sequence"]
            cleav_call = call_cleavability(seq, lib)
            p_rigid = float(rigidity_model.predict_proba(seq))
            flags = [
                "composite_rigid_spacer_plus_flexible_cassette",
                "requires_local_unfolding_at_cassette",
                "bridgeable" if item["bridgeable"] else "span_mismatch",
            ]
            if request.topology == "internal_insertion_two_tether":
                flags.append("requires_second_return_tether_for_internal_insertion")
            scored_candidates.append(
                DesignCandidate(
                    rank=idx,
                    candidate_id=item["candidate_id"],
                    authority_tier=item["authority_tier"],
                    sequence=seq,
                    length=item["length"],
                    construction_mode=item["construction_mode"],
                    geometric_verdict="bridgeable" if item["bridgeable"] else "outside_window",
                    target_length_deviation_residues=round(item["deviation"], 3),
                    flanking_context_similarity=0.0,
                    rigidity_predicted_label="rigid_spacer_plus_flexible_cassette",
                    rigidity_p_rigid=round(p_rigid, 4),
                    cleavability_label=cleav_call.label,
                    matched_proteases=cleav_call.proteases,
                    flags=flags,
                    provenance_note=item["provenance_note"],
                )
            )
    else:
        allowed_tiers = [
            "literature_engineered_template",
            "natural_cath_pdb_verified",
            "natural_disprot_curated_dfl",
        ]
        if request.include_mined_candidates:
            allowed_tiers.append("mined_high_confidence_dfl_call")
        records = load_linker_records(
            dataset_dir,
            tiers=allowed_tiers,
            min_length=request.min_length,
            max_length=request.max_length,
        )

        pool: list[tuple[tuple[Any, ...], DesignCandidate]] = []
        for row in records:
            seq = row["linker_sequence"]
            length = int(row["linker_length"])
            props = row["biophysical_properties"]
            pred_rig = props["rigidity_predicted_label"]
            curated_rig = props.get("curated_rigidity_label")
            eff_rig = (
                curated_rig
                if curated_rig in ("flexible", "rigid")
                else pred_rig
            )
            cleav_lbl = props["cleavability_label"]
            proteases = list(props.get("matched_proteases") or [])

            # Filter on requested rigidity and cleavability axes
            if request.rigidity != eff_rig:
                continue
            if request.cleavability == "cleavable":
                if cleav_lbl != "cleavable":
                    continue
                if request.protease and not any(
                    request.protease.lower() in p.lower() for p in proteases
                ):
                    continue
            else:
                if cleav_lbl == "cleavable":
                    continue

            verdict_info = verdict_for_requested_length(window, length)
            verdict_str = str(verdict_info["verdict"])
            if verdict_str != "bridgeable":
                continue

            flank_sim = flanking_context_similarity(
                request.left_flank_sequence,
                request.right_flank_sequence,
                request.left_domain_id,
                request.right_domain_id,
                row.get("flanking_context") or {},
            )
            len_dev = abs(length - target_length)
            flags = [verdict_str, f"tier:{row['authority_tier']}"]
            if request.rigidity == "flexible":
                strain = flexible_strain(request.distance_angstrom, length)
                if strain["near_full_extension"]:
                    flags.append("near_full_extension_strain")
                if not strain["requires_stretching_beyond_rms"]:
                    flags.append("unstrained_rms_reach")
            if request.topology == "internal_insertion_two_tether":
                flags.append("requires_second_return_tether_for_internal_insertion")

            cand = DesignCandidate(
                rank=0,
                candidate_id=row["record_id"],
                authority_tier=row["authority_tier"],
                sequence=seq,
                length=length,
                construction_mode="single_segment_template_or_authority_retrieval",
                geometric_verdict=verdict_str,
                target_length_deviation_residues=float(len_dev),
                flanking_context_similarity=round(flank_sim, 4),
                rigidity_predicted_label=eff_rig,
                rigidity_p_rigid=round(float(props["rigidity_p_rigid"]), 4),
                cleavability_label=cleav_lbl,
                matched_proteases=proteases,
                flags=flags,
                provenance_note=str(
                    props.get("citation")
                    or f"{row['authority_tier']} ({row['evidence_grade']})"
                ),
            )
            # Sort key: prefer unstrained reach, high flanking similarity (when query flanks provided),
            # small deviation from target length, and higher authority tier
            has_query_flank = bool(
                request.left_flank_sequence
                or request.right_flank_sequence
                or request.left_domain_id
                or request.right_domain_id
            )
            sort_key = (
                "near_full_extension_strain" in flags,
                -round(flank_sim, 3) if has_query_flank else 0.0,
                len_dev,
                TIER_PRIORITY.get(row["authority_tier"], 99),
                row["record_id"],
            )
            pool.append((sort_key, cand))

        pool.sort(key=lambda item: item[0])
        seen_seqs: set[str] = set()
        for _, cand in pool:
            if cand.sequence in seen_seqs:
                continue
            seen_seqs.add(cand.sequence)
            scored_candidates.append(
                DesignCandidate(
                    rank=len(scored_candidates) + 1,
                    candidate_id=cand.candidate_id,
                    authority_tier=cand.authority_tier,
                    sequence=cand.sequence,
                    length=cand.length,
                    construction_mode=cand.construction_mode,
                    geometric_verdict=cand.geometric_verdict,
                    target_length_deviation_residues=cand.target_length_deviation_residues,
                    flanking_context_similarity=cand.flanking_context_similarity,
                    rigidity_predicted_label=cand.rigidity_predicted_label,
                    rigidity_p_rigid=cand.rigidity_p_rigid,
                    cleavability_label=cand.cleavability_label,
                    matched_proteases=cand.matched_proteases,
                    flags=cand.flags,
                    provenance_note=cand.provenance_note,
                )
            )
            if len(scored_candidates) >= top_k:
                break

    return DesignReport(
        schema_version="1.0.0",
        q1_feasibility=q1_feasibility,
        q2_attachment_and_context=q2_attachment,
        q3_candidates=[asdict(c) for c in scored_candidates],
        scope_and_integrity_notes=[
            "Geometric bridgeability is a necessary physical condition (3.8 A/res C-alpha ceiling; 1.5 A/res helical rise), never a fusion expression or activity probability (SPEC.md section 5.8).",
            "Cleavable + rigid requests return a two-part composite construct (rigid spacer + flexible GS-flanked protease cassette) because proteolysis requires local backbone unfolding (Hubbard et al. 1994).",
            "Mined candidates (`mined_high_confidence_dfl_call`) are model-called inter-domain spans from Part 1 and carry explicit provenance flags distinct from curated natural or literature authorities.",
        ],
    )
