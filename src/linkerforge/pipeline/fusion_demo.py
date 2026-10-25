"""Fusion demonstration on real coordinates: sfGFP (PDB 2B3P) - linker - mCherry (PDB 2H5Q).

What this module does
---------------------
* reads the two experimental structures from `data/raw/pdb/` (downloaded from RCSB);
* takes the last / first *resolved* residues as attachment ports (the deposited constructs leave
  the terminal tails unresolved; their lengths are reported and are NOT modelled);
* places domain B by a seeded rigid-body search that keeps all heavy atoms >= 4 A from domain A and
  minimises the port-to-port distance.  This is ONE geometrically admissible arrangement chosen by
  a script; it is not a prediction of the fusion's true conformation;
* runs the repository's `evaluate_junction` on the placed pair for the real span / detour / exit
  geometry, proposes linkers (curated library, repeat templates, X2 retrieval), scores them with the
  scorecard and ranker, and writes a three-chain PDB (A = sfGFP, B = mCherry, L = linker C-alpha
  trace) plus a PyMOL script.

Peptide bonds at the junctions are not modelled; the linker trace is a schematic path.
"""

from __future__ import annotations

import math
import random
from pathlib import Path
from typing import Any

import numpy as np

from linkerforge.geometry.constants import DEFAULT_GS_COMPOSITION
from linkerforge.geometry.junction import evaluate_junction
from linkerforge.geometry.sasa import relative_accessibility, sphere_points
from linkerforge.geometry.span import required_length, verdict_for_requested_length
from linkerforge.geometry.structure import read_chains

from .evaluation import geometry_length
from .ranker import (
    AuthoritySequenceIndex,
    PairwiseCompatibilityRanker,
    evaluate_hard_rejection_scorecard,
    generate_synthetic_negatives,
)

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H",
    "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q",
    "ARG": "R", "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y", "MSE": "M",
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items() if k != "MSE"}
MIN_INTERDOMAIN_GAP_A = 4.0


def read_pdb_chain_lines(path: Path, chain: str) -> list[dict[str, Any]]:
    """First-model ATOM/HETATM lines of one chain (waters excluded; altloc blank or 'A')."""
    rows: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        if line[21] != chain or line[17:20].strip() == "HOH" or line[16] not in (" ", "A"):
            continue
        rows.append(
            {
                "line": line,
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "resseq": int(line[22:26]),
                "icode": line[26],
                "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                "element": (line[76:78].strip() or line[12:16].strip()[0]).upper(),
            }
        )
    return rows


def unresolved_residues(path: Path, chain: str) -> list[int]:
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("REMARK 465") and len(line) > 25:
            parts = line.split()
            if len(parts) == 5 and parts[3] == chain and parts[4].lstrip("-").isdigit():
                out.append(int(parts[4]))
    return sorted(out)


def _sequence_and_ca(rows: list[dict[str, Any]]) -> tuple[list[tuple[int, str, np.ndarray]], str]:
    cas = []
    seen = set()
    for r in rows:
        key = (r["resseq"], r["icode"])
        if r["name"] == "CA" and key not in seen:
            seen.add(key)
            cas.append((r["resseq"], THREE_TO_ONE.get(r["resname"], "X"), r["xyz"]))
    cas.sort(key=lambda x: x[0])
    return cas, "".join(c[1] for c in cas)


def _min_dist(a: np.ndarray, b: np.ndarray) -> float:
    d2 = (a * a).sum(1)[:, None] + (b * b).sum(1)[None, :] - 2.0 * a @ b.T
    return float(math.sqrt(max(float(d2.min()), 0.0)))


def _random_rotation(rng: random.Random) -> np.ndarray:
    m = np.array([[rng.gauss(0, 1) for _ in range(3)] for _ in range(3)])
    q, r = np.linalg.qr(m)
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def place_domain_b(
    a_xyz: np.ndarray,
    b_xyz: np.ndarray,
    port_a: np.ndarray,
    port_b: np.ndarray,
    *,
    n_rotations: int = 400,
    seed: int = 7,
) -> dict[str, Any]:
    """Seeded rigid-body search: clash-free placement minimising the port-to-port distance."""
    rng = random.Random(seed)
    c_a = a_xyz.mean(0)
    c_b = b_xyz.mean(0)
    u = port_a - c_a
    u = u / np.linalg.norm(u)
    best: dict[str, Any] | None = None
    for _ in range(n_rotations):
        rot = _random_rotation(rng)
        b_rot = (b_xyz - c_b) @ rot.T
        pb_rot = (port_b - c_b) @ rot.T

        def free_at(s: float) -> bool:
            return _min_dist(a_xyz, b_rot + port_a + s * u) >= MIN_INTERDOMAIN_GAP_A

        lo, hi = 0.0, 160.0
        if not free_at(hi):
            continue
        for _ in range(14):
            mid = 0.5 * (lo + hi)
            if free_at(mid):
                hi = mid
            else:
                lo = mid
        s = hi
        pb_placed = pb_rot + port_a + s * u
        dd = float(np.linalg.norm(pb_placed - port_a))
        if best is None or dd < best["port_distance_A"]:
            best = {"rot": rot, "s": s, "port_distance_A": dd}
    if best is None:
        raise RuntimeError("no clash-free placement found")
    rot, s = best["rot"], best["s"]
    shift = port_a + s * u

    def transform(x: np.ndarray) -> np.ndarray:
        return (x - c_b) @ rot.T + shift

    placed = transform(b_xyz)
    best.update(
        transform=transform,
        min_interatomic_distance_A=round(_min_dist(a_xyz, placed), 3),
        n_rotations_tried=n_rotations,
        seed=seed,
        translation_along_outward_axis_A=round(float(s), 3),
    )
    return best


def _fmt_line(line: str, chain: str, xyz: np.ndarray | None = None) -> str:
    out = line[:21] + chain + line[22:]
    if xyz is not None:
        out = out[:30] + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}" + out[54:]
    return out


def linker_ca_path(
    pa: np.ndarray, pb: np.ndarray, n: int, outward: np.ndarray, bond: float = 3.8
) -> tuple[np.ndarray, dict[str, Any]]:
    """Schematic C-alpha trace of an n-residue linker between two port C-alphas."""
    d = float(np.linalg.norm(pb - pa))
    axis = (pb - pa) / max(d, 1e-9)
    w = outward - axis * float(outward @ axis)
    if np.linalg.norm(w) < 1e-6:
        w = np.cross(axis, np.array([1.0, 0.0, 0.0]))
        if np.linalg.norm(w) < 1e-6:
            w = np.cross(axis, np.array([0.0, 1.0, 0.0]))
    w = w / np.linalg.norm(w)
    max_reach = (n + 1) * bond
    reachable = d <= max_reach
    ts = np.linspace(0.0, 1.0, 400)
    target_arc = min((n + 1) * 3.6, max_reach * 0.98)

    def curve(h: float, tt: np.ndarray) -> np.ndarray:
        return pa[None, :] + np.outer(tt, pb - pa) + h * np.outer(np.sin(np.pi * tt), w)

    def arc(h: float) -> float:
        pts = curve(h, ts)
        return float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())

    h = 0.0
    if reachable and arc(0.0) < target_arc:
        lo, hi = 0.0, max_reach
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if arc(mid) < target_arc:
                lo = mid
            else:
                hi = mid
        h = 0.5 * (lo + hi)
    if n > 0:
        grid = curve(h, ts)
        cum = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(grid, axis=0), axis=1))])
        wanted = cum[-1] * np.array([i / (n + 1) for i in range(1, n + 1)])
        t_equal = np.interp(wanted, cum, ts)
        pts = curve(h, t_equal)
    else:
        pts = np.zeros((0, 3))
    spacing = np.linalg.norm(np.diff(np.vstack([pa, pts, pb]), axis=0), axis=1) if n > 0 else np.array([d])
    return pts, {
        "reachable_at_3p8A_per_bond": reachable,
        "bulge_amplitude_A": round(h, 3),
        "ca_ca_spacing_min_A": round(float(spacing.min()), 3),
        "ca_ca_spacing_max_A": round(float(spacing.max()), 3),
    }


def build_assembly(
    a_rows: list[dict[str, Any]],
    b_rows: list[dict[str, Any]],
    b_placed: np.ndarray,
    linker_seq: str,
    linker_xyz: np.ndarray,
    *,
    header_remarks: list[str],
) -> tuple[str, dict[str, Any]]:
    lines = ["HEADER    LINKERFORGE SCHEMATIC ASSEMBLY (NOT A PREDICTED STRUCTURE)"]
    lines += [f"REMARK 999 {r}" for r in header_remarks]
    serial = 1
    for rows, chain, xyz in ((a_rows, "A", None), (b_rows, "B", b_placed)):
        for i, r in enumerate(rows):
            ln = _fmt_line(r["line"], chain, None if xyz is None else xyz[i])
            lines.append(ln[:6] + f"{serial:5d}" + ln[11:])
            serial += 1
        lines.append("TER")
    for j, aa in enumerate(linker_seq):
        x, y, z = linker_xyz[j]
        lines.append(
            f"ATOM  {serial:5d}  CA  {ONE_TO_THREE.get(aa, 'GLY')} L{j + 1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
        serial += 1
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n", {"total_atoms_written": serial - 1}


def pymol_script(pdb_name: str, linker_len: int, port_a_resi: int, port_b_resi: int) -> str:
    return "\n".join(
        [
            "# LinkerForge schematic assembly -- NOT a predicted structure.",
            f"load {pdb_name}",
            "bg_color white",
            "hide everything",
            "show cartoon, chain A or chain B",
            "color forest, chain A",
            "color tv_red, chain B",
            *(
                [
                    "show spheres, chain L",
                    "set sphere_scale, 0.45, chain L",
                    "color magenta, chain L",
                ]
                if linker_len > 0
                else []
            ),
            f"show spheres, chain A and resi {port_a_resi} and name CA",
            f"show spheres, chain B and resi {port_b_resi} and name CA",
            f"color gold, (chain A and resi {port_a_resi} and name CA) or (chain B and resi {port_b_resi} and name CA)",
            "orient",
            "",
        ]
    )


def _repeat(unit: str, n: int) -> str:
    return (unit * (n // len(unit) + 1))[:n] if n > 0 else ""


def run_fusion_demo(
    ranker: PairwiseCompatibilityRanker,
    authority_index: AuthoritySequenceIndex,
    library: Any,
    *,
    pdb_dir: Path,
    out_dir: Path,
    seed: int = 7,
) -> dict[str, Any]:
    pa_path, pb_path = pdb_dir / "2B3P.pdb", pdb_dir / "2H5Q.pdb"
    a_rows = read_pdb_chain_lines(pa_path, "A")
    b_rows = read_pdb_chain_lines(pb_path, "A")
    a_cas, a_seq = _sequence_and_ca(a_rows)
    b_cas, b_seq = _sequence_and_ca(b_rows)
    port_a_res, port_b_res = a_cas[-1], b_cas[0]
    unres_a = [r for r in unresolved_residues(pa_path, "A") if r > port_a_res[0]]
    unres_b = [r for r in unresolved_residues(pb_path, "A") if r < port_b_res[0]]

    a_xyz = np.array([r["xyz"] for r in a_rows])
    b_xyz = np.array([r["xyz"] for r in b_rows])
    placement = place_domain_b(a_xyz, b_xyz, port_a_res[2], port_b_res[2], seed=seed)
    b_placed = placement["transform"](b_xyz)
    pb_placed = placement["transform"](port_b_res[2])

    # Real junction geometry from the repository calculator on the placed pair.
    out_dir.mkdir(parents=True, exist_ok=True)
    pair_path = out_dir / "sfgfp_mcherry_pair_for_geometry.pdb"
    pair_lines = []
    for rows, chain, xyz in ((a_rows, "A", None), (b_rows, "B", b_placed)):
        for i, r in enumerate(rows):
            if r["line"].startswith("ATOM  "):
                pair_lines.append(_fmt_line(r["line"], chain, None if xyz is None else xyz[i]))
        pair_lines.append("TER")
    pair_path.write_text("\n".join(pair_lines + ["END"]) + "\n", encoding="utf-8")
    chains = read_chains(pair_path)
    junction = evaluate_junction(
        chains["A"], port_a_res[0], chains["B"], port_b_res[0], linker_type="flexible"
    )
    pair_path.unlink()
    path_A = float(junction["span"]["path_length_used_angstrom"])
    straight_A = float(junction["span"]["straight_ca_distance_angstrom"])

    # Port accessibility in the isolated domains (the usual per-domain definition).
    def port_rsa(path: Path, chain_key: str, number: int) -> float | None:
        ch = read_chains(path)[chain_key]
        atoms = tuple(a for r in ch.residues for a in r.atoms)
        res = ch.residue(number)
        return relative_accessibility(res, atoms, points=sphere_points(256)).get(
            "relative_solvent_accessibility"
        )

    rsa_a = port_rsa(pa_path, "A", port_a_res[0])
    rsa_b = port_rsa(pb_path, "A", port_b_res[0])

    left_flank, right_flank = a_seq[-30:], b_seq[:30]
    flex_win = required_length(path_A, linker_type="flexible")
    helix_win = required_length(path_A, linker_type="helical")
    l_flex = geometry_length(path_A, "unstrained_reach")
    h_min, h_max = int(helix_win["minimum_residues"]), helix_win.get("maximum_residues")
    l_rigid = int(round((h_min + int(h_max)) / 2)) if h_max is not None else h_min + 4

    entries = [e for e in library.entries if e.sequence and not e.excluded_from_model]

    def lib_candidates(pred, lo, hi):
        return [
            (f"library:{e.entry_id}", e.sequence, library.citation_for(e))
            for e in entries
            if pred(e) and lo <= len(e.sequence) <= hi
        ]

    modes: dict[str, dict[str, Any]] = {}
    mode_specs = {
        "flexible_non_cleavable": dict(
            rigidity="flexible", cleavable=False, target=l_flex, window=flex_win,
            cands=lib_candidates(lambda e: e.rigidity_label == "flexible", l_flex - 6, l_flex + 6)
            + [(f"template:(GGGGS)n_to_{l_flex}aa", _repeat("GGGGS", l_flex), "repeat of library entry flex-ggggs1 to target length")],
        ),
        "rigid_non_cleavable": dict(
            rigidity="rigid", cleavable=False, target=l_rigid, window=helix_win,
            cands=lib_candidates(lambda e: e.rigidity_label == "rigid", max(4, l_rigid - 8), l_rigid + 8)
            + [(f"template:A(EAAAK)nA_to_{l_rigid}aa", ("A" + _repeat("EAAAK", max(0, l_rigid - 2)) + "A")[:l_rigid], "A(EAAAK)nA pattern (library entries rigid-eaaak*, rigid-aeaaakeaaaka) truncated to target length")],
        ),
        "flexible_cleavable_factor_xa": dict(
            rigidity="flexible", cleavable=True, target=l_flex, window=flex_win,
            cands=[
                (
                    f"constructed:GGGGS-pad+library:cleav-factor-xa_to_{l_flex}aa",
                    (lambda pad: _repeat("GGGGS", pad) + "GGIEGRGS" + _repeat("GGGGS", l_flex - 8 - pad))(max(0, (l_flex - 8) // 2)),
                    "Factor Xa site GGIEGRGS from library entry cleav-factor-xa, flanked by GGGGS padding (padding is constructed, not from the library)",
                )
            ],
        ),
    }
    rng = random.Random(seed)
    for mode, spec in mode_specs.items():
        hits = authority_index.retrieve(
            target_length=spec["target"], rigidity=spec["rigidity"], cleavable=spec["cleavable"],
            left_flank_sequence=left_flank, right_flank_sequence=right_flank, top_k=4,
        )
        cands = list(spec["cands"]) + [
            (f"x2:{h.record_id}", h.sequence, f"X2 retrieval from {h.evidence_tier}") for h in hits
        ]
        cands += [
            ("baseline:direct_0aa_fusion", "", "baseline"),
            ("baseline:fixed_(GGGGS)3", "GGGGSGGGGSGGGGS", "baseline"),
            ("baseline:random_GS_target_length", "".join(rng.choice("GS") for _ in range(spec["target"])), "baseline (seeded random)"),
        ]
        negs = generate_synthetic_negatives(
            cands[0][1] or "GGGGS", target_length=spec["target"], requested_rigidity=spec["rigidity"],
            left_flank=left_flank, right_flank=right_flank,
        )
        cands += [(f"corruption:{k}", v, "constructed negative (by construction)") for k, v in negs.items()]

        rows_out = []
        seen = set()
        for cid, seq, prov in cands:
            if (cid.split(":")[0], seq) in seen:
                continue
            seen.add((cid.split(":")[0], seq))
            sc, card, feats = ranker.score_candidate(
                seq, target_length=spec["target"], required_path_length_A=path_A,
                requested_rigidity=spec["rigidity"], requested_cleavable=spec["cleavable"],
                left_flank=left_flank, right_flank=right_flank, port_a_sasa=rsa_a, port_b_sasa=rsa_b,
            )
            verdict = verdict_for_requested_length(spec["window"], len(seq))["verdict"] if seq else "no_linker"
            rows_out.append(
                {
                    "candidate_id": cid, "sequence": seq, "length": len(seq), "provenance": prov,
                    "ranker_score": round(sc, 4), "hard_rejection_pass": card.hard_rejection_pass,
                    "rejection_reasons": list(card.rejection_reasons),
                    "geometric_length_verdict": verdict,
                }
            )
        rows_out.sort(key=lambda r: (not r["hard_rejection_pass"], -r["ranker_score"]))
        for i, r in enumerate(rows_out, 1):
            r["rank"] = i
        top = next((r for r in rows_out if r["hard_rejection_pass"] and not r["candidate_id"].startswith(("baseline:", "corruption:"))), rows_out[0])

        outward = ((a_xyz.mean(0) + b_placed.mean(0)) / 2.0) * -1.0 + (port_a_res[2] + pb_placed) / 2.0
        best_pts, best_info, best_clash = None, None, None
        for sign in (1.0, -1.0):
            pts, info = linker_ca_path(port_a_res[2], pb_placed, len(top["sequence"]), sign * outward)
            allxyz = np.vstack([a_xyz, b_placed])
            clashes = int(sum(1 for p in pts if _min_dist(p[None, :], allxyz) < 3.0))
            if best_clash is None or clashes < best_clash:
                best_pts, best_info, best_clash = pts, dict(info, linker_ca_within_3A_of_domain_atoms=clashes), clashes
        stem = f"sfgfp_mcherry_{mode}"
        pdb_text, stats = build_assembly(
            a_rows, b_rows, b_placed, top["sequence"], best_pts,
            header_remarks=[
                "Domains: real coordinates from PDB 2B3P (chain A) and 2H5Q (chain A).",
                f"Domain B placed by seeded rigid-body search (seed {seed}); min interdomain heavy-atom distance {placement['min_interatomic_distance_A']} A.",
                "Linker is a schematic C-alpha trace; peptide bonds at the junctions are not modelled.",
                f"Unresolved terminal tails not modelled: sfGFP C-terminus {len(unres_a)} residues, mCherry N-terminus {len(unres_b)} residues.",
            ],
        )
        (out_dir / f"{stem}.pdb").write_text(pdb_text, encoding="utf-8")
        (out_dir / f"{stem}.pml").write_text(
            pymol_script(f"{stem}.pdb", len(top["sequence"]), port_a_res[0], port_b_res[0]), encoding="utf-8"
        )
        ca_trace: dict[str, list[list[float]]] = {"A": [], "B": [], "L": []}
        for ln in pdb_text.splitlines():
            if ln.startswith(("ATOM  ", "HETATM")) and ln[12:16].strip() == "CA" and ln[21] in ca_trace:
                ca_trace[ln[21]].append([round(float(ln[30:38]), 2), round(float(ln[38:46]), 2), round(float(ln[46:54]), 2)])
        modes[mode] = {
            "requested_rigidity": spec["rigidity"], "requested_cleavable": spec["cleavable"],
            "target_length": spec["target"], "ranked_candidates": rows_out,
            "top_non_baseline_candidate": top["candidate_id"], "linker_trace": best_info,
            "structure_files": [f"reports/structures/{stem}.pdb", f"reports/structures/{stem}.pml"],
            "ca_trace": ca_trace,
            **stats,
        }

    return {
        "status": "schematic_assembly_on_real_domain_coordinates_not_a_structure_prediction",
        "context": {
            "path_A": round(path_A, 3), "rsa_a": rsa_a, "rsa_b": rsa_b,
            "left_flank": left_flank, "right_flank": right_flank,
        },
        "domain_a": {"pdb": "2B3P", "chain": "A", "resolved_residues": len(a_cas),
                     "port_residue_number": port_a_res[0], "port_residue": port_a_res[1],
                     "port_relative_sasa_isolated_domain": None if rsa_a is None else round(rsa_a, 4),
                     "unresolved_c_terminal_tail_residues": len(unres_a)},
        "domain_b": {"pdb": "2H5Q", "chain": "A", "resolved_residues": len(b_cas),
                     "port_residue_number": port_b_res[0], "port_residue": port_b_res[1],
                     "port_relative_sasa_isolated_domain": None if rsa_b is None else round(rsa_b, 4),
                     "unresolved_n_terminal_tail_residues": len(unres_b)},
        "placement": {
            "method": "seeded random rotations, clash-free along the outward port axis, minimum port distance",
            "n_rotations_tried": placement["n_rotations_tried"], "seed": placement["seed"],
            "min_interdomain_heavy_atom_distance_A": placement["min_interatomic_distance_A"],
            "required_gap_A": MIN_INTERDOMAIN_GAP_A,
            "straight_port_ca_distance_A": round(straight_A, 3),
            "detour_path_length_A": round(path_A, 3),
        },
        "geometry_windows": {
            "flexible": {k: v for k, v in flex_win.items() if not k.endswith("basis") and k != "chain_model"},
            "helical": {k: v for k, v in helix_win.items() if not k.endswith("basis")},
        },
        "modes": modes,
    }
