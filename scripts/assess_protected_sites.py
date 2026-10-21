#!/usr/bin/env python3
"""Map protected residues and calculate solvent/contact evidence from a PDB.

Run with a Biopython environment.  This reports measurements without applying
an acceptance threshold; the project must freeze biological thresholds before
using them as hard rejection rules.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess

import Bio
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley


MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_sites(text: str, domain_a_length: int, linker_length: int, domain_b_length: int) -> list[dict]:
    sites: list[dict] = []
    for raw in text.replace(",", ";").split(";"):
        token = raw.strip()
        if not token:
            continue
        try:
            domain, span = [part.strip().upper() for part in token.split(":", 1)]
            start_text, *end_text = span.split("-", 1)
            start, end = int(start_text), int(end_text[0]) if end_text else int(start_text)
        except (ValueError, TypeError) as error:
            raise ValueError(f"Invalid protected-site token: {token}") from error
        limit = domain_a_length if domain == "A" else domain_b_length if domain == "B" else 0
        if start < 1 or end < start or end > limit:
            raise ValueError(f"Protected-site token is outside its domain: {token}")
        offset = 0 if domain == "A" else domain_a_length + linker_length
        sites.extend({"selector": token, "domain": domain, "domain_residue": value, "fusion_residue": offset + value} for value in range(start, end + 1))
    return sites


def heavy_atoms(residue) -> list:
    return [atom for atom in residue.get_atoms() if atom.element != "H"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--domain-a-length", type=int, required=True)
    parser.add_argument("--linker-length", type=int, required=True)
    parser.add_argument("--sites", required=True, help="Domain-relative syntax, for example A:48;B:27-32")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    structure = PDBParser(QUIET=True).get_structure("fusion", args.pdb)
    chains = list(structure.get_chains())
    if len(chains) != 1:
        raise ValueError("Protected-site assessment requires exactly one chain")
    residues = [residue for residue in chains[0] if residue.id[0] == " " and residue.resname in MAX_ASA]
    domain_b_length = len(residues) - args.domain_a_length - args.linker_length
    if domain_b_length < 1:
        raise ValueError("Domain/linker boundaries do not fit the PDB")
    sites = parse_sites(args.sites, args.domain_a_length, args.linker_length, domain_b_length)
    ShrakeRupley().compute(structure, level="R")
    linker_start = args.domain_a_length
    linker_end = linker_start + args.linker_length
    linker_atoms = [atom for residue in residues[linker_start:linker_end] for atom in heavy_atoms(residue)]
    measurements = []
    for site in sites:
        index = site["fusion_residue"] - 1
        residue = residues[index]
        atoms = heavy_atoms(residue)
        nonlocal_residues = [
            other for other_index, other in enumerate(residues)
            if abs(other_index - index) > 2
        ]
        nonlocal_atoms = [atom for other in nonlocal_residues for atom in heavy_atoms(other)]
        minimum_nonlocal = min(
            (float(np.linalg.norm(left.coord - right.coord)) for left in atoms for right in nonlocal_atoms),
            default=None,
        )
        nearby_residues = sum(
            any(float(np.linalg.norm(left.coord - right.coord)) < 5.0 for left in atoms for right in heavy_atoms(other))
            for other in nonlocal_residues
        )
        linker_minimum = min(
            (float(np.linalg.norm(left.coord - right.coord)) for left in atoms for right in linker_atoms),
            default=None,
        )
        sasa = float(getattr(residue, "sasa", 0.0))
        measurements.append({
            **site,
            "residue_name": residue.resname,
            "residue_sasa_angstrom2": round(sasa, 4),
            "relative_sasa": round(sasa / MAX_ASA[residue.resname], 4),
            "minimum_nonlocal_heavy_atom_distance_angstrom": round(minimum_nonlocal, 4) if minimum_nonlocal is not None else None,
            "nonlocal_residues_within_5_angstrom": nearby_residues,
            "minimum_linker_heavy_atom_distance_angstrom": round(linker_minimum, 4) if linker_minimum is not None else None,
            "linker_contact_under_4_angstrom": linker_minimum is not None and linker_minimum < 4.0,
        })
    config_record = {
        "domain_a_length": args.domain_a_length,
        "linker_length": args.linker_length,
        "domain_b_length": domain_b_length,
        "site_expression": args.sites,
        "sasa_algorithm": "ShrakeRupley default probe/radii",
    }
    config_sha = hashlib.sha256(
        json.dumps(config_record, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    root = Path(__file__).resolve().parents[1]
    commit_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    worktree_modified = bool(subprocess.run(
        ["git", "status", "--porcelain"], cwd=root, check=True,
        capture_output=True, text=True,
    ).stdout.strip())
    result = {
        "schema_version": "linkerforge.protected-site-assessment.v1",
        "pdb": str(args.pdb),
        "pdb_sha256": sha256(args.pdb),
        "source_ids": [f"pdb-sha256:{sha256(args.pdb)}"],
        "boundaries": {"domain_a_length": args.domain_a_length, "linker_length": args.linker_length, "domain_b_length": domain_b_length},
        "site_expression": args.sites,
        "measurements": measurements,
        "software": {"python": platform.python_version(), "biopython": Bio.__version__, "sasa_algorithm": "ShrakeRupley default probe/radii"},
        "provenance": {
            "commit_sha": commit_sha,
            "worktree_modified": worktree_modified,
            "config_sha256": config_sha,
        },
        "evidence_quality": "predicted_structure_geometric_sasa_pre_screen",
        "threshold_status": "not_frozen",
        "claim_boundary": "Geometric/SASA measurements only; ligand, partner, membrane and conformational-ensemble context are not represented.",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
