#!/usr/bin/env python3
"""CPU-only surface-burial measurements on sequence-verified fusion structures.

Compare four atom sets at identical coordinates: full fusion, both domains
without the linker, isolated Domain A, and isolated Domain B. No relaxation,
reference-state experiment, functional annotation, or success label is inferred.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import Bio
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schemas/fusion-surface-burial.v1.schema.json"
ALIASES = {"eaaak_15_rigid": "eaaak_15", "ggggs_15_fixed": "ggggs_15"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_protein(path: Path, sequence: str):
    structure = PDBParser(QUIET=False).get_structure("fusion", path)
    models = list(structure)
    if len(models) != 1 or len(list(models[0])) != 1:
        raise ValueError("Require exactly one model and one chain")
    residues = list(next(models[0].get_chains()))
    if any(r.id[0] != " " or seq1(r.resname) == "X" for r in residues):
        raise ValueError("Only canonical protein ATOM residues are supported")
    if "".join(seq1(r.resname) for r in residues) != sequence:
        raise ValueError("PDB sequence differs from the input manifest")
    for residue in residues:
        if not {"N", "CA", "C", "O"}.issubset(residue.child_dict):
            raise ValueError("Backbone atoms missing; cannot assess surface reliably")
        for atom in list(residue):
            if atom.is_disordered() or not np.isfinite(atom.coord).all():
                raise ValueError("Ambiguous alternate locations or nonfinite coordinates")
            if atom.element in {"H", "D"}:
                residue.detach_child(atom.id)
            elif atom.element not in {"C", "N", "O", "S"}:
                raise ValueError(f"Unsupported protein atom element: {atom.element}")
    return structure


def subset_sasa(structure, indices: set[int], n_points: int) -> dict[int, float]:
    selected = copy.deepcopy(structure)
    chain = next(selected.get_chains())
    mapping = {}
    for index, residue in enumerate(list(chain)):
        if index in indices:
            mapping[index] = residue
        else:
            chain.detach_child(residue.id)
    if not mapping:
        return {}
    ShrakeRupley(probe_radius=1.4, n_points=n_points).compute(selected, level="R")
    return {index: float(residue.sasa) for index, residue in mapping.items()}


def measure(structure, a_length: int, linker_length: int, n_points: int) -> list[dict]:
    residues = list(structure.get_residues())
    b_start = a_length + linker_length
    if a_length < 1 or linker_length < 0 or b_start >= len(residues):
        raise ValueError("Invalid domain/linker boundaries")
    all_indices = set(range(len(residues)))
    a_indices, b_indices = set(range(a_length)), set(range(b_start, len(residues)))
    full = subset_sasa(structure, all_indices, n_points)
    no_linker = subset_sasa(structure, a_indices | b_indices, n_points)
    isolated = subset_sasa(structure, a_indices, n_points) | subset_sasa(structure, b_indices, n_points)
    rows = []
    for index in sorted(a_indices | b_indices):
        linker_loss = no_linker[index] - full[index]
        total_loss = isolated[index] - full[index]
        if linker_loss < -1e-6 or total_loss < linker_loss - 1e-6:
            raise ValueError("Surface-subset monotonicity check failed")
        rows.append({
            "domain": "A" if index < a_length else "B",
            "domain_residue": index + 1 if index < a_length else index - b_start + 1,
            "fusion_residue": index + 1,
            "residue_name": residues[index].resname,
            "full_fusion_sasa_angstrom2": round(full[index], 4),
            "without_linker_sasa_angstrom2": round(no_linker[index], 4),
            "isolated_domain_same_coordinates_sasa_angstrom2": round(isolated[index], 4),
            "linker_attributed_burial_angstrom2": round(max(0.0, linker_loss), 4),
            "total_added_context_burial_angstrom2": round(max(0.0, total_loss), 4),
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--sequence-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--n-points", type=int, default=960)
    args = parser.parse_args()
    if args.n_points < 100:
        raise ValueError("Use at least 100 sphere points")
    comparison = json.loads(args.comparison.read_text(encoding="utf-8"))
    source = json.loads(args.sequence_manifest.read_text(encoding="utf-8"))
    if comparison.get("status") != "complete":
        raise ValueError("Require a complete comparison record")
    if comparison.get("schema_version") not in {
        "linkerforge.matched-structure-comparison.v1",
        "linkerforge.public-multipair-structure-comparison.v1",
    }:
        raise ValueError("Unsupported comparison schema")
    source_rows = {row["construct_id"]: row for row in source["records"]}
    if len(source_rows) != len(source["records"]):
        raise ValueError("Duplicate source construct IDs")
    # Verify the entire intake before creating any outputs.
    intake = []
    for row in comparison["records"]:
        if "domain_pair" not in row and comparison.get("analysis_config", {}).get("domain_a_length") != 76:
            raise ValueError("Legacy single-pair records must describe the documented 1UBQ/2GB1 panel")
        pair = row.get("domain_pair", "1UBQ_2GB1")
        linker = row.get("linker_id", row.get("candidate_id"))
        identifier = f"{pair}__{ALIASES.get(linker, linker)}"
        record = source_rows[identifier]
        a_length = record.get("domain_a_length")
        linker_length = record.get("linker_length")
        if not isinstance(a_length, int) or not isinstance(linker_length, int):
            raise ValueError(f"Non-integer boundaries: {identifier}")
        b_length = len(record["sequence"]) - a_length - linker_length
        expected_start = a_length + 1 if linker_length else None
        expected_end = a_length + linker_length if linker_length else None
        if a_length < 1 or linker_length < 0 or b_length < 1:
            raise ValueError(f"Invalid domain/linker boundaries: {identifier}")
        if record.get("linker_start_1based") != expected_start or record.get("linker_end_1based") != expected_end:
            raise ValueError(f"Linker coordinates disagree with boundaries: {identifier}")
        pdb = Path(row["pdb"])
        if not pdb.is_absolute():
            pdb = ROOT / pdb
        if sha256(pdb) != row["pdb_sha256"]:
            raise ValueError(f"PDB hash mismatch: {identifier}")
        if hashlib.sha256(record["sequence"].encode()).hexdigest() != record["sequence_sha256"]:
            raise ValueError(f"Sequence hash mismatch: {identifier}")
        if record["linker_length"] != row["linker_length"]:
            raise ValueError(f"Linker boundary mismatch: {identifier}")
        intake.append((identifier, record, pdb, read_protein(pdb, record["sequence"])))
    if not intake or len({item[0] for item in intake}) != len(intake):
        raise ValueError("Comparison must contain unique nonempty constructs")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    config = {"n_points": args.n_points, "probe_radius_angstrom": 1.4,
              "atom_policy": "canonical protein heavy atoms; reject alternate locations",
              "comparison_sha256": sha256(args.comparison),
              "sequence_manifest_sha256": sha256(args.sequence_manifest)}
    rows = []
    for identifier, record, pdb, structure in intake:
        measurements = measure(structure, record["domain_a_length"], record["linker_length"], args.n_points)
        domain_b_length = len(record["sequence"]) - record["domain_a_length"] - record["linker_length"]
        rows.append({"construct_id": identifier, "pdb": str(pdb), "pdb_sha256": sha256(pdb),
                     "sequence_sha256": record["sequence_sha256"], "domain_sources": record["domain_sources"],
                     "boundaries": {
                         "domain_a_length": record["domain_a_length"],
                         "linker_length": record["linker_length"],
                         "domain_b_length": domain_b_length,
                         "linker_start_1based": record.get("linker_start_1based"),
                         "linker_end_1based": record.get("linker_end_1based"),
                     },
                     "measurements": measurements,
                     "domain_surface_buried_by_linker_angstrom2": round(sum(m["linker_attributed_burial_angstrom2"] for m in measurements), 4),
                     "domain_surface_buried_by_added_context_angstrom2": round(sum(m["total_added_context_burial_angstrom2"] for m in measurements), 4)})
        print(f"Measured {identifier}", flush=True)
    result = {
        "schema_version": "linkerforge.fusion-surface-burial.v1", "status": "complete",
        "created_utc": datetime.now(timezone.utc).isoformat(), "records": rows,
        "evidence_quality": "predicted_structure_fixed_coordinate_surface_measurement",
        "source_ids": [
            f"matched-comparison-sha256:{sha256(args.comparison)}",
            f"sequence-manifest-sha256:{sha256(args.sequence_manifest)}",
        ],
        "config": config,
        "provenance": {
            "commit_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "worktree_modified": bool(subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, text=True).strip()),
            "config_sha256": hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest(),
            "script_sha256": sha256(Path(__file__)), "executable": sys.executable, "arguments": sys.argv,
            "versions": {"python": platform.python_version(), "biopython": Bio.__version__, "numpy": np.__version__},
        },
        "threshold_status": "not_frozen",
        "limitations": [
            "The isolated domain uses the fusion conformation; no reference structure or relaxation is implied.",
            "Surface burial describes this predicted pose; uncertainty in domain placement remains.",
            "No functional/protected-site labels are inferred. Residue tables require researcher annotations.",
            "Heavy atoms only; missing sidechains, ligand, partner, solvent and conformational ensembles affect interpretation.",
            "Domain totals sum per-residue lost surface; they are not divided by two as an interface-area convention.",
        ],
        "claim_boundary": "Geometric surface exposure and burial only; no function, stability or compatibility probability.",
    }
    Draft202012Validator(json.loads(SCHEMA.read_text(encoding="utf-8"))).validate(result)
    (args.output_dir / "surface_burial.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    report = ["# Fixed-coordinate surface-burial measurements", "", result["claim_boundary"], "",
              "| Construct | Domain surface buried by linker (Å²) | By linker + other domain (Å²) |", "| --- | ---: | ---: |"]
    report.extend(f"| {row['construct_id']} | {row['domain_surface_buried_by_linker_angstrom2']:.2f} | {row['domain_surface_buried_by_added_context_angstrom2']:.2f} |" for row in rows)
    report += ["", "These are descriptive measurements with no acceptance cutoff. Lower burial is not automatically a better fusion.", ""]
    report += [f"- {limitation}" for limitation in result["limitations"]]
    (args.output_dir / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
