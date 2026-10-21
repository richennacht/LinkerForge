#!/usr/bin/env python3
"""Run a pinned Reduce/Probe hydrogen-aware severe-overlap screen.

This is a MolProbity-component geometry cross-check on fixed predicted
coordinates.  It is deliberately not labelled as the canonical MolProbity
clashscore and does not perform linker reconstruction or conformational
sampling.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from jsonschema import Draft202012Validator

try:
    from scripts.assess_fusion_surface_burial import ALIASES, read_protein
except ModuleNotFoundError:  # Direct execution adds scripts/, not the repository root.
    from assess_fusion_surface_burial import ALIASES, read_protein


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schemas/reduce-probe-overlap-audit.v1.schema.json"
ATOM_RE = re.compile(r"^\s*(\S)\s+(-?\d+)\s+([A-Za-z0-9]{3})\s+(\S+)\s*$")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        raise ValueError(f"A Windows drive path is required: {resolved}")
    tail = resolved.as_posix().split(":", 1)[1].lstrip("/")
    return f"/mnt/{drive}/{tail}"


def tool_command(manifest: Path, tool: str, *arguments: str) -> list[str]:
    return [
        "wsl.exe", "--exec", "/home/linkerforge/.pixi/bin/pixi", "run",
        "--manifest-path", wsl_path(manifest), tool, *arguments,
    ]


def run(command: list[str], *, stdin: str | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, input=stdin, capture_output=True, text=True, check=False)


def parse_atom(text: str) -> dict[str, object]:
    match = ATOM_RE.match(text)
    if not match:
        raise ValueError(f"Unrecognized Probe atom identifier: {text!r}")
    return {
        "chain": match.group(1),
        "residue_number": int(match.group(2)),
        "residue_name": match.group(3),
        "atom_name": match.group(4),
    }


def atom_key(atom: dict[str, object]) -> tuple[str, int, str, str]:
    return (str(atom["chain"]), int(atom["residue_number"]), str(atom["residue_name"]), str(atom["atom_name"]))


def region(residue_number: int, a_length: int, linker_length: int, total_length: int) -> str:
    if residue_number < 1 or residue_number > total_length:
        raise ValueError(f"Probe residue number outside the verified fusion sequence: {residue_number}")
    if residue_number <= a_length:
        return "domain_a"
    if residue_number <= a_length + linker_length:
        return "linker"
    return "domain_b"


def parse_probe_records(text: str, a_length: int, linker_length: int, total_length: int) -> list[dict[str, object]]:
    pairs: dict[tuple[tuple[str, int, str, str], tuple[str, int, str, str]], dict[str, object]] = {}
    for line in text.splitlines():
        if not line.startswith(":"):
            continue
        fields = line.split(":")
        if len(fields) < 9 or fields[2] not in {"bo", "wo"}:
            raise ValueError(f"Unrecognized Probe record: {line}")
        source, target = parse_atom(fields[3]), parse_atom(fields[4])
        source_key, target_key = atom_key(source), atom_key(target)
        key = tuple(sorted((source_key, target_key)))
        surface_gap = float(fields[6])
        record = {
            "contact_type": fields[2],
            "source_atom": source,
            "target_atom": target,
            "condensed_dot_count": int(fields[5]),
            "minimum_surface_gap_angstrom": surface_gap,
            "overlap_magnitude_angstrom": round(max(0.0, -surface_gap), 4),
            "region_pair": "--".join(sorted((
                region(int(source["residue_number"]), a_length, linker_length, total_length),
                region(int(target["residue_number"]), a_length, linker_length, total_length),
            ))),
        }
        previous = pairs.get(key)
        if previous is None or surface_gap < float(previous["minimum_surface_gap_angstrom"]):
            pairs[key] = record
    return sorted(pairs.values(), key=lambda item: (float(item["minimum_surface_gap_angstrom"]), str(item["region_pair"])))


def count_explicit_atoms(pdb_text: str) -> int:
    seen: set[tuple[str, str, str]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM  ") or line[16:17] not in {" ", "A"}:
            continue
        seen.add((line[21:22], line[22:27].strip(), line[12:16].strip()))
    if not seen:
        raise ValueError("Reduce output contains no explicit protein atoms")
    return len(seen)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--sequence-manifest", required=True, type=Path)
    parser.add_argument("--pixi-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite an existing run directory: {args.output_dir}")
    comparison = json.loads(args.comparison.read_text(encoding="utf-8"))
    source = json.loads(args.sequence_manifest.read_text(encoding="utf-8"))
    if comparison.get("status") != "complete" or comparison.get("schema_version") not in {
        "linkerforge.matched-structure-comparison.v1",
        "linkerforge.public-multipair-structure-comparison.v1",
    }:
        raise ValueError("Require a supported complete structure comparison")
    source_rows = {row["construct_id"]: row for row in source["records"]}
    if len(source_rows) != len(source["records"]):
        raise ValueError("Duplicate source construct IDs")
    intake: list[tuple[str, dict, Path]] = []
    for row in comparison["records"]:
        pair = row.get("domain_pair", "1UBQ_2GB1")
        linker = row.get("linker_id", row.get("candidate_id"))
        identifier = f"{pair}__{ALIASES.get(linker, linker)}"
        record = source_rows[identifier]
        pdb = Path(row["pdb"])
        if not pdb.is_absolute():
            pdb = ROOT / pdb
        if sha256(pdb) != row["pdb_sha256"]:
            raise ValueError(f"PDB hash mismatch: {identifier}")
        if hashlib.sha256(record["sequence"].encode()).hexdigest() != record["sequence_sha256"]:
            raise ValueError(f"Sequence hash mismatch: {identifier}")
        if record["linker_length"] != row["linker_length"]:
            raise ValueError(f"Linker boundary mismatch: {identifier}")
        read_protein(pdb, record["sequence"])
        intake.append((identifier, record, pdb))
    if not intake or len({row[0] for row in intake}) != len(intake):
        raise ValueError("Comparison must contain unique nonempty constructs")
    version_commands = {
        "reduce": tool_command(args.pixi_manifest, "reduce", "-version"),
        "probe": tool_command(args.pixi_manifest, "probe", "-version"),
    }
    versions: dict[str, str] = {}
    for tool, command in version_commands.items():
        completed = run(command)
        version_text = (completed.stdout or completed.stderr).strip()
        # Reduce 4.15 prints a valid version to stderr and exits 2 for -version.
        # Validate the identity string instead of assuming Unix success semantics.
        if not version_text.lower().startswith(f"{tool}."):
            raise RuntimeError(
                f"{tool} version check failed (exit {completed.returncode}): "
                f"{completed.stderr.strip()}"
            )
        versions[tool] = version_text

    args.output_dir.mkdir(parents=True, exist_ok=False)
    reduced_dir, raw_dir = args.output_dir / "reduced", args.output_dir / "probe_raw"
    reduced_dir.mkdir()
    raw_dir.mkdir()
    records = []
    invocation_records = []
    for identifier, record, pdb in intake:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", identifier)
        reduce_command = tool_command(args.pixi_manifest, "reduce", "-quiet", "-build", wsl_path(pdb))
        reduced = run(reduce_command)
        if reduced.returncode != 0 or not any(line.startswith("ATOM  ") for line in reduced.stdout.splitlines()):
            raise RuntimeError(f"Reduce failed for {identifier}: {reduced.stderr.strip()}")
        reduced_path = reduced_dir / f"{safe}.pdb"
        reduced_path.write_text(reduced.stdout, encoding="utf-8", newline="\n")
        probe_command = tool_command(
            args.pixi_manifest, "probe", "-quiet", "-Unformated", "-condense", "-onlybadout",
            "-self", "-mc", "-nohet", "-nowater", "-stdbonds", "all", wsl_path(reduced_path),
        )
        probed = run(probe_command)
        if probed.returncode != 0:
            raise RuntimeError(f"Probe failed for {identifier}: {probed.stderr.strip()}")
        raw_path = raw_dir / f"{safe}.txt"
        raw_path.write_text(probed.stdout, encoding="utf-8", newline="\n")
        overlaps = parse_probe_records(
            probed.stdout, record["domain_a_length"], record["linker_length"], len(record["sequence"])
        )
        atom_count = count_explicit_atoms(reduced.stdout)
        region_counts: dict[str, int] = {}
        for overlap in overlaps:
            region_counts[str(overlap["region_pair"])] = region_counts.get(str(overlap["region_pair"]), 0) + 1
        records.append({
            "construct_id": identifier,
            "input_pdb": str(pdb),
            "input_pdb_sha256": sha256(pdb),
            "reduced_pdb": str(reduced_path),
            "reduced_pdb_sha256": sha256(reduced_path),
            "probe_raw": str(raw_path),
            "probe_raw_sha256": sha256(raw_path),
            "sequence_sha256": record["sequence_sha256"],
            "boundaries": {
                "domain_a_length": record["domain_a_length"],
                "linker_length": record["linker_length"],
                "domain_b_length": len(record["sequence"]) - record["domain_a_length"] - record["linker_length"],
            },
            "explicit_protein_atom_count": atom_count,
            "unique_severe_overlap_pair_count": len(overlaps),
            "unique_severe_overlap_pairs_per_1000_explicit_atoms": round(1000 * len(overlaps) / atom_count, 4),
            "worst_overlap_magnitude_angstrom": round(max((float(row["overlap_magnitude_angstrom"]) for row in overlaps), default=0.0), 4),
            "region_pair_counts": region_counts,
            "overlaps": overlaps,
        })
        invocation_records.append({
            "construct_id": identifier,
            "reduce": {"arguments": reduce_command, "exit_code": reduced.returncode, "stderr": reduced.stderr},
            "probe": {"arguments": probe_command, "exit_code": probed.returncode, "stderr": probed.stderr},
        })
        print(f"Screened {identifier}: {len(overlaps)} unique severe pairs", flush=True)

    config = {
        "selection": "all protein atoms after Reduce; no HETATM or waters",
        "probe_flags": ["-Unformated", "-condense", "-onlybadout", "-self", "-mc", "-nohet", "-nowater", "-stdbonds"],
        "comparison_sha256": sha256(args.comparison),
        "sequence_manifest_sha256": sha256(args.sequence_manifest),
        "pixi_manifest_sha256": sha256(args.pixi_manifest),
        "pixi_lock_sha256": sha256(args.pixi_manifest.with_name("pixi.lock")),
    }
    result = {
        "schema_version": "linkerforge.reduce-probe-overlap-audit.v1",
        "status": "complete",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "records": records,
        "source_ids": [
            f"comparison-sha256:{config['comparison_sha256']}",
            f"sequence-manifest-sha256:{config['sequence_manifest_sha256']}",
        ],
        "config": config,
        "provenance": {
            "commit_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "worktree_modified": bool(subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, text=True).strip()),
            "config_sha256": hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest(),
            "script_sha256": sha256(Path(__file__)),
            "executable": sys.executable,
            "arguments": sys.argv,
            "platform": platform.platform(),
            "versions": versions,
            "invocations": invocation_records,
        },
        "evidence_quality": "predicted_structure_hydrogen_aware_geometry_crosscheck",
        "threshold_status": "descriptive_no_linkerforge_acceptance_threshold",
        "limitations": [
            "This reports Probe-defined severe overlaps after Reduce hydrogen placement on fixed predicted coordinates.",
            "The transparent pair rate is not the canonical MolProbity clashscore and must not be presented as one.",
            "Hydrogen placement and side-chain flips can change contacts; ligands, partners and waters are excluded.",
            "No linker reconstruction, conformational sampling, Rosetta KIC, energetic ranking or experimental validation is performed.",
            "No frozen LinkerForge pass/fail threshold or candidate winner is inferred.",
        ],
        "claim_boundary": "Complementary hydrogen-aware fixed-coordinate overlap screen only; not full MolProbity, Rosetta, stability or compatibility evidence.",
    }
    Draft202012Validator(json.loads(SCHEMA.read_text(encoding="utf-8"))).validate(result)
    output = args.output_dir / "reduce_probe_overlap_audit.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    report = [
        "# Reduce/Probe hydrogen-aware severe-overlap screen", "", result["claim_boundary"], "",
        "| Construct | Explicit atoms | Unique severe pairs | Pairs / 1000 explicit atoms | Worst overlap (Å) |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    report.extend(
        f"| {row['construct_id']} | {row['explicit_protein_atom_count']} | {row['unique_severe_overlap_pair_count']} | "
        f"{row['unique_severe_overlap_pairs_per_1000_explicit_atoms']:.2f} | {row['worst_overlap_magnitude_angstrom']:.3f} |"
        for row in records
    )
    report += ["", "These descriptive counts have no frozen LinkerForge acceptance threshold and are not canonical MolProbity clashscores.", ""]
    report += [f"- {item}" for item in result["limitations"]]
    (args.output_dir / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
