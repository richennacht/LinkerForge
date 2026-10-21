#!/usr/bin/env python3
"""Fail-closed Rosetta linker-remodelling assessment adapter.

Rosetta is separately licensed and is not downloaded by this repository.  This
adapter records exact inputs and parses score terms only after a real Rosetta
execution succeeds; it never substitutes a synthetic score.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_scorefile(path: Path) -> list[dict]:
    header: list[str] | None = None
    records = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("SCORE:"):
            continue
        values = line.split()[1:]
        if "description" in values:
            header = values
            continue
        if header is None or len(values) != len(header):
            continue
        record = {}
        for key, value in zip(header, values, strict=True):
            try:
                record[key] = float(value)
            except ValueError:
                record[key] = value
        records.append(record)
    if not records:
        raise ValueError("Rosetta produced no parseable SCORE records")
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rosetta-scripts", type=Path, required=True)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--linker-start", type=int, required=True)
    parser.add_argument("--linker-end", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not args.rosetta_scripts.is_file():
        raise SystemExit("Rosetta executable is unavailable; obtain it under the appropriate licence before running")
    if args.output_dir.exists():
        raise SystemExit("Rosetta output directory already exists; refusing to overwrite it")
    if args.linker_start < 2 or args.linker_end < args.linker_start:
        raise SystemExit("Invalid linker residue range")
    config = json.loads(args.config.read_text(encoding="utf-8"))
    args.output_dir.mkdir(parents=True)
    loop_file = args.output_dir / "linker.loop"
    cut = (args.linker_start + args.linker_end) // 2
    loop_file.write_text(f"LOOP {args.linker_start} {args.linker_end} {cut} 0 1\n", encoding="utf-8")
    command = [
        str(args.rosetta_scripts), "-s", str(args.pdb),
        "-parser:protocol", str(args.protocol),
        "-loops:loop_file", str(loop_file),
        "-nstruct", str(int(config["nstruct"])),
        "-out:path:all", str(args.output_dir),
        "-out:file:scorefile", str(args.output_dir / "score.sc"),
        "-constant_seed", "-jran", str(int(config["seed"])),
        "-overwrite", "false",
    ]
    if args.database:
        command.extend(["-database", str(args.database)])
    started = datetime.now(timezone.utc).isoformat()
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    (args.output_dir / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (args.output_dir / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise SystemExit(f"Rosetta failed with exit code {completed.returncode}; inspect retained logs")
    records = parse_scorefile(args.output_dir / "score.sc")
    required_terms = config["required_score_terms"]
    missing = sorted({term for term in required_terms if any(term not in record for record in records)})
    if missing:
        raise SystemExit(f"Rosetta score file omits required terms: {', '.join(missing)}")
    thresholds = config["hard_rejection_maxima"]
    evaluated = []
    for record in records:
        reasons = [f"{term}>{maximum}" for term, maximum in thresholds.items() if float(record[term]) > float(maximum)]
        evaluated.append({"description": record["description"], "scores": record, "hard_rejection_pass": not reasons, "rejection_reasons": reasons})
    root = Path(__file__).resolve().parents[1]
    commit_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    worktree_modified = bool(subprocess.run(
        ["git", "status", "--porcelain"], cwd=root, check=True,
        capture_output=True, text=True,
    ).stdout.strip())
    manifest = {
        "schema_version": "linkerforge.rosetta-assessment.v1",
        "status": "complete",
        "started_at": started,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "input_pdb_sha256": sha256(args.pdb),
        "source_ids": [
            f"pdb-sha256:{sha256(args.pdb)}",
            f"protocol-sha256:{sha256(args.protocol)}",
            f"rosetta-executable-sha256:{sha256(args.rosetta_scripts)}",
        ],
        "protocol_sha256": sha256(args.protocol),
        "config_sha256": sha256(args.config),
        "linker_range": [args.linker_start, args.linker_end],
        "models": evaluated,
        "passing_models": sum(item["hard_rejection_pass"] for item in evaluated),
        "command": command,
        "provenance": {
            "commit_sha": commit_sha,
            "worktree_modified": worktree_modified,
            "config_sha256": sha256(args.config),
            "versions": {
                "python": platform.python_version(),
                "adapter": "1.0.0",
                "rosetta_executable_sha256": sha256(args.rosetta_scripts),
            },
        },
        "evidence_quality": "rosetta_computational_conformational_energy_screen",
        "claim_boundary": "Rosetta conformational/energy screen under frozen thresholds; not experimental validation.",
    }
    (args.output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
