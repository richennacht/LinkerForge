#!/usr/bin/env python3
"""Summarize the complete public multi-pair ColabFold sensitivity panel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from jsonschema import Draft202012Validator


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schemas" / "public-multipair-structure-comparison.v1.schema.json"
ALIASES = {"eaaak_15_rigid": "eaaak_15", "ggggs_15_fixed": "ggggs_15"}
AA3 = dict(zip(
    "ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCEQGHILKMFPSTWYV",
))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def pdb_ca_sequence(path: Path) -> str:
    residues: dict[tuple[str, str], str] = {}
    chains: set[str] = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM  ") or line[12:16].strip() != "CA":
            continue
        residue = AA3.get(line[17:20].strip())
        if residue is None:
            raise ValueError(f"Unknown residue in {path.name}")
        chain = line[21:22]
        chains.add(chain)
        residues[(chain, line[22:27])] = residue
    if len(chains) != 1 or not residues:
        raise ValueError(f"Expected one nonempty chain in {path}")
    return "".join(residues.values())


def matching_file(directory: Path, prediction_id: str, marker: str, suffix: str) -> Path:
    matches = [path for path in directory.glob(f"*{suffix}") if prediction_id in path.name and marker in path.name]
    if len(matches) != 1:
        raise ValueError(f"Expected one {prediction_id} {marker} {suffix}; found {len(matches)}")
    return matches[0]


def metrics(record: dict, prediction_dir: Path) -> dict:
    score_path = matching_file(prediction_dir, record["prediction_id"], "scores_rank_001", ".json")
    pdb_path = matching_file(prediction_dir, record["prediction_id"], "unrelaxed_rank_001", ".pdb")
    if pdb_ca_sequence(pdb_path) != record["sequence"]:
        raise ValueError(f"Rank-1 PDB sequence mismatch for {record['construct_id']}")
    score = json.loads(score_path.read_text(encoding="utf-8"))
    plddt = [float(value) for value in score["plddt"]]
    domain_a_end = record["domain_a_length"]
    domain_b_start = domain_a_end + record["linker_length"]
    if len(plddt) != len(record["sequence"]):
        raise ValueError(f"pLDDT length mismatch for {record['construct_id']}")
    pae = score.get("pae")
    cross = [] if not pae else [
        float(pae[left][right])
        for left in range(domain_a_end)
        for right in range(domain_b_start, len(plddt))
    ] + [
        float(pae[left][right])
        for left in range(domain_b_start, len(plddt))
        for right in range(domain_a_end)
    ]
    return {
        "construct_id": record["construct_id"],
        "domain_pair": record["domain_pair"],
        "domain_sources": record["domain_sources"],
        "linker_id": record["linker_id"],
        "linker_sequence": record["linker_sequence"],
        "linker_length": record["linker_length"],
        "full_sequence_length": len(record["sequence"]),
        "mean_plddt": round(mean(plddt) or 0.0, 4),
        "domain_a_mean_plddt": round(mean(plddt[:domain_a_end]) or 0.0, 4),
        "linker_mean_plddt": None if not record["linker_length"] else round(mean(plddt[domain_a_end:domain_b_start]) or 0.0, 4),
        "domain_b_mean_plddt": round(mean(plddt[domain_b_start:]) or 0.0, 4),
        "ptm": score.get("ptm"),
        "mean_cross_domain_pae_angstrom": None if not cross else round(mean(cross) or 0.0, 4),
        "rank_1_model": pdb_path.stem.split("rank_001_")[-1],
        "pdb": str(pdb_path),
        "pdb_sha256": sha256(pdb_path),
        "score_json": str(score_path),
        "score_json_sha256": sha256(score_path),
    }


def existing_rows(existing: dict, linker_sequences: dict[str, str]) -> list[dict]:
    rows = []
    for row in existing["records"]:
        linker_id = ALIASES.get(row["candidate_id"], row["candidate_id"])
        rows.append({
            "construct_id": f"1UBQ_2GB1__{linker_id}",
            "domain_pair": "1UBQ_2GB1",
            "domain_sources": ["PDB 1UBQ chain A", "PDB 2GB1 chain A"],
            "linker_id": linker_id,
            "linker_sequence": linker_sequences[linker_id],
            "linker_length": row["linker_length"],
            "full_sequence_length": 132 + row["linker_length"],
            "mean_plddt": row["mean_plddt"],
            "domain_a_mean_plddt": row["domain_a_mean_plddt"],
            "linker_mean_plddt": row["linker_mean_plddt"],
            "domain_b_mean_plddt": row["domain_b_mean_plddt"],
            "ptm": row["ptm"],
            "mean_cross_domain_pae_angstrom": row["mean_cross_domain_pae_angstrom"],
            "rank_1_model": Path(row["pdb"]).stem.split("rank_001_")[-1],
            "pdb": row["pdb"],
            "pdb_sha256": row["pdb_sha256"],
            "score_json": row["score_json"],
            "score_json_sha256": row["score_json_sha256"],
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-dir", type=Path, required=True)
    parser.add_argument("--existing-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--csv-output", type=Path, required=True)
    parser.add_argument("--progress", type=Path, required=True)
    args = parser.parse_args()

    panel_dir = args.panel_dir.resolve()
    manifest_path = panel_dir / "panel_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    existing_path = args.existing_summary.resolve()
    existing = json.loads(existing_path.read_text(encoding="utf-8"))
    prediction_dir = panel_dir / "prediction"
    existing_manifest_path = existing_path.parent / "panel_manifest.json"
    existing_manifest = json.loads(existing_manifest_path.read_text(encoding="utf-8"))
    linker_sequences = {
        ALIASES.get(record["candidate_id"], record["candidate_id"]): record["linker_sequence"]
        for record in existing_manifest["records"]
    }
    # The 10-aa public acceptance construct predates the matched-panel manifest.
    linker_sequences["ggggs_10"] = "GGGGSGGGGS"
    expected_linker_ids = {
        ALIASES.get(row["candidate_id"], row["candidate_id"])
        for row in existing["records"]
    }
    if set(linker_sequences) != expected_linker_ids:
        raise ValueError("Existing comparison and panel manifest linker IDs do not match")
    rows = existing_rows(existing, linker_sequences) + [metrics(record, prediction_dir) for record in manifest["records"]]
    if len(rows) != 21:
        raise ValueError(f"Expected 21 constructs across three pairs; found {len(rows)}")
    model_files = list(prediction_dir.glob("*scores_rank_*.json"))
    new_models = len(model_files)
    if new_models != manifest["planned_rankable_models"]:
        raise ValueError(f"Expected {manifest['planned_rankable_models']} new rankable models; found {new_models}")

    ranked = sorted(rows, key=lambda row: row["mean_plddt"], reverse=True)
    for index, row in enumerate(ranked, start=1):
        row["descriptive_mean_plddt_rank"] = index
    best_by_pair = []
    for pair in sorted({row["domain_pair"] for row in rows}):
        best_by_pair.append(max(
            (row for row in rows if row["domain_pair"] == pair and row["linker_length"] > 0),
            key=lambda row: row["mean_plddt"],
        ))
    result = {
        "schema_version": "linkerforge.public-multipair-structure-comparison.v1",
        "status": "complete",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "evidence_quality": "matched_three_public_domain_pair_descriptive_sensitivity",
        "domain_pair_count": 3,
        "construct_count": len(rows),
        "models_per_construct": 5,
        "rankable_model_count": 35 + new_models,
        "records": sorted(rows, key=lambda row: (row["domain_pair"], row["linker_id"])),
        "descriptive_top_by_mean_plddt": ranked[:8],
        "descriptive_top_per_domain_pair": sorted(best_by_pair, key=lambda row: row["domain_pair"]),
        "selection_policy": "Post-hoc descriptive sorting by rank-1 mean pLDDT. No winner threshold was preregistered; all constructs remain reported.",
        "metric_interpretation": {
            "mean_plddt": "AlphaFold local-coordinate confidence averaged across the complete fusion; not compatibility or function.",
            "linker_mean_plddt": "Local confidence for inserted linker residues; low values may be consistent with flexibility but do not prove it.",
            "ptm": "AlphaFold predicted global topology confidence.",
            "mean_cross_domain_pae_angstrom": "Mean bidirectional uncertainty between Domain A and Domain B; lower values indicate more confident relative placement.",
        },
        "provenance": {
            "panel_manifest": str(manifest_path),
            "panel_manifest_sha256": sha256(manifest_path),
            "existing_summary": str(existing_path),
            "existing_summary_sha256": sha256(existing_path),
            "existing_panel_manifest": str(existing_manifest_path),
            "existing_panel_manifest_sha256": sha256(existing_manifest_path),
            "commit_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "worktree_modified": bool(subprocess.check_output(
                ["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, text=True
            ).strip()),
            "versions": {"python": platform.python_version(), "summarizer": "1.0.0"},
        },
        "limitations": [
            "Three deliberately short public domain pairs do not estimate generalization to arbitrary proteins.",
            "Sorting after observing AlphaFold confidence is descriptive presentation selection, not an independent benchmark decision.",
            "AlphaFold confidence cannot establish expression, activity, stability, or an optimal linker.",
            "All structure models come from the AlphaFold2 family through one ColabFold environment and one seed.",
        ],
        "claim_boundary": "This panel demonstrates repeatable structure prediction and confidence extraction across public constructs only.",
    }
    Draft202012Validator(json.loads(SCHEMA.read_text(encoding="utf-8"))).validate(result)
    args.output.resolve().write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    fields = [
        "descriptive_mean_plddt_rank", "construct_id", "domain_pair", "linker_id", "linker_length",
        "full_sequence_length", "mean_plddt", "domain_a_mean_plddt", "linker_mean_plddt",
        "domain_b_mean_plddt", "ptm", "mean_cross_domain_pae_angstrom", "rank_1_model",
    ]
    with args.csv_output.resolve().open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ranked)
    progress = json.loads(args.progress.resolve().read_text(encoding="utf-8"))
    progress.update(state="summarized", updated_utc=datetime.now(timezone.utc).isoformat(), overall_percent=95.0)
    args.progress.resolve().write_text(json.dumps(progress, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"constructs": len(rows), "rankable_models": result["rankable_model_count"], "top": ranked[0]["construct_id"]}))


if __name__ == "__main__":
    main()
