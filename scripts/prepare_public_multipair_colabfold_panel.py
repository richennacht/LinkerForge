#!/usr/bin/env python3
"""Prepare a public, matched, short-fusion ColabFold sensitivity panel.

This selects public sequences from the existing detector-sensitivity manifest.
It never reads protected benchmark labels and does not create compatibility or
experimental-success labels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from jsonschema import Draft202012Validator


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schemas" / "public-multipair-colabfold-panel.v1.schema.json"
CANONICAL = set("ACDEFGHIKLMNPQRSTVWY")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_value(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "not-available"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pairs", nargs="+", required=True)
    parser.add_argument("--max-residues", type=int, default=220)
    args = parser.parse_args()

    source_path = args.source_manifest.resolve()
    source = json.loads(source_path.read_text(encoding="utf-8"))
    selected = [record for record in source["records"] if record["domain_pair"] in args.pairs]
    observed_pairs = {record["domain_pair"] for record in selected}
    missing_pairs = sorted(set(args.pairs) - observed_pairs)
    if missing_pairs:
        raise ValueError(f"Requested domain pairs are absent from the source manifest: {missing_pairs}")
    for pair in args.pairs:
        pair_records = [record for record in selected if record["domain_pair"] == pair]
        if len(pair_records) != 7:
            raise ValueError(f"Expected seven matched constructs for {pair}; found {len(pair_records)}")

    records = []
    fasta_parts = []
    for source_record in sorted(selected, key=lambda row: row["construct_id"]):
        sequence = source_record["sequence"]
        if not sequence or len(sequence) > args.max_residues or set(sequence) - CANONICAL:
            raise ValueError(f"Sequence gate failed for {source_record['construct_id']}")
        prediction_id = f"LinkerForge_{source_record['domain_pair']}_{source_record['linker_id']}"
        records.append({
            "construct_id": source_record["construct_id"],
            "prediction_id": prediction_id,
            "domain_pair": source_record["domain_pair"],
            "domain_sources": source_record["domain_sources"],
            "domain_a_length": source_record["domain_a_length"],
            "domain_b_length": len(sequence) - source_record["domain_a_length"] - source_record["linker_length"],
            "linker_id": source_record["linker_id"],
            "linker_sequence": source_record["linker_sequence"],
            "linker_length": source_record["linker_length"],
            "linker_class": source_record["linker_class"],
            "full_sequence_length": len(sequence),
            "sequence": sequence,
            "sequence_sha256": source_record["sequence_sha256"],
        })
        fasta_parts.append(f">{prediction_id}\n{sequence}\n")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    fasta = "".join(fasta_parts)
    fasta_path = output_dir / "panel.fasta"
    fasta_path.write_text(fasta, encoding="utf-8")
    config = {
        "pairs": args.pairs,
        "max_residues": args.max_residues,
        "predictor": "ColabFold/AlphaFold2-ptm",
        "models_per_construct": 5,
        "seeds": [0],
        "templates": False,
        "msa_mode": "mmseqs2_uniref_env",
        "rank_by": "plddt",
        "fast_kernels": True,
    }
    manifest = {
        "schema_version": "linkerforge.public-multipair-colabfold-panel.v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "Matched public-sequence structural sensitivity across short domain pairs",
        "evidence_quality": "descriptive_computational_sensitivity_not_biological_validation",
        "records": records,
        "record_count": len(records),
        "planned_rankable_models": len(records) * config["models_per_construct"],
        "planned_settings": config,
        "data_firewall": {
            "protected_benchmark_labels_read": False,
            "public_sequences_only": True,
            "used_for_model_training": False,
            "used_for_threshold_selection": False,
        },
        "provenance": {
            "source_manifest": str(source_path),
            "source_manifest_sha256": sha256(source_path),
            "panel_fasta_sha256": hashlib.sha256(fasta.encode("utf-8")).hexdigest(),
            "commit_sha": git_value("rev-parse", "HEAD"),
            "worktree_modified": bool(git_value("status", "--porcelain", "--untracked-files=no")),
            "config_sha256": hashlib.sha256(json.dumps(config, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
            "versions": {"python": platform.python_version(), "panel_builder": "1.0.0"},
        },
        "selection_policy": "Report every construct; presentation examples may be sorted post hoc by AlphaFold confidence only when labelled descriptive.",
        "claim_boundary": "Confidence ranks structure-model outputs; it does not identify a biologically successful or optimal linker.",
    }
    Draft202012Validator(json.loads(SCHEMA.read_text(encoding="utf-8"))).validate(manifest)
    (output_dir / "panel_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output_dir), "constructs": len(records), "planned_models": manifest["planned_rankable_models"]}))


if __name__ == "__main__":
    main()
