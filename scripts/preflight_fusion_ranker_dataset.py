#!/usr/bin/env python3
"""Fail-closed eligibility check for a future fusion compatibility ranker.

This gate is deliberately separate from the active DFL detector firewall.  It
does not turn mined linker calls or synthetic negatives into experimental
fusion-success labels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import platform
import subprocess


ALLOWED_LABEL_CLASSES = {"natural_observation", "synthetic_negative", "engineered_experimental_outcome"}
ALLOWED_SPLITS = {"train", "development", "blind"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    failures: list[str] = []
    manifest_path = args.dataset_dir / "manifest.json"
    records_path = args.dataset_dir / "records.jsonl"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    inventory_counts: dict[str, int] = {}
    source_inventory = manifest.get("dataset_id") == "mined_and_natural_linkers_v1"
    if source_inventory:
        records_path = args.dataset_dir / str(manifest.get("files", {}).get("all_records_jsonl", {}).get("path", "all_records.jsonl"))
        if not records_path.is_file():
            failures.append("published source inventory is missing all_records.jsonl")
            records = []
        else:
            expected = manifest["files"]["all_records_jsonl"]
            if sha256(records_path) != expected["sha256"]:
                failures.append("source-inventory all_records.jsonl hash differs from manifest")
            records = [json.loads(line) for line in records_path.read_text(encoding="utf-8").splitlines() if line.strip()]
            if len(records) != expected["records"]:
                failures.append("source-inventory record count differs from manifest")
            inventory_counts = dict(Counter(str(record["authority_tier"]) for record in records))
        failures.extend([
            "published composite is a provenance-stratified source inventory, not a curator-reviewed fusion-ranker dataset",
            "model-mined candidate tiers are predictions and cannot be promoted to natural or experimental success labels",
            "ranker-specific family-disjoint train/development/blind assignments are absent",
            "distinctly labelled context-matched synthetic negatives are absent",
        ])
    elif not manifest_path.is_file() or not records_path.is_file():
        failures.append("dataset must contain manifest.json and records.jsonl")
        records = []
    else:
        records = [json.loads(line) for line in records_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if manifest.get("status") != "curator_reviewed_family_disjoint_fusion_ranker_dataset":
        failures.append("manifest status is not curator_reviewed_family_disjoint_fusion_ranker_dataset")
    if manifest.get("records_sha256") and manifest["records_sha256"] != sha256(records_path):
        failures.append("records.jsonl hash differs from manifest")
    clusters: dict[str, set[str]] = defaultdict(set)
    label_counts: dict[str, int] = defaultdict(int)
    required = {
        "record_id", "sequence", "domain_a_source_id", "domain_b_source_id",
        "linker_sequence", "label_class", "evidence_quality", "homology_cluster_id",
        "split", "source_release", "source_ids", "provenance", "config_hash",
    }
    for index, record in enumerate([] if source_inventory else records, start=1):
        missing = sorted(required - record.keys())
        if missing:
            failures.append(f"record {index} missing: {', '.join(missing)}")
            continue
        if record["label_class"] not in ALLOWED_LABEL_CLASSES:
            failures.append(f"record {record['record_id']} has invalid label_class")
        if record["split"] not in ALLOWED_SPLITS:
            failures.append(f"record {record['record_id']} has invalid split")
        clusters[str(record["homology_cluster_id"])].add(str(record["split"]))
        label_counts[str(record["label_class"])] += 1
        if record["label_class"] == "synthetic_negative" and record["evidence_quality"] == "experimental_failure":
            failures.append(f"synthetic record {record['record_id']} is mislabeled as an experimental failure")
    leaked = sorted(cluster for cluster, splits in clusters.items() if len(splits) > 1)
    if leaked:
        failures.append(f"homology clusters cross splits: {', '.join(leaked[:10])}")
    if not source_inventory and label_counts["natural_observation"] == 0:
        failures.append("no natural observations")
    if not source_inventory and label_counts["synthetic_negative"] == 0:
        failures.append("no distinctly labelled synthetic negatives")
    gate_config = {
        "allowed_label_classes": sorted(ALLOWED_LABEL_CLASSES),
        "allowed_splits": sorted(ALLOWED_SPLITS),
        "required_record_fields": sorted(required),
    }
    config_hash = hashlib.sha256(
        json.dumps(gate_config, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    root = Path(__file__).resolve().parents[1]
    try:
        commit_sha = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, check=True,
            capture_output=True, text=True,
        ).stdout.strip()
        worktree_modified = bool(subprocess.run(
            ["git", "status", "--porcelain"], cwd=root, check=True,
            capture_output=True, text=True,
        ).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        commit_sha = "unavailable"
        worktree_modified = True
    report = {
        "schema_version": "linkerforge.fusion-ranker-preflight.v1",
        "status": "passed" if not failures else "source_inventory_only" if source_inventory else "failed_closed",
        "dataset_dir": str(args.dataset_dir),
        "manifest_sha256": sha256(manifest_path) if manifest_path.is_file() else None,
        "records_sha256": sha256(records_path) if records_path.is_file() else None,
        "record_count": len(records),
        "label_counts": dict(label_counts),
        "source_inventory_counts": inventory_counts,
        "cluster_count": len(clusters),
        "failures": failures,
        "training_authorized": not failures,
        "evidence_quality": "dataset_integrity_and_eligibility_preflight_only",
        "source_ids": [str(manifest.get("dataset_id", args.dataset_dir.name))],
        "provenance": {
            "commit_sha": commit_sha,
            "worktree_modified": worktree_modified,
            "config_sha256": config_hash,
            "versions": {"python": platform.python_version(), "preflight": "1.0.0"},
        },
        "claim_boundary": "A passing data-integrity preflight does not establish biological performance; blind evaluation and experimental evidence remain separate.",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if failures:
        raise SystemExit("Fusion-ranker dataset preflight failed closed")


if __name__ == "__main__":
    main()
