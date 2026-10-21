#!/usr/bin/env python3
"""Summarize replicated OpenMM structural metrics without inventing a pass threshold."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from jsonschema import Draft202012Validator


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schemas/openmm-convergence-summary.v1.schema.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mean(values: list[float]) -> float:
    return statistics.fmean(values)


def summarize_rows(rows: list[dict[str, str]], blocks: int) -> dict:
    if len(rows) < blocks or len(rows) % blocks:
        raise ValueError(f"Metric frame count {len(rows)} is not divisible into {blocks} blocks")
    steps = [int(row["step"]) for row in rows]
    if steps != sorted(set(steps)):
        raise ValueError("Structural metric steps must be unique and increasing")
    keys = [key for key in rows[0] if key != "step"]
    width = len(rows) // blocks
    metrics = {}
    for key in keys:
        values = [float(row[key]) for row in rows]
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"Non-finite metric: {key}")
        block_means = [mean(values[index * width:(index + 1) * width]) for index in range(blocks)]
        metrics[key] = {
            "overall_mean": round(mean(values), 6),
            "first_half_mean": round(mean(values[:len(values) // 2]), 6),
            "second_half_mean": round(mean(values[len(values) // 2:]), 6),
            "half_mean_absolute_drift": round(abs(mean(values[len(values) // 2:]) - mean(values[:len(values) // 2])), 6),
            "block_means": [round(value, 6) for value in block_means],
            "block_mean_range": round(max(block_means) - min(block_means), 6),
        }
    return {"frame_count": len(rows), "first_step": steps[0], "last_step": steps[-1], "metrics": metrics}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=5)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite convergence summary: {args.output}")
    manifest_path = args.run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete" or len(manifest.get("replicates", [])) < 2:
        raise ValueError("Require a complete replicated OpenMM run")
    replicate_summaries = []
    for replicate in manifest["replicates"]:
        number = int(replicate["replicate"])
        metrics_path = args.run_dir / f"replicate-{number:03d}" / "structural-metrics.csv"
        if sha256(metrics_path) != replicate["structural_metrics_sha256"]:
            raise ValueError(f"Structural-metrics hash mismatch for replicate {number}")
        with metrics_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        summary = summarize_rows(rows, args.blocks)
        summary.update(replicate=number, seed=int(replicate["seed"]), metrics_path=str(metrics_path), metrics_sha256=sha256(metrics_path))
        replicate_summaries.append(summary)
    metric_names = list(replicate_summaries[0]["metrics"])
    across = {}
    for name in metric_names:
        replicate_means = [float(row["metrics"][name]["overall_mean"]) for row in replicate_summaries]
        grand_mean = mean(replicate_means)
        across[name] = {
            "replicate_means": replicate_means,
            "mean_of_replicate_means": round(grand_mean, 6),
            "replicate_mean_range": round(max(replicate_means) - min(replicate_means), 6),
            "replicate_mean_coefficient_of_variation": round(statistics.pstdev(replicate_means) / abs(grand_mean), 6) if grand_mean else None,
        }
    config = {"blocks_per_replicate": args.blocks, "acceptance_thresholds": "not_frozen"}
    result = {
        "schema_version": "linkerforge.openmm-convergence-summary.v1",
        "status": "complete_descriptive_only",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_manifest": str(manifest_path),
        "run_manifest_sha256": sha256(manifest_path),
        "config": config,
        "replicates": replicate_summaries,
        "across_replicates": across,
        "threshold_status": "not_frozen_no_convergence_verdict",
        "provenance": {
            "commit_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "script_sha256": sha256(Path(__file__)),
            "config_sha256": hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest(),
        },
        "limitations": [
            "Block drift and between-replicate dispersion are descriptive because no acceptance thresholds were preregistered.",
            "Stored coordinates have lower temporal resolution than structural metrics under the storage-aware protocol.",
            "A completed finite simulation does not establish equilibrium, biological stability, activity or experimental success.",
        ],
        "claim_boundary": "Descriptive replicated trajectory diagnostics only; no convergence pass, stability claim or biological validation.",
    }
    Draft202012Validator(json.loads(SCHEMA.read_text(encoding="utf-8"))).validate(result)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
