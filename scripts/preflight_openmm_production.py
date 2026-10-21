#!/usr/bin/env python3
"""Estimate OpenMM production runtime/storage from a completed smoke run.

This is an operational gate, not a scientific result.  It deliberately uses
observed files and timestamps from the same construct instead of claiming a
hardware-independent performance estimate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path


GIB = 1024**3


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_time(value: str) -> datetime:
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def build_preflight(
    config_path: Path,
    smoke_manifest_path: Path,
    smoke_run_dir: Path,
    target_parent: Path,
    *,
    storage_multiplier: float = 1.5,
    reserve_gib: float = 5.0,
) -> dict:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    smoke = json.loads(smoke_manifest_path.read_text(encoding="utf-8"))
    if smoke.get("status") != "complete":
        raise ValueError("Smoke manifest is not complete")
    if not smoke.get("replicates"):
        raise ValueError("Smoke manifest contains no completed replicates")

    target_lengths = (
        int(config["domain_a_length"]),
        int(config["linker_length"]),
        int(config["domain_b_length"]),
    )
    smoke_config = smoke["config"]
    smoke_lengths = (
        int(smoke_config["domain_a_length"]),
        int(smoke_config["linker_length"]),
        int(smoke_config["domain_b_length"]),
    )
    if target_lengths != smoke_lengths:
        raise ValueError(
            "Production and smoke configurations describe different construct lengths: "
            f"{target_lengths} != {smoke_lengths}"
        )

    observed_trajectory_bytes = 0
    observed_frames = 0
    observed_production_seconds = 0.0
    for replicate in smoke["replicates"]:
        replicate_number = int(replicate["replicate"])
        trajectory = smoke_run_dir / f"replicate-{replicate_number:03d}" / "trajectory.dcd"
        if not trajectory.is_file():
            raise FileNotFoundError(f"Missing smoke trajectory: {trajectory}")
        observed_trajectory_bytes += trajectory.stat().st_size
        observed_frames += int(replicate.get("trajectory_frames", replicate["structural_metric_frames"]))
        observed_production_seconds += (
            _parse_time(replicate["completed_at"])
            - _parse_time(replicate["production_started_at"])
        ).total_seconds()
    if observed_frames <= 0 or observed_production_seconds <= 0:
        raise ValueError("Smoke run has no usable trajectory frames or elapsed production time")

    production_steps = int(config["production_steps"])
    report_interval = int(config["report_interval_steps"])
    trajectory_report_interval = int(config.get("trajectory_report_interval_steps", report_interval))
    replicates = int(config["replicates"])
    if production_steps % report_interval or production_steps % trajectory_report_interval:
        raise ValueError("Production steps must be divisible by metric and trajectory report intervals")
    projected_frames = (production_steps // trajectory_report_interval) * replicates
    projected_metric_frames = (production_steps // report_interval) * replicates
    bytes_per_frame = observed_trajectory_bytes / observed_frames
    projected_trajectory_bytes = int(bytes_per_frame * projected_frames)

    smoke_steps = int(smoke_config["production_steps"]) * len(smoke["replicates"])
    seconds_per_step = observed_production_seconds / smoke_steps
    target_total_steps = (
        int(config["production_steps"]) + int(config["equilibration_steps"])
    ) * replicates
    projected_runtime_seconds = seconds_per_step * target_total_steps

    usage = shutil.disk_usage(target_parent)
    reserve_bytes = int(reserve_gib * GIB)
    required_free_bytes = int(projected_trajectory_bytes * storage_multiplier) + reserve_bytes
    authorized = usage.free >= required_free_bytes
    status = "pass" if authorized else "blocked_insufficient_storage_headroom"

    return {
        "schema_version": "linkerforge.openmm-production-preflight.v1",
        "status": status,
        "run_authorized": authorized,
        "inputs": {
            "production_config": str(config_path.resolve()),
            "production_config_sha256": sha256(config_path),
            "smoke_manifest": str(smoke_manifest_path.resolve()),
            "smoke_manifest_sha256": sha256(smoke_manifest_path),
            "smoke_run_directory": str(smoke_run_dir.resolve()),
            "target_parent": str(target_parent.resolve()),
        },
        "observations": {
            "smoke_atom_count": int(smoke["atom_count"]),
            "smoke_replicates": len(smoke["replicates"]),
            "smoke_production_steps_total": smoke_steps,
            "smoke_production_seconds_total": round(observed_production_seconds, 3),
            "smoke_trajectory_frames_total": observed_frames,
            "smoke_trajectory_bytes_total": observed_trajectory_bytes,
            "observed_bytes_per_frame": round(bytes_per_frame, 3),
            "observed_seconds_per_step": round(seconds_per_step, 9),
        },
        "projection": {
            "replicates": replicates,
            "production_nanoseconds_per_replicate": (
                production_steps * float(config["timestep_fs"]) / 1_000_000
            ),
            "trajectory_frames_total": projected_frames,
            "trajectory_report_interval_steps": trajectory_report_interval,
            "structural_metric_frames_total": projected_metric_frames,
            "structural_metric_report_interval_steps": report_interval,
            "trajectory_bytes": projected_trajectory_bytes,
            "trajectory_gib": round(projected_trajectory_bytes / GIB, 3),
            "runtime_seconds": round(projected_runtime_seconds, 3),
            "runtime_hours": round(projected_runtime_seconds / 3600, 2),
        },
        "storage_gate": {
            "free_bytes": usage.free,
            "free_gib": round(usage.free / GIB, 3),
            "trajectory_storage_multiplier": storage_multiplier,
            "fixed_reserve_gib": reserve_gib,
            "required_free_bytes": required_free_bytes,
            "required_free_gib": round(required_free_bytes / GIB, 3),
            "shortfall_bytes": max(0, required_free_bytes - usage.free),
            "shortfall_gib": round(max(0, required_free_bytes - usage.free) / GIB, 3),
        },
        "claim_boundary": (
            "This preflight estimates operations from one completed same-construct smoke run. "
            "It is not evidence of trajectory convergence, stability, or biological function."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--smoke-manifest", type=Path, required=True)
    parser.add_argument("--smoke-run-dir", type=Path, required=True)
    parser.add_argument("--target-parent", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--storage-multiplier", type=float, default=1.5)
    parser.add_argument("--reserve-gib", type=float, default=5.0)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output already exists; refusing overwrite: {args.output}")
    result = build_preflight(
        args.config,
        args.smoke_manifest,
        args.smoke_run_dir,
        args.target_parent,
        storage_multiplier=args.storage_multiplier,
        reserve_gib=args.reserve_gib,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["run_authorized"] else 2)


if __name__ == "__main__":
    main()
