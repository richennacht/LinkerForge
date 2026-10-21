#!/usr/bin/env python3
"""Run a durable WSL ColabFold panel with GPU locking and progress records."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

from jsonschema import Draft202012Validator


ROOT = Path(__file__).resolve().parents[1]
RUN_SCHEMA = ROOT / "schemas" / "colabfold-panel-run.v1.schema.json"
CREATE_NO_WINDOW = 0x08000000 if os.name == "nt" else 0


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_manifest(value: dict) -> None:
    Draft202012Validator(json.loads(RUN_SCHEMA.read_text(encoding="utf-8"))).validate(value)


def fasta_headers(path: Path) -> list[str]:
    headers = [line[1:].strip() for line in path.read_text(encoding="utf-8").splitlines() if line.startswith(">")]
    if not headers or len(headers) != len(set(headers)):
        raise ValueError("Panel FASTA must contain unique nonempty headers")
    return headers


def wsl_path(path: Path) -> str:
    result = subprocess.run(
        ["wsl.exe", "--exec", "wslpath", "-a", str(path.resolve())],
        check=True,
        capture_output=True,
        text=True,
        creationflags=CREATE_NO_WINDOW,
    )
    translated = result.stdout.strip()
    if not translated:
        raise RuntimeError(f"WSL returned an empty path for {path}")
    return translated


def count_outputs(output: Path) -> tuple[int, int]:
    completed = len(list(output.glob("*.done.txt"))) if output.exists() else 0
    rankable = len(list(output.glob("*scores_rank_*.json"))) if output.exists() else 0
    return completed, rankable


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--progress", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--owner", required=True)
    parser.add_argument("--timeout-hours", type=float, default=18.0)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume a previously interrupted, hash-matched panel and retain completed .done.txt queries.",
    )
    args = parser.parse_args()

    fasta = args.fasta.resolve()
    output = args.output.resolve()
    progress = args.progress.resolve()
    lock = args.lock.resolve()
    headers = fasta_headers(fasta)
    total_constructs = len(headers)
    models_per_construct = 5
    total_models = total_constructs * models_per_construct
    manifest_path = output.parent / "colabfold_run_manifest.json"
    prior_manifest: dict | None = None
    if output.exists() and any(output.iterdir()):
        if not args.resume:
            raise FileExistsError(f"Refusing to overwrite nonempty prediction directory: {output}")
        if not manifest_path.is_file():
            raise FileNotFoundError(f"Resume requires the prior run manifest: {manifest_path}")
        prior_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        validate_manifest(prior_manifest)
        if prior_manifest["status"] == "complete":
            raise RuntimeError("Refusing to resume a completed ColabFold panel")
        if prior_manifest["input_fasta_sha256"] != sha256(fasta):
            raise RuntimeError("Resume FASTA hash does not match the prior run")
        if Path(prior_manifest["output"]).resolve() != output:
            raise RuntimeError("Resume output path does not match the prior run")
        if prior_manifest["construct_count"] != total_constructs:
            raise RuntimeError("Resume construct count does not match the prior run")
        if prior_manifest["models_per_construct"] != models_per_construct:
            raise RuntimeError("Resume model count does not match the prior run")
    elif args.resume:
        raise FileNotFoundError("Resume was requested but the prediction directory is empty or absent")
    output.mkdir(parents=True, exist_ok=True)
    progress.parent.mkdir(parents=True, exist_ok=True)
    lock.parent.mkdir(parents=True, exist_ok=True)

    lock_contents = (
        f"pid={os.getpid()}\n"
        f"owner={args.owner}\n"
        f"output={output.parent}\n"
        f"created_utc={utc_now()}\n"
    )
    try:
        with lock.open("x", encoding="utf-8") as handle:
            handle.write(lock_contents)
    except FileExistsError as exc:
        raise RuntimeError(
            f"Another LinkerForge GPU job may be running. Inspect {lock} and remove it only after confirming no prediction or simulation process is active."
        ) from exc

    command = [
        "wsl.exe", "--exec", "env",
        "TF_FORCE_UNIFIED_MEMORY=1",
        "TF_FORCE_GPU_ALLOW_GROWTH=true",
        "XLA_PYTHON_CLIENT_MEM_FRACTION=4.0",
        "XLA_PYTHON_CLIENT_PREALLOCATE=false",
        "XLA_PYTHON_CLIENT_ALLOCATOR=platform",
        "colabfold_batch",
        "--model-type", "alphafold2_ptm",
        "--num-models", str(models_per_construct),
        "--num-seeds", "1",
        "--random-seed", "0",
        "--rank", "plddt",
        "--use-fast-kernels",
        wsl_path(fasta),
        wsl_path(output),
    ]
    started_now = utc_now()
    if prior_manifest is None:
        manifest = {
            "schema_version": "linkerforge.colabfold-panel-run.v1",
            "status": "running",
            "started_utc": started_now,
            "input_fasta": str(fasta),
            "input_fasta_sha256": sha256(fasta),
            "output": str(output),
            "construct_count": total_constructs,
            "models_per_construct": models_per_construct,
            "planned_rankable_models": total_models,
            "command": command,
            "environment_overrides": {
                "TF_FORCE_UNIFIED_MEMORY": "1",
                "TF_FORCE_GPU_ALLOW_GROWTH": "true",
                "XLA_PYTHON_CLIENT_MEM_FRACTION": "4.0",
                "XLA_PYTHON_CLIENT_PREALLOCATE": "false",
                "XLA_PYTHON_CLIENT_ALLOCATOR": "platform",
            },
            "provenance": {
                "commit_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "runner_sha256": sha256(Path(__file__)),
                "versions": {"python": platform.python_version()},
            },
            "claim_boundary": "Matched public-sequence structure predictions; not compatibility, function, or experimental validation.",
        }
    else:
        completed_before_resume, rankable_before_resume = count_outputs(output)
        manifest = prior_manifest
        manifest.setdefault("resume_events", []).append({
            "resumed_utc": started_now,
            "prior_status": prior_manifest["status"],
            "completed_constructs_before_resume": completed_before_resume,
            "rankable_models_before_resume": rankable_before_resume,
            "provenance": {
                "commit_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "runner_sha256": sha256(Path(__file__)),
                "versions": {"python": platform.python_version()},
            },
        })
        manifest.update(
            status="running",
            latest_started_utc=started_now,
            command=command,
        )
        manifest.pop("completed_utc", None)
        manifest.pop("failure_reason", None)
    validate_manifest(manifest)
    atomic_json(manifest_path, manifest)
    completed_at_start, rankable_at_start = count_outputs(output)
    atomic_json(progress, {
        "schema_version": "linkerforge.panel-progress.v1",
        "state": "resuming" if args.resume else "starting",
        "updated_utc": utc_now(),
        "completed_constructs": completed_at_start,
        "total_constructs": total_constructs,
        "rankable_models": rankable_at_start,
        "planned_rankable_models": total_models,
        "prediction_percent": round(100 * rankable_at_start / total_models, 1),
        "overall_percent": round(90 * rankable_at_start / total_models, 1),
    })
    log_path = output.parent / "colabfold.log"
    process: subprocess.Popen[bytes] | None = None
    deadline = time.monotonic() + args.timeout_hours * 3600
    try:
        with log_path.open("ab" if args.resume else "wb") as log:
            if args.resume:
                log.write(f"\n=== resume {started_now} ===\n".encode("utf-8"))
                log.flush()
            process = subprocess.Popen(
                command,
                cwd=ROOT,
                stdout=log,
                stderr=subprocess.STDOUT,
                creationflags=CREATE_NO_WINDOW,
            )
            while process.poll() is None:
                if time.monotonic() > deadline:
                    process.kill()
                    raise TimeoutError(f"ColabFold panel exceeded {args.timeout_hours:g} hours")
                completed, rankable = count_outputs(output)
                percent = min(99.0, round(100 * rankable / total_models, 1))
                atomic_json(progress, {
                    "schema_version": "linkerforge.panel-progress.v1",
                    "state": "predicting",
                    "updated_utc": utc_now(),
                    "completed_constructs": completed,
                    "total_constructs": total_constructs,
                    "rankable_models": rankable,
                    "planned_rankable_models": total_models,
                    "prediction_percent": percent,
                    "overall_percent": round(percent * 0.9, 1),
                    "process_id": process.pid,
                })
                time.sleep(15)
        if process.returncode != 0:
            raise RuntimeError(f"ColabFold exited {process.returncode}; see {log_path}")
        completed, rankable = count_outputs(output)
        if completed != total_constructs or rankable != total_models:
            raise RuntimeError(
                f"Incomplete panel: {completed}/{total_constructs} constructs and {rankable}/{total_models} rankable models"
            )
        manifest.update(status="complete", completed_utc=utc_now(), rankable_models=rankable)
        atomic_json(progress, {
            "schema_version": "linkerforge.panel-progress.v1",
            "state": "prediction_complete",
            "updated_utc": utc_now(),
            "completed_constructs": completed,
            "total_constructs": total_constructs,
            "rankable_models": rankable,
            "planned_rankable_models": total_models,
            "prediction_percent": 100.0,
            "overall_percent": 90.0,
        })
    except Exception as exc:
        manifest.update(status="failed", completed_utc=utc_now(), failure_reason=str(exc))
        completed, rankable = count_outputs(output)
        atomic_json(progress, {
            "schema_version": "linkerforge.panel-progress.v1",
            "state": "failed",
            "updated_utc": utc_now(),
            "completed_constructs": completed,
            "total_constructs": total_constructs,
            "rankable_models": rankable,
            "planned_rankable_models": total_models,
            "prediction_percent": round(100 * rankable / total_models, 1),
            "overall_percent": round(90 * rankable / total_models, 1),
            "error": str(exc),
        })
        raise
    finally:
        validate_manifest(manifest)
        atomic_json(manifest_path, manifest)
        if lock.exists() and lock.read_text(encoding="utf-8") == lock_contents:
            lock.unlink()


if __name__ == "__main__":
    main()
