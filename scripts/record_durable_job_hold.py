#!/usr/bin/env python3
"""Record an intentional hold after a durable job has stopped."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from jsonschema import Draft202012Validator


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schemas" / "durable-job-hold.v1.schema.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-directory", type=Path, required=True)
    parser.add_argument("--reason", required=True)
    parser.add_argument("--resume-policy", required=True)
    args = parser.parse_args()

    run_directory = args.run_directory.resolve()
    request_path = run_directory / "request.json"
    status_path = run_directory / "status.json"
    hold_path = run_directory / "hold.json"
    if hold_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing hold record: {hold_path}")
    if not request_path.is_file() or not status_path.is_file():
        raise FileNotFoundError("Durable run request.json and status.json are required")
    status = json.loads(status_path.read_text(encoding="utf-8-sig"))
    if status.get("state") == "running":
        raise RuntimeError("Refusing to record a hold while status.json still reports running")

    record = {
        "schema_version": "linkerforge.durable-job-hold.v1",
        "state": "on_hold",
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "reason": args.reason,
        "resume_policy": args.resume_policy,
        "run_directory": str(run_directory),
        "request_sha256": sha256(request_path),
        "status_sha256": sha256(status_path),
        "stopped_state": status.get("state"),
        "stopped_stage": status.get("current_stage"),
        "claim_boundary": "An operational hold record; not a completed scientific result.",
    }
    Draft202012Validator(json.loads(SCHEMA.read_text(encoding="utf-8"))).validate(record)
    hold_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
