#!/usr/bin/env python3
"""Release the known matched-panel GPU lock after its completion marker exists."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lock", required=True, type=Path)
    parser.add_argument("--completion-marker", required=True, type=Path)
    parser.add_argument("--expected-owner", required=True)
    parser.add_argument("--expected-output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    marker = args.completion_marker.resolve()
    if not marker.is_file():
        raise FileNotFoundError(f"Completion marker is missing: {marker}")

    lock = args.lock.resolve()
    if not lock.is_file():
        print(f"GPU lock already absent: {lock}")
        return 0
    fields: dict[str, str] = {}
    for line in lock.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            fields[key.strip()] = value.strip()
    if fields.get("owner") != args.expected_owner:
        raise ValueError(f"Unexpected GPU-lock owner: {fields.get('owner')!r}")
    expected_output = str(args.expected_output.resolve())
    actual_output = str(Path(fields.get("output", "")).resolve())
    if actual_output.casefold() != expected_output.casefold():
        raise ValueError(
            f"GPU-lock output mismatch: expected {expected_output!r}, got {actual_output!r}"
        )
    lock.unlink()
    print(f"Released completed GPU lock: {lock}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
