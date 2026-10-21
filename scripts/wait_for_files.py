#!/usr/bin/env python3
"""Wait for durable upstream completion markers without relying on a chat session."""

from __future__ import annotations

import argparse
import time
from datetime import UTC, datetime
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", action="append", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=int, default=14_400)
    parser.add_argument("--poll-seconds", type=int, default=15)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.timeout_seconds < 1 or args.poll_seconds < 1:
        raise ValueError("Timeout and poll interval must be positive")
    targets = [path.resolve() for path in args.path]
    deadline = time.monotonic() + args.timeout_seconds
    last_missing: tuple[str, ...] | None = None
    while True:
        missing = tuple(str(path) for path in targets if not path.is_file())
        if not missing:
            print(
                f"{datetime.now(UTC).isoformat()} all {len(targets)} completion markers exist",
                flush=True,
            )
            return 0
        if missing != last_missing:
            print(
                f"{datetime.now(UTC).isoformat()} waiting for: {', '.join(missing)}",
                flush=True,
            )
            last_missing = missing
        if time.monotonic() >= deadline:
            raise TimeoutError(f"Timed out waiting for: {', '.join(missing)}")
        time.sleep(min(args.poll_seconds, max(0.0, deadline - time.monotonic())))


if __name__ == "__main__":
    raise SystemExit(main())
