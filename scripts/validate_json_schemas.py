#!/usr/bin/env python3
"""Fail when a checked-in JSON Schema is malformed or uses an unknown dialect."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from jsonschema.validators import validator_for


def validate_schema(path: Path) -> None:
    schema = json.loads(path.read_text(encoding="utf-8"))
    validator = validator_for(schema)
    validator.check_schema(schema)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--schemas", type=Path, default=Path("schemas"))
    args = parser.parse_args()
    paths = sorted(args.schemas.glob("*.json"))
    if not paths:
        raise SystemExit(f"No JSON schemas found in {args.schemas}")
    for path in paths:
        validate_schema(path)
        print(f"valid {path}")
    print(f"validated {len(paths)} schemas")


if __name__ == "__main__":
    main()
