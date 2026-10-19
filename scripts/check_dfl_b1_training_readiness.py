"""Fail-closed B1 census for the deduplicated DFL positive/unknown pool.

This does not open benchmark labels or authorize training. It makes the
observed-label and family-split requirements visible before allocating GPU.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


EXPECTED_STATUS = "deduplicated_positive_unknown_candidates_not_B1"
EXPECTED_RECORD_TYPE = "deduplicated_dfl_positive_candidate"
EXPECTED_ELIGIBILITY = "candidate_not_B1_no_adjudicated_negatives_or_frozen_split"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(pool: Path) -> dict:
    pool = pool.resolve()
    manifest_path = pool / "manifest.json"
    records_path = pool / "positive_unknown_candidates.jsonl"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != EXPECTED_STATUS:
        raise ValueError("Unexpected pool status; this audit accepts only the candidate pool")
    if sha(records_path) != manifest.get("records_sha256"):
        raise ValueError("Candidate records hash mismatch")
    rows = [json.loads(line) for line in records_path.read_text(encoding="utf-8").splitlines() if line]
    if len(rows) != manifest.get("counts", {}).get("unique_positive_sequences") or not rows:
        raise ValueError("Candidate records count mismatch or empty pool")
    counts = Counter()
    sequence_hashes = set()
    component_ids = set()
    for row in rows:
        sequence, labels = row["sequence"], row["labels"]
        sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
        if (row.get("record_type") != EXPECTED_RECORD_TYPE
                or row.get("training_eligibility") != EXPECTED_ELIGIBILITY
                or row.get("sequence_sha256") != sequence_hash
                or sequence_hash in sequence_hashes
                or not sequence or len(sequence) != len(labels)
                or set(labels) - set("1x") or "1" not in labels
                or not row.get("source_observations")
                or not row.get("family_component_id")):
            raise ValueError("Invalid or duplicate candidate record")
        sequence_hashes.add(sequence_hash)
        component_ids.add(row["family_component_id"])
        counts.update(labels)
    if (counts["1"] != manifest["counts"].get("merged_positive_residues")
            or counts["x"] != manifest["counts"].get("merged_unknown_residues")):
        raise ValueError("Manifest mask counts mismatch")
    gates = {
        "verified_negative_residues": counts["0"] > 0,
        "frozen_family_disjoint_train_development_split": False,
        "independently_adjudicated_no_linker_panel": False,
        "benchmark_firewall_preflight_for_split": False,
    }
    # A candidate-only pool cannot satisfy split-level or panel-level gates.
    # These cannot be inferred from component IDs or unlabelled residues.
    return {
        "status": "ready" if all(gates.values()) else "blocked_not_training_ready",
        "manifest_sha256": sha(manifest_path),
        "records_sha256": sha(records_path),
        "candidate_sequences": len(rows),
        "exploratory_component_ids": len(component_ids),
        "residue_counts": {label: counts[label] for label in "10x"},
        "gates": gates,
        "reason": "The candidate pool has no verified negatives, frozen family split, or no-linker calibration panel; unknown is not negative.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pool", type=Path)
    args = parser.parse_args()
    try:
        result = audit(args.pool)
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        result = {"status": "blocked_invalid_pool", "reason": str(error)}
    print(json.dumps(result, indent=2, sort_keys=True))
    raise SystemExit(0 if result["status"] == "ready" else 2)
