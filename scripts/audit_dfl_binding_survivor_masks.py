"""Audit DFL source-mask semantics for B1 capacity survivors, without test labels.

The input is the hash-pinned capacity manifest, not a training split. A
published zero is counted as unverified; a DisProt nonpositive remains x.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


SOURCE_PATHS = {
    "published_dfl_train": "data/processed/linker_pred2_style_pool_v1/published_dfl_train.jsonl",
    "disprot_positive_only": "data/processed/linker_pred2_style_pool_v1/disprot_positive_only.jsonl",
    "disprot_candidate_catalog": "data/processed/dfl_positive_candidates_2026_10_12.jsonl",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_rows(root: Path) -> tuple[dict, dict]:
    pool = json.loads((root / "data/processed/linker_pred2_style_pool_v1/manifest.json").read_text(encoding="utf-8"))
    catalog = json.loads((root / "data/processed/dfl_positive_candidates_2026_10_12.manifest.json").read_text(encoding="utf-8"))
    expected = {
        "published_dfl_train": pool["files"]["published_dfl_train"]["sha256"],
        "disprot_positive_only": pool["files"]["disprot_positive_only"]["sha256"],
        "disprot_candidate_catalog": catalog["output_sha256"],
    }
    by_key = {}
    for source, relative_path in SOURCE_PATHS.items():
        path = root / relative_path
        if sha(path) != expected[source]:
            raise ValueError(f"Source hash mismatch: {source}")
        for line in path.read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            key = (source, row["protein_id"])
            if key in by_key:
                raise ValueError(f"Duplicate source key: {key}")
            if len(row["sequence"]) != len(row["labels"]) or set(row["labels"]) - set("01x"):
                raise ValueError(f"Invalid mask: {key}")
            if source != "published_dfl_train" and "0" in row["labels"]:
                raise ValueError(f"Positive-only source contains zero: {key}")
            by_key[key] = row
    return by_key, expected


def audit(root: Path, capacity_path: Path, expected_capacity_hash: str) -> dict:
    root, capacity_path = root.resolve(), capacity_path.resolve()
    if sha(capacity_path) != expected_capacity_hash:
        raise ValueError("Capacity manifest hash mismatch")
    capacity = json.loads(capacity_path.read_text(encoding="utf-8"))
    if capacity.get("status") != "screened_capacity_not_training_split":
        raise ValueError("Unexpected capacity status")
    source, source_hashes = source_rows(root)
    counts = Counter()
    seen_sequences = set()
    source_counts = Counter()
    evidence_counts = Counter()
    for candidate in capacity["candidates"]:
        masks = []
        for identifier in candidate["source_ids"]:
            name, protein_id = identifier["source"], identifier["protein_id"]
            row = source[(name, protein_id)]
            sequence_hash = hashlib.sha256(row["sequence"].upper().encode()).hexdigest()
            if sequence_hash != candidate["sequence_sha256"]:
                raise ValueError(f"Source sequence changed: {name}/{protein_id}")
            masks.append(row["labels"])
            source_counts[name] += 1
            if name == "disprot_candidate_catalog":
                for region in row["positive_regions"]:
                    evidence_counts[region.get("eco_id") or "missing_eco_id"] += 1
                    counts["catalog_regions_with_reference_id"] += bool(region.get("reference_id"))
        if candidate["sequence_sha256"] in seen_sequences:
            raise ValueError("Duplicate survivor sequence")
        seen_sequences.add(candidate["sequence_sha256"])
        if not masks or not any("1" in mask for mask in masks):
            raise ValueError("Survivor has no source-observed positive")
        counts["survivor_sequences"] += 1
        counts["observed_positive_residues_union"] += sum(any(mask[i] == "1" for mask in masks)
                                                           for i in range(len(masks[0])))
        counts["published_unverified_zero_residues"] += sum(mask.count("0") for mask in masks)
        counts["cross_source_positive_zero_conflicts"] += sum(
            any(mask[i] == "1" for mask in masks) and any(mask[i] == "0" for mask in masks)
            for i in range(len(masks[0]))
        )
        counts["multiple_source_sequences"] += len(masks) > 1
    return {"status": "binding_survivor_mask_audit_not_training_authorization",
            "capacity_manifest_sha256": sha(capacity_path), "source_sha256": source_hashes,
            "counts": dict(counts), "source_record_appearances": dict(source_counts),
            "catalog_eco_region_counts": dict(evidence_counts),
            "limitations": ["Published zeros remain unverified negatives.",
                            "DisProt unknown residues remain unknown, never zero.",
                            "ECO/region source metadata were read only for candidate training records, never benchmark labels.",
                            "No train/development assignment or independent no-linker set is created here."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--capacity-manifest", type=Path, required=True)
    parser.add_argument("--capacity-manifest-sha256", required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.source_root, args.capacity_manifest,
                           args.capacity_manifest_sha256), indent=2, sort_keys=True))
