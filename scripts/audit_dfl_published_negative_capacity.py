"""Screen published non-DFL training records without treating unknown as zero.

This counts a potential negative source, not adjudicated no-linker proteins.
The author-provided training masks and DisProt family metadata are source
data. Benchmark annotations are never opened or parsed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

if __package__:
    from scripts.audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
    from scripts.audit_dfl_disprot_family_metadata import DISPROT_SHA256, project_entry
    from scripts.preflight_dfl_family_training import local_hit
else:
    from audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
    from audit_dfl_disprot_family_metadata import DISPROT_SHA256, project_entry
    from preflight_dfl_family_training import local_hit


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def is_nonpositive_published(row: dict) -> bool:
    sequence, labels = row["sequence"].upper(), row["labels"]
    if (row.get("source") != "TransDFL_TR166" or row.get("task") != "DFL"
            or not sequence or len(sequence) != len(labels) or set(labels) - set("01x")):
        raise ValueError(f"Invalid published DFL record: {row.get('protein_id')}")
    return "1" not in labels and "0" in labels


def audit(root: Path) -> dict:
    root = root.resolve()
    pool = root / "data/processed/linker_pred2_style_pool_v1"
    pool_manifest = json.loads((pool / "manifest.json").read_text(encoding="utf-8"))
    source_path = pool / "published_dfl_train.jsonl"
    if sha(source_path) != pool_manifest["files"]["published_dfl_train"]["sha256"]:
        raise ValueError("Published DFL train source changed")
    records = [json.loads(line) for line in source_path.read_text(encoding="utf-8").splitlines()]
    if len(records) != pool_manifest["files"]["published_dfl_train"]["records"]:
        raise ValueError("Published DFL train record count changed")
    caid_path = root / "artifacts/caid3-linker.fasta"
    if sha(caid_path) != pool_manifest["source_sha256"]["caid3"]:
        raise ValueError("CAID3 source changed")
    caid3 = caid3_sequences(caid_path)
    dp_root = root / "data/processed/prospective_dfl_holdout_dp94_v1"
    dp_manifest_path = dp_root / "manifest.json"
    dp_manifest = json.loads(dp_manifest_path.read_text(encoding="utf-8"))
    if dp_manifest.get("status") != "candidate_not_certified_independent":
        raise ValueError("DP94 registry status changed")
    dp_fasta = dp_root / "sequences.fasta"
    canonical = dp_fasta.read_text(encoding="utf-8").replace("\r\n", "\n")
    if hashlib.sha256(canonical.encode()).hexdigest() != dp_manifest["sequence_fasta_sha256"]:
        raise ValueError("DP94 sequence roster changed")
    dp94 = two_line_fasta(dp_fasta)
    if [key for key, _ in dp94] != dp_manifest["holdout_ids"]:
        raise ValueError("DP94 IDs changed")
    raw_path = root / "data/raw/disprot/disprot_2026_06_regions.json"
    if sha(raw_path) != DISPROT_SHA256:
        raise ValueError("DisProt family source changed")
    candidates = [row for row in records if is_nonpositive_published(row)]
    candidate_ids = {row["protein_id"] for row in candidates}
    benchmark_ids = {key for key, _ in caid3 + dp94}
    needed_ids = candidate_ids | benchmark_ids
    by_id = {}
    for entry in json.loads(raw_path.read_text(encoding="utf-8"))["data"]:
        if entry["disprot_id"] in needed_ids:
            try:
                by_id[entry["disprot_id"]] = project_entry(entry)
            except ValueError:
                continue
    benchmark_families = {name: set() for name in ("pfam", "gene3D", "uniref50")}
    for key, sequence in caid3 + dp94:
        entry = by_id.get(key)
        if entry and entry["sequence"] == sequence:
            for name, values in entry["families"].items():
                benchmark_families[name].update(values)
    counts = Counter()
    for row in candidates:
        key, sequence, labels = row["protein_id"], row["sequence"].upper(), row["labels"]
        counts["nonpositive_published_records"] += 1
        counts["published_zero_residues_unverified"] += labels.count("0")
        counts["unknown_residues"] += labels.count("x")
        metadata = by_id.get(key)
        exact_metadata = metadata is not None and metadata["sequence"] == sequence
        counts["exact_disprot_family_metadata"] += exact_metadata
        if key in benchmark_ids or any(local_hit(sequence, ref) for _, ref in caid3 + dp94):
            counts["binding_direct_hit_either_registry"] += 1
            continue
        counts["binding_direct_survivors"] += 1
        if not exact_metadata:
            counts["binding_survivors_without_exact_family_metadata"] += 1
            continue
        if any(set(metadata["families"][name]) & benchmark_families[name]
               for name in benchmark_families):
            counts["binding_survivors_sharing_direct_family_id"] += 1
        else:
            counts["binding_survivors_without_direct_family_id_hit"] += 1
    return {"status": "published_negative_capacity_not_adjudicated_or_split",
            "source_sha256": sha(source_path), "disprot_family_source_sha256": sha(raw_path),
            "caid3_reference_sha256": sha(caid_path),
            "dp94_sequence_manifest_sha256": sha(dp_manifest_path), "counts": dict(counts),
            "limitations": ["Published zeros are not independently verified per-residue negatives.",
                            "Direct family IDs omit transitive and remote-homology links.",
                            "Benchmark labels are never parsed. No training or calibration is authorized."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.source_root), indent=2, sort_keys=True))
