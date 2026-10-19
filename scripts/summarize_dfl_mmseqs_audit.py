"""Summarize sequence-only MMseqs2 clusters against sealed registries.

Requires an explicit expected TSV hash and verifies every exported record is
present exactly once. This is not a validated domain-family split.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(input_dir: Path, cluster_tsv: Path, expected_tsv_sha256: str) -> dict:
    input_dir, cluster_tsv = input_dir.resolve(), cluster_tsv.resolve()
    manifest_path = input_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["status"] != "sequence_only_family_audit_not_split":
        raise ValueError("Unexpected source manifest status")
    records_path = input_dir / "records.jsonl"
    fasta_path = input_dir / "sequences.fasta"
    for path in (records_path, fasta_path):
        if sha(path) != manifest["files_sha256"][path.name]:
            raise ValueError(f"Input hash mismatch: {path}")
    if sha(cluster_tsv) != expected_tsv_sha256:
        raise ValueError("MMseqs cluster TSV hash mismatch")
    records = {row["record_id"]: row for row in
               map(json.loads, records_path.read_text(encoding="utf-8").splitlines())}
    if len(records) != manifest["candidate_sequences"] + manifest["benchmark_sequences"]:
        raise ValueError("Input record count mismatch")
    clusters = defaultdict(list)
    seen = set()
    for line in cluster_tsv.read_text(encoding="utf-8").splitlines():
        rep, member = line.split("\t")
        if member not in records or rep not in records or member in seen:
            raise ValueError(f"Unexpected or duplicate cluster member: {member}")
        clusters[rep].append(member)
        seen.add(member)
    if seen != set(records):
        raise ValueError("MMseqs cluster output omits exported records")
    counts = Counter()
    counts["clusters"] = len(clusters)
    for members in clusters.values():
        benchmark = {records[key].get("registry") for key in members
                     if records[key]["role"] == "sealed_benchmark_sequence_only"}
        benchmark.discard(None)
        candidates = [key for key in members if records[key]["role"] == "candidate_positive_not_training_cleared"]
        counts["benchmark_touching_clusters"] += bool(benchmark)
        counts["candidate_only_clusters"] += bool(candidates) and not benchmark
        counts["candidate_sequences_in_benchmark_clusters"] += len(candidates) if benchmark else 0
        counts["candidate_sequences_outside_benchmark_clusters"] += len(candidates) if not benchmark else 0
        if candidates and benchmark:
            for registry in benchmark:
                counts[f"candidate_sequences_in_{registry}_clusters"] += len(candidates)
    return {"status": "sequence_only_mmseqs_sensitivity_not_B1",
            "input_manifest_sha256": sha(manifest_path), "cluster_tsv_sha256": sha(cluster_tsv),
            "counts": dict(counts),
            "limitations": ["MMseqs full-length clusters alone do not establish domain-family independence.",
                            "The binding local-similarity preflight is unchanged and still fails.",
                            "No benchmark labels or scores were read."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--cluster-tsv", type=Path, required=True)
    parser.add_argument("--cluster-tsv-sha256", required=True)
    args = parser.parse_args()
    print(json.dumps(summarize(args.input_dir, args.cluster_tsv, args.cluster_tsv_sha256), indent=2, sort_keys=True))
