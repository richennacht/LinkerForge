"""Create a hash-pinned, label-free MMseqs2 input for B1 family audit.

Only observed-positive training candidates and sealed benchmark sequences are
exported. This is not a split or a training authorization. Output is ignored.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

if __package__:
    from scripts.audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
else:
    from audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha(path: Path) -> str:
    return digest(path.read_bytes())


def prepare(root: Path, output_dir: Path) -> dict:
    root, output_dir = root.resolve(), output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit input: {output_dir}")
    pool_root = root / "data/processed/linker_pred2_style_pool_v1"
    pool_manifest_path = pool_root / "manifest.json"
    pool_manifest = json.loads(pool_manifest_path.read_text(encoding="utf-8"))
    catalog_path = root / "data/processed/dfl_positive_candidates_2026_10_12.jsonl"
    catalog_manifest_path = root / "data/processed/dfl_positive_candidates_2026_10_12.manifest.json"
    catalog_manifest = json.loads(catalog_manifest_path.read_text(encoding="utf-8"))
    sources = []
    for name in ("published_dfl_train", "disprot_positive_only"):
        path = pool_root / f"{name}.jsonl"
        if sha(path) != pool_manifest["files"][name]["sha256"]:
            raise ValueError(f"Source hash mismatch: {name}")
        sources.append((name, path, pool_manifest["files"][name]["records"]))
    if sha(catalog_path) != catalog_manifest["output_sha256"]:
        raise ValueError("Candidate catalog hash mismatch")
    sources.append(("disprot_candidate_catalog", catalog_path, catalog_manifest["records"]))
    candidate_by_sequence = defaultdict(list)
    for name, path, count in sources:
        rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
        if len(rows) != count:
            raise ValueError(f"Source record count mismatch: {name}")
        for row in rows:
            sequence, labels = row["sequence"].upper(), row["labels"]
            if not sequence or len(sequence) != len(labels) or set(labels) - set("01x"):
                raise ValueError(f"Invalid source mask: {name}/{row['protein_id']}")
            if "1" in labels:
                candidate_by_sequence[sequence].append({"source": name, "protein_id": row["protein_id"]})
    caid_path = root / "artifacts/caid3-linker.fasta"
    if sha(caid_path) != pool_manifest["source_sha256"]["caid3"]:
        raise ValueError("CAID3 reference hash mismatch")
    caid3 = caid3_sequences(caid_path)
    dp_root = root / "data/processed/prospective_dfl_holdout_dp94_v1"
    dp_manifest_path = dp_root / "manifest.json"
    dp_manifest = json.loads(dp_manifest_path.read_text(encoding="utf-8"))
    if dp_manifest["status"] != "candidate_not_certified_independent":
        raise ValueError("DP94 registry status changed")
    dp_path = dp_root / "sequences.fasta"
    canonical = dp_path.read_text(encoding="utf-8").replace("\r\n", "\n")
    if digest(canonical.encode()) != dp_manifest["sequence_fasta_sha256"]:
        raise ValueError("DP94 FASTA hash mismatch")
    dp94 = two_line_fasta(dp_path)
    if [identifier for identifier, _ in dp94] != dp_manifest["holdout_ids"]:
        raise ValueError("DP94 ID roster mismatch")
    records = []
    for index, sequence in enumerate(sorted(candidate_by_sequence)):
        records.append({"record_id": f"C{index:04d}", "role": "candidate_positive_not_training_cleared",
                        "sequence": sequence, "sequence_sha256": digest(sequence.encode()),
                        "source_ids": sorted(candidate_by_sequence[sequence], key=lambda x: (x["source"], x["protein_id"]))})
    for registry, entries in (("CAID3", caid3), ("DP94", dp94)):
        for index, (identifier, sequence) in enumerate(entries):
            records.append({"record_id": f"{registry}_{index:03d}", "role": "sealed_benchmark_sequence_only",
                            "registry": registry, "source_id": identifier, "sequence": sequence,
                            "sequence_sha256": digest(sequence.encode())})
    output_dir.mkdir(parents=True)
    fasta_path = output_dir / "sequences.fasta"
    records_path = output_dir / "records.jsonl"
    fasta_path.write_text("".join(f">{row['record_id']}\n{row['sequence']}\n" for row in records), encoding="utf-8")
    records_path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in records), encoding="utf-8")
    manifest = {"schema_version": "1.0.0", "status": "sequence_only_family_audit_not_split",
                "candidate_sequences": len(candidate_by_sequence), "benchmark_sequences": len(caid3) + len(dp94),
                "input_hashes": {name: sha(path) for name, path, _ in sources},
                "caid3_reference_sha256": sha(caid_path),
                "dp94_sequence_manifest_sha256": sha(dp_manifest_path),
                "files_sha256": {fasta_path.name: sha(fasta_path), records_path.name: sha(records_path)},
                "label_seal": "benchmark_annotation_tracks_not_exported_or_inspected"}
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source_root, args.output_dir), indent=2, sort_keys=True))
