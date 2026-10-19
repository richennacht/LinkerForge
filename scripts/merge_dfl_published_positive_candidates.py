"""Exact-sequence deduplicate pinned DFL sources under the 24-sequence capacity screen.

All source zeros are dropped unless independently adjudicated; this merged
positive/unknown candidate pool is not a B1 training split. Benchmark labels
are never opened. The pinned capacity manifest already encodes the binding
sequence and exploratory family screens, which this builder checks again.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from jsonschema import Draft202012Validator

if __package__:
    from scripts.audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
    from scripts.preflight_dfl_family_training import local_hit
else:
    from audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
    from preflight_dfl_family_training import local_hit


CAPACITY_SHA256 = "bca74cae3932d8b086d5d292d18adb41f0c9ce39e9744cc4b05d63fea9912368"
SOURCE_SHA256 = {
    "published_dfl_train": "a1a5f3b3b30548ad65f6f378535c7c732272b8c9b60048b38f6500fca9bb1aed",
    "disprot_positive_only": "0265b062ca9ff88cde2576789b613860355f33842fa9714788ba1105a19df42a",
    "disprot_candidate_catalog": "2bc6726545f47256bcec70376bd51d689db0b5c80a64f45a221b9d88bb0c4c8d",
    "disprot_curated_strict": "07e054bbf8d3be9d17438cd843d0610bda84cdd7baa0d1e5b6083e4ae2b98845",
}
CAID_SHA256 = "4c204ef49ffa12fe166883e6b2cff0012cabee21260b8f32f997081d0e2fa0f8"
DP94_MANIFEST_SHA256 = "fe64e7c688cfac48ed0dac73bcf09daae5cb0adfd0e03e9d530824d3d92c0051"
CONFIG = {"schema_version": "1.0.0", "merge_key": "sha256_uppercase_full_sequence",
          "positive_policy": "union_only_explicit_source_ones_except_one_zero_disagreement_masked_x",
          "zero_policy": "discard_all_unverified_source_zeros_as_unknown",
          "eligibility": "capacity_manifest_24_sequence_allowlist_plus_direct_registry_rescreen"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha(value: dict) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def source_paths(root: Path) -> dict[str, Path]:
    pool = root / "data/processed/linker_pred2_style_pool_v1"
    return {"published_dfl_train": pool / "published_dfl_train.jsonl",
            "disprot_positive_only": pool / "disprot_positive_only.jsonl",
            "disprot_candidate_catalog": root / "data/processed/dfl_positive_candidates_2026_10_12.jsonl",
            "disprot_curated_strict": root / "data/processed/dfl_curated_evidence_pool_v6/screened.jsonl"}


def union_positive_masks(sequence: str, source_labels: list[str]) -> str:
    if not source_labels:
        raise ValueError("No observations to merge")
    if any(len(labels) != len(sequence) or set(labels) - set("01x")
           for labels in source_labels):
        raise ValueError("Incompatible source mask")
    return "".join("1" if (any(labels[i] == "1" for labels in source_labels)
                            and not any(labels[i] == "0" for labels in source_labels)) else "x"
                   for i in range(len(sequence)))


def merge(root: Path, capacity_path: Path, output_dir: Path) -> dict:
    root, capacity_path, output_dir = root.resolve(), capacity_path.resolve(), output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite merged candidate pool: {output_dir}")
    if sha(capacity_path) != CAPACITY_SHA256:
        raise ValueError("Pinned family/binding capacity manifest changed")
    capacity = json.loads(capacity_path.read_text(encoding="utf-8"))
    if capacity.get("status") != "screened_capacity_not_training_split":
        raise ValueError("Capacity manifest status changed")
    allow = {row["sequence_sha256"]: row["component_id"] for row in capacity["candidates"]}
    if len(allow) != 24 or len(capacity["candidates"]) != 24:
        raise ValueError("Expected 24 unique sequence-only survivor hashes")
    caid_path = root / "artifacts/caid3-linker.fasta"
    dp_manifest_path = root / "data/processed/prospective_dfl_holdout_dp94_v1/manifest.json"
    if sha(caid_path) != CAID_SHA256 or sha(dp_manifest_path) != DP94_MANIFEST_SHA256:
        raise ValueError("Protected sequence registry changed")
    dp_root = dp_manifest_path.parent
    dp_manifest = json.loads(dp_manifest_path.read_text(encoding="utf-8"))
    if dp_manifest.get("status") != "candidate_not_certified_independent":
        raise ValueError("DP94 registry status changed")
    dp_fasta = dp_root / "sequences.fasta"
    if hashlib.sha256(dp_fasta.read_text(encoding="utf-8").replace("\r\n", "\n").encode()).hexdigest() != dp_manifest["sequence_fasta_sha256"]:
        raise ValueError("DP94 registry sequence hash changed")
    references = caid3_sequences(caid_path) + two_line_fasta(dp_fasta)
    paths = source_paths(root)
    for name, path in paths.items():
        if sha(path) != SOURCE_SHA256[name]:
            raise ValueError(f"Source hash changed: {name}")
    schema_path = Path(__file__).resolve().parents[1] / "schemas/dfl_merged_positive_candidate.v1.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)
    observations = defaultdict(list)
    sequences = {}
    counts = Counter()
    for name, path in paths.items():
        for line in path.read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            sequence = row["sequence"].upper()
            labels = row["labels"]
            if len(sequence) != len(labels) or set(labels) - set("01x"):
                raise ValueError(f"Invalid source mask: {name}/{row.get('protein_id')}")
            if name == "disprot_candidate_catalog":
                if row.get("record_type") != "disprot_dfl_positive_candidate" or "0" in labels:
                    raise ValueError("Catalog is not positive/unknown DFL evidence")
            elif name == "disprot_curated_strict":
                if (row.get("record_type") != "disprot_dfl_evidence_candidate"
                        or row.get("training_eligibility") != "candidate_not_B1"):
                    raise ValueError("Strict curated source status changed")
            elif row.get("task") not in {"DFL", "DFL_positive_unlabeled"}:
                raise ValueError("Substituted task in published/DisProt source")
            sequence_sha = hashlib.sha256(sequence.encode()).hexdigest()
            if sequence_sha not in allow or "1" not in labels:
                continue
            protein_id = row["protein_id"]
            if any(protein_id == ref_id or local_hit(sequence, ref_sequence)
                   for ref_id, ref_sequence in references):
                raise ValueError(f"Capacity allowlist contradicted by registry: {name}/{protein_id}")
            if sequence_sha in sequences and sequences[sequence_sha] != sequence:
                raise ValueError("Sequence hash collision")
            sequences[sequence_sha] = sequence
            observations[sequence_sha].append((name, protein_id, labels))
            counts[f"{name}_retained_source_rows"] += 1
    if set(observations) != set(allow):
        raise ValueError(f"Capacity/source mismatch: missing={len(set(allow)-set(observations))}")
    config_hash = canonical_sha(CONFIG)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    rows = []
    for sequence_sha in sorted(allow):
        sequence = sequences[sequence_sha]
        source_rows = []
        for name, protein_id, labels in sorted(observations[sequence_sha]):
            source_rows.append({"source": name, "source_id": protein_id,
                                "source_file_sha256": SOURCE_SHA256[name],
                                "positive_residues": labels.count("1"),
                                "discarded_unverified_zero_residues": labels.count("0")})
        labels = union_positive_masks(sequence, [labels for _, _, labels in observations[sequence_sha]])
        disagreement = sum(any(src_labels[i] == "1" for _, _, src_labels in observations[sequence_sha])
                           and any(src_labels[i] == "0" for _, _, src_labels in observations[sequence_sha])
                           for i in range(len(sequence)))
        row = {"schema_version": "1.0.0", "record_type": "deduplicated_dfl_positive_candidate",
               "sequence_sha256": sequence_sha, "sequence": sequence, "labels": labels,
               "masked_disagreement_residues": disagreement,
               "family_component_id": allow[sequence_sha], "source_observations": source_rows,
               "provenance": {"capacity_manifest_sha256": CAPACITY_SHA256,
                              "config_hash": config_hash, "schema_sha256": sha(schema_path),
                              "pipeline_commit": commit},
               "training_eligibility": "candidate_not_B1_no_adjudicated_negatives_or_frozen_split"}
        validator.validate(row)
        rows.append(row)
        counts["merged_positive_residues"] += labels.count("1")
        counts["merged_unknown_residues"] += labels.count("x")
        counts["masked_disagreement_residues"] += disagreement
        counts["discarded_unverified_zero_residue_observations"] += sum(
            source["discarded_unverified_zero_residues"] for source in source_rows)
    counts["unique_positive_sequences"] = len(rows)
    counts["duplicate_source_rows_removed"] = sum(counts[f"{name}_retained_source_rows"]
                                                  for name in SOURCE_SHA256) - len(rows)
    output_dir.mkdir(parents=True)
    records_path = output_dir / "positive_unknown_candidates.jsonl"
    records_path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
                            encoding="utf-8", newline="\n")
    manifest = {"schema_version": "1.0.0", "status": "deduplicated_positive_unknown_candidates_not_B1",
                "config": CONFIG, "config_hash": config_hash,
                "source_sha256": SOURCE_SHA256, "capacity_manifest_sha256": CAPACITY_SHA256,
                "record_schema_sha256": sha(schema_path),
                "benchmark_sequence_only": {"caid3_sha256": CAID_SHA256,
                                            "dp94_manifest_sha256": DP94_MANIFEST_SHA256},
                "records_sha256": sha(records_path), "counts": dict(counts),
                "limitations": ["No verified whole-protein no-linker negatives or calibrated precision set.",
                                "Family graph is exploratory; no development allocation or B1 permission.",
                                "Source zeros were discarded because their negative semantics are not verified.",
                                "Only source papers/datasets with pinned local row-level evidence are represented."]}
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8", newline="\n")
    return {"manifest": str(manifest_path), "manifest_sha256": sha(manifest_path),
            "counts": dict(counts)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--capacity-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(merge(args.source_root, args.capacity_manifest, args.output_dir),
                     indent=2, sort_keys=True))
