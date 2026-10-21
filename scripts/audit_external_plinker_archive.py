#!/usr/bin/env python3
"""Audit the P-Linker ZIP as untrusted data without extracting or executing it."""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
from pathlib import Path, PurePosixPath
from zipfile import ZipFile


SILVER_ENTRY = "Linker - Copy/data/processed/dataset_v2.jsonl"
GOLD_ENTRY = "Linker - Copy/data/processed/dataset_v1.jsonl"
RAW_ENTRY = "Linker - Copy/data/processed/mined_linkers.jsonl"
DATA_CARD_ENTRY = "Linker - Copy/docs/data_card.md"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def jsonl_from_zip(archive: ZipFile, entry: str) -> tuple[bytes, list[dict]]:
    payload = archive.read(entry)
    rows = [json.loads(line) for line in payload.decode("utf-8").splitlines() if line.strip()]
    return payload, rows


def unsafe_member_names(names: list[str]) -> list[str]:
    unsafe: list[str] = []
    for name in names:
        normalized = name.replace("\\", "/")
        path = PurePosixPath(normalized)
        if path.is_absolute() or ".." in path.parts or (path.parts and ":" in path.parts[0]):
            unsafe.append(name)
    return unsafe


def build_audit(zip_path: Path, existing_corpus_path: Path) -> dict:
    with ZipFile(zip_path) as archive:
        infos = archive.infolist()
        names = [info.filename for info in infos]
        required = [SILVER_ENTRY, GOLD_ENTRY, RAW_ENTRY, DATA_CARD_ENTRY]
        missing = [entry for entry in required if entry not in names]
        if missing:
            raise ValueError(f"Archive is missing required audit entries: {missing}")
        silver_bytes, silver = jsonl_from_zip(archive, SILVER_ENTRY)
        gold_bytes, gold = jsonl_from_zip(archive, GOLD_ENTRY)
        raw_bytes, raw = jsonl_from_zip(archive, RAW_ENTRY)
        data_card = archive.read(DATA_CARD_ENTRY).decode("utf-8")

    existing = [
        json.loads(line)
        for line in existing_corpus_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    existing_proteins = {str(row.get("protein_id", "")) for row in existing}
    existing_linkers = {str(row.get("linker_sequence", "")) for row in existing}

    silver_accessions = {row["protein_accession"] for row in silver}
    silver_linkers = {row["sequences"]["linker"] for row in silver}
    domain_pairs = collections.Counter(
        (row["domain_A_acc"], row["domain_B_acc"]) for row in silver
    )
    exact_coordinate_errors = 0
    full_length_errors = 0
    ambiguous_full_sequence_rows = 0
    ambiguous_linker_rows = 0
    for row in silver:
        sequences = row["sequences"]
        full = sequences["full"]
        if len(full) != int(row["protein_length"]):
            full_length_errors += 1
        coordinate_ok = (
            int(row["linker_length"]) == int(row["linker_end"]) - int(row["linker_start"]) + 1
            and sequences["domain_A"] == full[int(row["domain_A_start"]) - 1 : int(row["domain_A_end"])]
            and sequences["linker"] == full[int(row["linker_start"]) - 1 : int(row["linker_end"])]
            and sequences["domain_B"] == full[int(row["domain_B_start"]) - 1 : int(row["domain_B_end"])]
        )
        if not coordinate_ok:
            exact_coordinate_errors += 1
        ambiguous_full_sequence_rows += bool(set(full) - STANDARD_AA)
        ambiguous_linker_rows += bool(set(sequences["linker"]) - STANDARD_AA)

    gold_boundary_mismatches = []
    for row in gold:
        domain_a = row.get("domains", {}).get("domain_A", {})
        domain_b = row.get("domains", {}).get("domain_B", {})
        linker = row.get("linker", {})
        consistent = (
            isinstance(domain_a.get("end"), int)
            and isinstance(domain_b.get("start"), int)
            and isinstance(linker.get("start"), int)
            and isinstance(linker.get("end"), int)
            and linker["start"] == domain_a["end"] + 1
            and domain_b["start"] == linker["end"] + 1
        )
        if not consistent:
            gold_boundary_mismatches.append(
                {
                    "pdb_id": row.get("pdb_id"),
                    "declared_domain_a_end": domain_a.get("end"),
                    "found_linker_start": linker.get("start"),
                    "found_linker_end": linker.get("end"),
                    "declared_domain_b_start": domain_b.get("start"),
                }
            )

    same_family = sum(row["domain_A_acc"] == row["domain_B_acc"] for row in silver)
    seven_tm_ids = {"PF00001", "PF00002", "PF00003"}
    seven_tm_rows = sum(
        row["domain_A_acc"] in seven_tm_ids or row["domain_B_acc"] in seven_tm_ids
        for row in silver
    )
    source_metadata_rows = sum(
        any(key in row for key in ("source_ids", "retrieved_at", "release", "license", "citation"))
        for row in silver
    )
    project_license_names = {
        "linker - copy/license",
        "linker - copy/license.md",
        "linker - copy/license.txt",
        "linker - copy/copying",
    }
    has_license = any(name.rstrip("/").lower() in project_license_names for name in names)

    return {
        "schema_version": "linkerforge.external-linker-dataset-audit.v1",
        "status": "quarantined_not_training_authorized",
        "source_artifact": {
            "filename": zip_path.name,
            "sha256": sha256_file(zip_path),
            "entry_count": len(infos),
            "uncompressed_bytes": sum(info.file_size for info in infos),
            "unsafe_member_names": unsafe_member_names(names),
            "embedded_virtual_environment_entries": sum("/.venv/" in name for name in names),
            "embedded_git_entries": sum("/.git/" in name for name in names),
            "license_file_present": has_license,
        },
        "files": {
            "silver": {"entry": SILVER_ENTRY, "sha256": sha256_bytes(silver_bytes), "records": len(silver)},
            "gold": {"entry": GOLD_ENTRY, "sha256": sha256_bytes(gold_bytes), "records": len(gold)},
            "raw": {"entry": RAW_ENTRY, "sha256": sha256_bytes(raw_bytes), "records": len(raw)},
        },
        "data_card_consistency": {
            "claims_7903_total": "7,903" in data_card,
            "claims_11_gold": "Gold | 11" in data_card,
            "observed_gold_records": len(gold),
            "observed_silver_records": len(silver),
            "observed_total_gold_plus_silver": len(gold) + len(silver),
            "consistent": len(gold) == 11 and len(gold) + len(silver) == 7903,
        },
        "silver_profile": {
            "unique_accessions": len(silver_accessions),
            "unique_linker_sequences": len(silver_linkers),
            "linker_length_min": min(int(row["linker_length"]) for row in silver),
            "linker_length_max": max(int(row["linker_length"]) for row in silver),
            "within_linkerforge_5_to_30_residues": sum(5 <= int(row["linker_length"]) <= 30 for row in silver),
            "outside_linkerforge_5_to_30_residues": sum(not 5 <= int(row["linker_length"]) <= 30 for row in silver),
            "same_family_domain_pairs": same_family,
            "cross_family_domain_pairs": len(silver) - same_family,
            "seven_transmembrane_family_rows": seven_tm_rows,
            "ambiguous_full_sequence_rows": ambiguous_full_sequence_rows,
            "ambiguous_linker_rows": ambiguous_linker_rows,
            "coordinate_or_slice_errors": exact_coordinate_errors,
            "full_length_errors": full_length_errors,
            "rows_with_record_level_source_or_license_metadata": source_metadata_rows,
            "domain_pair_counts": [
                {"domain_a": pair[0], "domain_b": pair[1], "records": count}
                for pair, count in domain_pairs.most_common()
            ],
        },
        "gold_profile": {
            "records": len(gold),
            "declared_boundary_mismatch_records": len(gold_boundary_mismatches),
            "boundary_mismatches": gold_boundary_mismatches,
        },
        "overlap_with_existing_27766_corpus": {
            "exact_protein_identifier_overlap": len(silver_accessions & existing_proteins),
            "unique_linker_sequence_overlap": len(silver_linkers & existing_linkers),
        },
        "fitness": {
            "compatibility_ranker_training": False,
            "curator_validated_linker_corpus": False,
            "automatic_merge_into_existing_corpus": False,
            "possible_use": "quarantined mining-method audit and reconstruction source after licensing, provenance, topology, homology, and boundary review",
        },
        "blockers": [
            "No fusion compatibility outcomes, matched negatives, or success/failure labels",
            "No homology-aware split or family-disjoint evaluation contract",
            "Nearly all silver records join two copies of the same Pfam family",
            "Seven-transmembrane families dominate and conflict with the frozen soluble-domain scope",
            "Gold record count disagrees with the data card and declared boundaries disagree with located linker motifs",
            "No archive licence file and no record-level source release, retrieval date, or licence fields",
        ],
        "claim_boundary": (
            "This archive audit does not validate linker identity, compatibility, or biological function. "
            "Archived code and environments were not executed."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", type=Path, required=True)
    parser.add_argument("--existing-corpus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output already exists; refusing overwrite: {args.output}")
    result = build_audit(args.zip, args.existing_corpus)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
