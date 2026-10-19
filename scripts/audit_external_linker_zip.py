"""Read-only audit of Ram's P-Linker archive; never promotes gap candidates to DFL.

The archive is streamed directly without extraction or executing bundled code.
Only compact aggregate statistics are printed; protein sequences are not copied.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import zipfile

PREFIX = "Linker - Copy/data/processed/"
FILES = {
    "gold": PREFIX + "dataset_v1.jsonl",
    "gold_review": PREFIX + "dataset_v1_needs_review.jsonl",
    "silver": PREFIX + "dataset_v2.jsonl",
    "raw_candidates": PREFIX + "mined_linkers.jsonl",
}


def read_jsonl(archive: zipfile.ZipFile, name: str) -> tuple[list[dict], str]:
    digest = hashlib.sha256()
    rows: list[dict] = []
    with archive.open(name) as stream:
        for number, raw in enumerate(stream, 1):
            digest.update(raw)
            try:
                rows.append(json.loads(raw))
            except (UnicodeDecodeError, json.JSONDecodeError) as error:
                raise ValueError(f"Invalid JSON at {name}:{number}") from error
    return rows, digest.hexdigest()


def audit(archive_path: Path, split_manifest: Path | None = None) -> dict:
    with archive_path.open("rb") as archive_stream:
        archive_sha = hashlib.file_digest(archive_stream, "sha256").hexdigest()
    with zipfile.ZipFile(archive_path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate ZIP member names")
        for name in names:
            parts = Path(name.replace("\\", "/")).parts
            if name.startswith("/") or ".." in parts or ":" in name:
                raise ValueError("Unsafe ZIP member name")
        rows, hashes = {}, {}
        for role, name in FILES.items():
            rows[role], hashes[role] = read_jsonl(archive, name)

    gold, silver = rows["gold"], rows["silver"]
    gold_full = Counter(row["sequences"]["full"] for row in gold)
    gold_motif = Counter(row["linker"]["sequence"] for row in gold)
    gold_length = Counter(row["linker"]["length"] for row in gold)
    gold_issues = Counter()
    for row in gold:
        full, link = row["sequences"]["full"], row["linker"]
        if full[link["start"] - 1:link["end"]] != link["sequence"]:
            gold_issues["linker_sequence_coordinate_mismatch"] += 1
        if row["sequence_length"] != len(full):
            gold_issues["declared_full_length_mismatch"] += 1
        domains = row["domains"]
        if domains["domain_A"]["end"] != link["start"] - 1:
            gold_issues["domain_A_boundary_mismatch"] += 1
        if domains["domain_B"]["start"] != link["end"] + 1:
            gold_issues["domain_B_boundary_mismatch"] += 1
        if domains["domain_B"]["end"] > len(full):
            gold_issues["domain_B_outside_full_sequence"] += 1

    silver_full = Counter(row["sequences"]["full"] for row in silver)
    silver_acc = Counter(row["protein_accession"] for row in silver)
    silver_pair = Counter((row["domain_A_acc"], row["domain_B_acc"]) for row in silver)
    silver_motif = Counter(row["sequences"]["linker"] for row in silver)
    silver_issues = Counter()
    silver_lengths = Counter()
    for row in silver:
        sequence = row["sequences"]
        full = sequence["full"]
        a0, a1 = row["domain_A_start"], row["domain_A_end"]
        b0, b1 = row["domain_B_start"], row["domain_B_end"]
        l0, l1 = row["linker_start"], row["linker_end"]
        link = sequence["linker"]
        silver_lengths[row["linker_length"]] += 1
        if len(full) != row["protein_length"]:
            silver_issues["interpro_uniprot_protein_length_mismatch"] += 1
        if not (1 <= a0 <= a1 < l0 <= l1 < b0 <= b1 <= len(full)):
            silver_issues["invalid_coordinate_order_or_bounds"] += 1
        if l0 != a1 + 1 or l1 != b0 - 1 or len(link) != row["linker_length"]:
            silver_issues["gap_geometry_mismatch"] += 1
        if sequence["domain_A"] != full[a0 - 1:a1] or link != full[l0 - 1:l1] or sequence["domain_B"] != full[b0 - 1:b1]:
            silver_issues["sequence_slice_mismatch"] += 1
        if any(aa not in "ACDEFGHIKLMNPQRSTVWY" for aa in link):
            silver_issues["noncanonical_linker_residue"] += 1
        if len(link) and sum(aa in "AILMFWVY" for aa in link) / len(link) >= 0.5:
            silver_issues["at_least_half_hydrophobic_linker"] += 1
        if row["domain_A_acc"] == "PF00001" or row["domain_B_acc"] == "PF00001":
            silver_issues["PF00001_seven_tm_receptor_family"] += 1

    overlap = None
    if split_manifest is not None:
        manifest = json.loads(split_manifest.read_text(encoding="utf-8"))
        base = split_manifest.parent
        known_accessions, known_sequences = set(), set()
        for role in ("train", "development", "auxiliary"):
            with (base / manifest["files"][role]["path"]).open(encoding="utf-8") as stream:
                for line in stream:
                    row = json.loads(line)
                    known_accessions.add(row["protein_id"])
                    known_sequences.add(row["sequence"])
        overlap = {
            "split_manifest_sha256": hashlib.sha256(split_manifest.read_bytes()).hexdigest(),
            "exact_accession_overlap": sum(row["protein_accession"] in known_accessions for row in silver),
            "exact_full_sequence_overlap": sum(row["sequences"]["full"] in known_sequences for row in silver),
            "limitation": "No domain/family or local-similarity screen; not a training firewall",
        }
    return {
        "archive_sha256": archive_sha, "archive_bytes": archive_path.stat().st_size,
        "zip_members": len(names), "member_sha256": hashes,
        "counts": {role: len(value) for role, value in rows.items()},
        "gold": {
            "unique_full_sequences": len(gold_full),
            "duplicate_full_sequence_rows": sum(n - 1 for n in gold_full.values()),
            "unique_linker_motifs": len(gold_motif),
            "linker_motif_counts": dict(gold_motif),
            "length_counts": dict(sorted(gold_length.items())),
            "issues": dict(gold_issues),
            "label_scope": "engineered_fusion_literature_candidates_not_curated_natural_DFL",
        },
        "silver": {
            "unique_accessions": len(silver_acc),
            "duplicate_accession_rows": sum(n - 1 for n in silver_acc.values()),
            "unique_full_sequences": len(silver_full),
            "duplicate_full_sequence_rows": sum(n - 1 for n in silver_full.values()),
            "unique_linker_sequences": len(silver_motif),
            "unique_ordered_pfam_pairs": len(silver_pair),
            "top_pfam_pairs": [{"a": a, "b": b, "rows": n} for (a, b), n in silver_pair.most_common(8)],
            "length_counts": dict(sorted(silver_lengths.items())),
            "issues": dict(silver_issues),
            "label_scope": "two_Pfam_match_gap_3_to_40_candidate_not_observed_DFL",
        },
        "overlap_with_b1_v2": overlap,
        "training_eligible_DFL_positive_rows": 0,
        "training_eligible_explanation": "No row supplies independently verified natural disordered-flexible-linker function, 1/0/x mask, and full benchmark/family firewall; this does not assert every gap is biologically non-linker.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--split-manifest", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.archive, args.split_manifest), indent=2, sort_keys=True))
