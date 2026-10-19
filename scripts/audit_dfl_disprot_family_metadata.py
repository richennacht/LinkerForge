"""Audit Pfam/Gene3D/UniRef family evidence without using benchmark labels.

The pinned DisProt release JSON contains region annotations. This script
loads it but projects only entry ID, accession, sequence and existing family
metadata; it never accesses the regions or consensus fields. Missing or
changed-sequence entry metadata is unresolved, not guessed from a similar ID.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

if __package__:
    from scripts.preflight_dfl_family_training import local_hit
else:
    from preflight_dfl_family_training import local_hit


DISPROT_SHA256 = "1fdfce26a3ea9b626e63559f36abd86e023b50b612dc2c8756acbb08edead140"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def project_entry(entry: dict) -> dict:
    # Deliberately do not access entry['regions'] or entry['disprot_consensus'].
    sequence = entry["sequence"].upper()
    features = entry.get("features") or {}
    families = {}
    for source in ("pfam", "gene3D"):
        intervals = features.get(source) or []
        if not isinstance(intervals, list):
            raise ValueError(f"Invalid {source} feature list: {entry['disprot_id']}")
        values = set()
        for item in intervals:
            if not (isinstance(item, dict) and isinstance(item.get("id"), str)
                    and isinstance(item.get("start"), int) and isinstance(item.get("end"), int)
                    and 1 <= item["start"] <= item["end"] <= len(sequence)):
                raise ValueError(f"Invalid {source} feature coordinates: {entry['disprot_id']}")
            values.add(item["id"])
        families[source] = sorted(values)
    families["uniref50"] = [entry["uniref50"]] if entry.get("uniref50") else []
    return {"disprot_id": entry["disprot_id"], "accession": entry.get("acc"),
            "sequence": sequence, "families": families}


def audit(source_root: Path, input_dir: Path, cluster_tsv: Path, expected_cluster_hash: str,
          output_manifest: Path | None = None) -> dict:
    source_root, input_dir, cluster_tsv = source_root.resolve(), input_dir.resolve(), cluster_tsv.resolve()
    raw_path = source_root / "data/raw/disprot/disprot_2026_06_regions.json"
    if sha(raw_path) != DISPROT_SHA256:
        raise ValueError("Pinned DisProt source hash changed")
    input_manifest_path = input_dir / "manifest.json"
    input_manifest = json.loads(input_manifest_path.read_text(encoding="utf-8"))
    if input_manifest["status"] != "sequence_only_family_audit_not_split":
        raise ValueError("Unexpected MMseqs input status")
    records_path = input_dir / "records.jsonl"
    if sha(records_path) != input_manifest["files_sha256"]["records.jsonl"]:
        raise ValueError("MMseqs input records hash changed")
    if sha(cluster_tsv) != expected_cluster_hash:
        raise ValueError("Cluster TSV hash changed")
    records = [json.loads(line) for line in records_path.read_text(encoding="utf-8").splitlines()]
    by_record = {row["record_id"]: row for row in records}
    if len(by_record) != len(records):
        raise ValueError("Duplicate audit input record ID")
    needed_ids = {source["protein_id"] for row in records for source in row.get("source_ids", [])}
    needed_ids.update(row["source_id"] for row in records if "source_id" in row)
    raw = json.loads(raw_path.read_text(encoding="utf-8"))
    by_disprot = {}
    invalid_ids = set()
    for entry in raw["data"]:
        if entry["disprot_id"] not in needed_ids:
            continue
        try:
            projected = project_entry(entry)
        except ValueError:
            invalid_ids.add(entry["disprot_id"])
            continue
        if projected["disprot_id"] in by_disprot:
            raise ValueError("Duplicate DisProt ID")
        by_disprot[projected["disprot_id"]] = projected
    metadata = {}
    counts = Counter()
    counts["source_ids_with_invalid_family_coordinates"] = len(invalid_ids)
    for row in records:
        identifiers = ([source["protein_id"] for source in row["source_ids"]]
                       if "source_ids" in row else [row["source_id"]])
        exact = [by_disprot[key] for key in identifiers
                 if key in by_disprot and by_disprot[key]["sequence"] == row["sequence"]]
        counts["records"] += 1
        if not exact:
            counts["unresolved_id_or_sequence"] += 1
            continue
        families = {name: set() for name in ("pfam", "gene3D", "uniref50")}
        for entry in exact:
            for name, ids in entry["families"].items():
                families[name].update(ids)
        metadata[row["record_id"]] = families
        counts["verified_exact_sequence_metadata"] += 1
        for name in families:
            counts[f"records_with_{name}"] += bool(families[name])
    cluster_reps = {}
    for line in cluster_tsv.read_text(encoding="utf-8").splitlines():
        rep, member = line.split("\t")
        if rep not in by_record or member not in by_record or member in cluster_reps:
            raise ValueError("Malformed or duplicate MMseqs cluster entry")
        cluster_reps[member] = rep
    if set(cluster_reps) != set(by_record):
        raise ValueError("MMseqs clusters omit input records")
    protected = {key for key, row in by_record.items() if row["role"] == "sealed_benchmark_sequence_only"}
    candidate = set(by_record) - protected
    by_family = {name: defaultdict(set) for name in ("pfam", "gene3D", "uniref50")}
    for key, families in metadata.items():
        for name, ids in families.items():
            for family_id in ids:
                by_family[name][family_id].add(key)
    for name, index in by_family.items():
        touched = set()
        for members in index.values():
            if members & protected:
                touched.update(members & candidate)
        counts[f"candidate_sequences_sharing_{name}_with_benchmark"] = len(touched)
    mmseq_benchmark_reps = {cluster_reps[key] for key in protected}
    mmseq_flagged = {key for key in candidate if cluster_reps[key] in mmseq_benchmark_reps}
    counts["candidate_sequences_in_mmseq_benchmark_clusters"] = len(mmseq_flagged)
    mapped_candidate = candidate & set(metadata)
    counts["candidate_sequences_without_exact_metadata"] = len(candidate - mapped_candidate)
    union_family_flags = set(mmseq_flagged)
    for name, index in by_family.items():
        for members in index.values():
            if members & protected:
                union_family_flags.update(members & candidate)
    counts["candidate_sequences_directly_flagged_by_mmseq_or_any_family"] = len(union_family_flags)
    counts["candidate_sequences_mapped_and_not_directly_flagged"] = len(mapped_candidate - union_family_flags)
    parent = {key: key for key in by_record}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(left, right):
        parent[find(left)] = find(right)

    by_mmseq_rep = defaultdict(list)
    for member, rep in cluster_reps.items():
        by_mmseq_rep[rep].append(member)
    for groups in [by_mmseq_rep, *by_family.values()]:
        for members in groups.values():
            members = list(members)
            for other in members[1:]:
                union(members[0], other)
    benchmark_components = {find(key) for key in protected}
    mapped_survivors = {key for key in mapped_candidate if find(key) not in benchmark_components}
    survivor_components = Counter(find(key) for key in mapped_survivors)
    counts["candidate_sequences_mapped_outside_transitive_benchmark_components"] = len(mapped_survivors)
    counts["candidate_only_mapped_components"] = len(survivor_components)
    counts["largest_candidate_only_component"] = max(survivor_components.values(), default=0)
    protected_rows = [by_record[key] for key in sorted(protected)]
    binding_survivors = set()
    for key in sorted(mapped_survivors):
        row = by_record[key]
        source_ids = {source["protein_id"] for source in row["source_ids"]}
        if any(ref["source_id"] in source_ids or local_hit(row["sequence"], ref["sequence"])
               for ref in protected_rows):
            continue
        binding_survivors.add(key)
    binding_components = Counter(find(key) for key in binding_survivors)
    counts["candidate_sequences_outside_family_components_and_binding_direct_local"] = len(binding_survivors)
    counts["candidate_components_after_both_binding_and_family_screens"] = len(binding_components)
    counts["candidate_sequences_family_clear_but_binding_local_collision"] = len(mapped_survivors - binding_survivors)
    by_source = defaultdict(set)
    for key in binding_survivors:
        for source in by_record[key]["source_ids"]:
            by_source[source["source"]].add(key)
    counts["candidate_sequences_after_both_screens_by_source"] = {
        source: len(keys) for source, keys in sorted(by_source.items())
    }
    for registry in ("CAID3", "DP94"):
        registry_keys = {key for key in protected if by_record[key]["registry"] == registry}
        registry_components = {find(key) for key in registry_keys}
        registry_rows = [by_record[key] for key in sorted(registry_keys)]
        registry_survivors = set()
        for key in sorted(mapped_candidate):
            if find(key) in registry_components:
                continue
            row = by_record[key]
            source_ids = {source["protein_id"] for source in row["source_ids"]}
            if any(ref["source_id"] in source_ids or local_hit(row["sequence"], ref["sequence"])
                   for ref in registry_rows):
                continue
            registry_survivors.add(key)
        counts[f"scenario_{registry.lower()}_only_mapped_binding_and_family_survivors"] = len(registry_survivors)
    result = {"status": "metadata_family_sensitivity_not_B1",
            "disprot_source_sha256": sha(raw_path), "input_manifest_sha256": sha(input_manifest_path),
            "cluster_tsv_sha256": sha(cluster_tsv), "counts": dict(counts),
            "limitations": ["Shared Pfam/Gene3D/UniRef50 IDs are a conservative family screen, not a complete remote-homology certificate.",
                            "Unresolved IDs or changed sequences are not eligible for automatic split assignment.",
                            "The source JSON contains annotations, but this script never accesses its region/consensus fields.",
                            "The binding local preflight remains unchanged; this is not training authorization."]}
    if output_manifest is not None:
        output_manifest = output_manifest.resolve()
        if output_manifest.exists():
            raise FileExistsError(f"Refusing to overwrite capacity manifest: {output_manifest}")
        if not output_manifest.parent.is_dir():
            raise ValueError("Capacity manifest parent directory does not exist")
        component_members = defaultdict(list)
        for key in binding_survivors:
            component_members[find(key)].append(key)
        candidates = []
        for key in sorted(binding_survivors):
            row = by_record[key]
            members = sorted(component_members[find(key)])
            component_id = hashlib.sha256("|".join(members).encode()).hexdigest()[:20]
            candidates.append({"record_id": key, "component_id": component_id,
                               "sequence_sha256": row["sequence_sha256"],
                               "source_ids": row["source_ids"]})
        artifact = {"schema_version": "1.0.0", "status": "screened_capacity_not_training_split",
                    "screening": result, "candidates": candidates,
                    "warning": "No label masks, training allocation or negative semantics are certified here."}
        output_manifest.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        result["capacity_manifest_sha256"] = sha(output_manifest)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--cluster-tsv", type=Path, required=True)
    parser.add_argument("--cluster-tsv-sha256", required=True)
    parser.add_argument("--output-manifest", type=Path,
                        help="Write an immutable, sequence-hash-only capacity manifest; never a training split")
    args = parser.parse_args()
    print(json.dumps(audit(args.source_root, args.input_dir, args.cluster_tsv,
                           args.cluster_tsv_sha256, args.output_manifest), indent=2, sort_keys=True))
