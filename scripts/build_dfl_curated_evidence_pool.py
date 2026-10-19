"""Build a screened, provenance-bearing DisProt DFL evidence *candidate* pool.

The strict tier includes only validated, non-flagged DisProt flexible-linker
and experimental-order regions. A separate non-trainable review tier includes
validated, cited and released regions carrying an unresolved `unpublished`
flag. All other residues remain x. Benchmark labels are never opened.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from jsonschema import Draft202012Validator

if __package__:
    from scripts.audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
    from scripts.audit_dfl_disprot_family_metadata import DISPROT_SHA256, project_entry
    from scripts.preflight_dfl_family_training import RULE, local_hit
else:
    from audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
    from audit_dfl_disprot_family_metadata import DISPROT_SHA256, project_entry
    from preflight_dfl_family_training import RULE, local_hit


CONFIG = {"schema_version": "1.0.0", "source_release": "DisProt_2026_06",
          "positive_term": "IDPO:0000033", "negative_term": "IDPO:0000006",
          "unknown_policy": "all_unannotated_and_other_function_residues_are_x",
          "region_policy": "validated_cited_released_nonflagged_only_conflicting_residues_x",
          "review_policy": "validated_cited_released_unpublished_flag_separate_not_trainable",
          "benchmark_rule": RULE,
          "family_edge_types": ["pfam", "gene3D", "uniref50"]}
CAID3_SHA256 = "4c204ef49ffa12fe166883e6b2cff0012cabee21260b8f32f997081d0e2fa0f8"
DP94_MANIFEST_SHA256 = "fe64e7c688cfac48ed0dac73bcf09daae5cb0adfd0e03e9d530824d3d92c0051"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(value: dict) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def mask_entry(entry: dict, *, include_unpublished_flag: bool = False) -> dict | None:
    sequence = entry["sequence"].upper()
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
        return None
    positive, ordered = set(), set()
    evidence = {"positive": [], "ordered": []}
    for region in entry.get("regions", []):
        term = region.get("term_id")
        if term not in (CONFIG["positive_term"], CONFIG["negative_term"]):
            continue
        if not region.get("validated"):
            continue
        if region.get("unpublished") and not include_unpublished_flag:
            continue
        if (not region.get("reference_id")
                or region.get("reference_source") not in {"pmid", "doi"}
                or not region.get("released")):
            continue
        start, end = region.get("start"), region.get("end")
        if not isinstance(start, int) or not isinstance(end, int) or not 1 <= start <= end <= len(sequence):
            raise ValueError(f"Invalid curated coordinates: {entry['disprot_id']}/{region.get('region_id')}")
        target = positive if term == CONFIG["positive_term"] else ordered
        target.update(range(start - 1, end))
        key = "positive" if term == CONFIG["positive_term"] else "ordered"
        evidence[key].append({"region_id": region["region_id"], "start": start, "end": end,
                              "term_id": term, "reference_id": region.get("reference_id"),
                              "ec_id": region.get("ec_id"), "released": region.get("released"),
                              "unpublished_flag": bool(region.get("unpublished"))})
    if not positive and not ordered:
        return None
    conflict = positive & ordered
    labels = ["x"] * len(sequence)
    for index in positive - conflict:
        labels[index] = "1"
    for index in ordered - conflict:
        labels[index] = "0"
    if "1" not in labels and "0" not in labels:
        return None
    return {"protein_id": entry["disprot_id"], "uniprot_accession": entry.get("acc"),
            "sequence": sequence, "labels": "".join(labels),
            "task": "DFL_curated_positive_and_experimental_order",
            "source": "DisProt_2026_06", "evidence_quality": "curated_nonflagged_regions_not_adjudicated_no_linker",
            "evidence": evidence, "conflicting_residues_masked_x": len(conflict)}


def family_components(metadata: dict[str, dict]) -> dict[str, str]:
    parent = {key: key for key in metadata}

    def find(key: str) -> str:
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(left: str, right: str) -> None:
        parent[find(left)] = find(right)

    for kind in CONFIG["family_edge_types"]:
        groups = defaultdict(list)
        for key, entry in metadata.items():
            for family_id in entry["families"][kind]:
                groups[family_id].append(key)
        for members in groups.values():
            for member in members[1:]:
                union(members[0], member)
    components = defaultdict(list)
    for key in metadata:
        components[find(key)].append(key)
    return {key: hashlib.sha256("|".join(sorted(components[find(key)])).encode()).hexdigest()[:20]
            for key in metadata}


def build(source_root: Path, output_dir: Path) -> dict:
    source_root, output_dir = source_root.resolve(), output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite dataset: {output_dir}")
    schema_path = Path(__file__).resolve().parents[1] / "schemas/dfl_curated_evidence_candidate.v1.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)
    schema_hash = sha(schema_path)
    raw_path = source_root / "data/raw/disprot/disprot_2026_06_regions.json"
    if sha(raw_path) != DISPROT_SHA256:
        raise ValueError("Pinned DisProt source changed")
    caid_path = source_root / "artifacts/caid3-linker.fasta"
    if sha(caid_path) != CAID3_SHA256:
        raise ValueError("Pinned CAID3 sequence reference changed")
    caid = caid3_sequences(caid_path)  # Discards the annotation line without parsing it.
    dp_root = source_root / "data/processed/prospective_dfl_holdout_dp94_v1"
    dp_manifest_path = dp_root / "manifest.json"
    if sha(dp_manifest_path) != DP94_MANIFEST_SHA256:
        raise ValueError("Pinned DP94 sequence registry manifest changed")
    dp_manifest = json.loads(dp_manifest_path.read_text(encoding="utf-8"))
    if dp_manifest.get("status") != "candidate_not_certified_independent":
        raise ValueError("DP94 sequence registry status changed")
    dp_path = dp_root / "sequences.fasta"
    if hashlib.sha256(dp_path.read_text(encoding="utf-8").replace("\r\n", "\n").encode()).hexdigest() != dp_manifest["sequence_fasta_sha256"]:
        raise ValueError("DP94 sequence registry hash changed")
    dp = two_line_fasta(dp_path)
    if [key for key, _ in dp] != dp_manifest["holdout_ids"]:
        raise ValueError("DP94 ID roster changed")
    references = caid + dp
    raw = json.loads(raw_path.read_text(encoding="utf-8"))
    metadata, by_id = {}, {}
    annotation_census = Counter()
    for entry in raw["data"]:
        key = entry["disprot_id"]
        if key in by_id:
            raise ValueError(f"Duplicate DisProt ID: {key}")
        by_id[key] = entry
        for region in entry.get("regions", []):
            term = region.get("term_id")
            if term in (CONFIG["positive_term"], CONFIG["negative_term"]):
                name = "linker" if term == CONFIG["positive_term"] else "order"
                annotation_census[f"{name}_regions_total"] += 1
                if region.get("validated"):
                    annotation_census[f"{name}_regions_validated"] += 1
                if (region.get("validated") and not region.get("unpublished")
                        and region.get("reference_id")
                        and region.get("reference_source") in {"pmid", "doi"}
                        and region.get("released")):
                    annotation_census[f"{name}_regions_strict_included"] += 1
                if region.get("validated") and region.get("unpublished"):
                    annotation_census[f"{name}_regions_unpublished_flag_with_reference"] += bool(region.get("reference_id"))
        try:
            metadata[key] = project_entry(entry)
        except ValueError:
            continue
    components = family_components(metadata)
    benchmark_components = {components[key] for key, sequence in references
                            if key in metadata and metadata[key]["sequence"] == sequence}
    source_hashes = {"disprot": sha(raw_path), "caid3_reference": sha(caid_path),
                     "dp94_sequence_manifest": sha(dp_manifest_path), "dp94_fasta": sha(dp_path)}
    config_hash = canonical_hash(CONFIG)
    pipeline_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    counts = Counter()
    accepted, quarantined, review_accepted, review_quarantined = [], [], [], []
    reference_ids = {ref_key for ref_key, _ in references}
    for key, entry in sorted(by_id.items()):
        strict = mask_entry(entry)
        expanded = mask_entry(entry, include_unpublished_flag=True)
        for tier, row in (("strict", strict), ("review", expanded)):
            if row is None or (tier == "review" and strict is not None
                               and row["labels"] == strict["labels"]):
                continue
            prefix = "" if tier == "strict" else "review_"
            if tier == "review":
                row["evidence_quality"] = "curator_validated_cited_but_unpublished_flag_requires_review"
            counts[f"{prefix}candidate_proteins"] += 1
            for label in "10x":
                counts[f"{prefix}candidate_{label}_residues"] += row["labels"].count(label)
            row["provenance"] = {"source_id": key, "input_sha256": source_hashes["disprot"],
                                 "pipeline_commit": pipeline_commit, "config_hash": config_hash,
                                 "schema_sha256": schema_hash,
                                 "python_version": platform.python_version()}
            row["schema_version"] = "1.0.0"
            row["record_type"] = "disprot_dfl_evidence_candidate"
            if key not in metadata or metadata[key]["sequence"] != row["sequence"]:
                reason = "unresolved_family_metadata"
            elif key in reference_ids or any(
                    local_hit(row["sequence"], ref_seq) for _, ref_seq in references):
                reason = "benchmark_id_or_direct_sequence_match"
            elif components[key] in benchmark_components:
                reason = "benchmark_connected_family_component"
            else:
                reason = None
            if reason:
                row["quarantine_reason"] = reason
                row["training_eligibility"] = "excluded_benchmark_or_unresolved"
                (quarantined if tier == "strict" else review_quarantined).append(row)
                counts[f"{prefix}quarantine_{reason}"] += 1
            else:
                row["family_component_id"] = components[key]
                row["training_eligibility"] = ("candidate_not_B1" if tier == "strict"
                                               else "review_only_not_trainable")
                (accepted if tier == "strict" else review_accepted).append(row)
                counts[f"{prefix}screened_proteins"] += 1
                for label in "10x":
                    counts[f"{prefix}screened_{label}_residues"] += row["labels"].count(label)
            if len(row["sequence"]) != len(row["labels"]):
                raise ValueError(f"Sequence/label length mismatch: {key}")
            validator.validate(row)
    output_dir.mkdir(parents=True)
    files = {}
    for name, rows in (("screened.jsonl", accepted), ("quarantine.jsonl", quarantined),
                       ("review_screened.jsonl", review_accepted),
                       ("review_quarantine.jsonl", review_quarantined)):
        path = output_dir / name
        path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows), encoding="utf-8", newline="\n")
        files[name] = {"records": len(rows), "sha256": sha(path)}
    manifest = {"schema_version": "1.0.0", "status": "screened_curated_evidence_not_B1",
                "config": CONFIG, "config_hash": config_hash, "source_sha256": source_hashes,
                "record_schema_sha256": schema_hash,
                "source_rights": "DisProt CC BY 4.0; attribution required",
                "pipeline_commit": pipeline_commit, "python_version": platform.python_version(),
                "benchmark_sequence_only": {"caid3_count": len(caid), "dp94_count": len(dp)},
                "family_method": "transitive_shared_Pfam_Gene3D_UniRef50_in_DisProt_2026_06; incomplete_remote_homology_screen",
                "files": files, "counts": dict(counts), "annotation_census": dict(annotation_census),
                "limitations": ["Not an independently certified family split or B1 training authorization.",
                                "Order regions are explicit non-DFL residue evidence, not proof of a no-linker protein.",
                                "Review tier includes curator-validated cited regions with an unresolved unpublished flag and is not trainable.",
                                "No benchmark labels were read; CAID3 has a discarded annotation line.",
                                "Full-length local similarity may flag unrelated low-complexity segments."]}
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    return {"manifest": str(manifest_path), "manifest_sha256": sha(manifest_path),
            "status": manifest["status"], "counts": dict(counts), "files": files}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.source_root, args.output_dir), indent=2, sort_keys=True))
