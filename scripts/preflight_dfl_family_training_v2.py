"""Sequence-only, fail-closed V2 gate for a preregistered family-disjoint DFL split.

Enforces `literature/DFL_FIREWALL_V2_PREREGISTRATION_2026-10-18.md`:
1. True ungapped reciprocal residue coverage (`RULE_V2`) fixing V1's unbounded
   gapped-length ratio (`hit.length / min(len)` > 1.0).
2. Mandatory shuffle-null specificity verification (`<= 0.005` false-hit rate).
3. Transitive domain-family (`Pfam`, `Gene3D`, `UniRef50`, optional `MMseqs2`)
   connected-component disjunction against sealed `caid3` and `dp94` registries
   and across `train`, `development`, and `auxiliary` (`order_calibration_panel`).
4. Explicit `1` (`IDPO:0000033`) and `0` (`IDPO:0000006` experimental-order)
   supervision in both `train` and `development`, with `x` never relabelled `0`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from jsonschema import Draft202012Validator


ROLES = ("train", "development", "auxiliary")
REQUIRED_REGISTRIES = frozenset({"caid3", "dp94"})
REGISTRY_STATES = {
    "caid3": ("sequence_only_protected", "retrospective_spent_exclusion_only"),
    "dp94": ("sequence_only_labels_sealed", "unopened"),
}
RULE_V2 = {
    "rule_version": "2.0.0",
    "matrix": "BLOSUM62",
    "gap_open": 10,
    "gap_extend": 1,
    "minimum_identity": 0.30,
    "minimum_aligned_residues": 50,
    "minimum_each_sequence_residue_coverage": 0.80,
    "coverage_definition": "ungapped_aligned_residues_divided_by_sequence_length",
    "maximum_shuffle_null_rate": 0.005,
}


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def anchored(base: Path, name: str) -> Path:
    path = (base / name).resolve()
    if not path.is_relative_to(base.resolve()):
        raise ValueError(f"Path escapes manifest directory: {name}")
    return path


def pinned_json(base: Path, info: dict) -> tuple[Path, dict]:
    path = anchored(base, info["path"])
    raw = path.read_bytes()
    if digest(raw) != info["sha256"]:
        raise ValueError(f"Hash mismatch: {path}")
    return path, json.loads(raw)


def fasta(path: Path, expected_hash: str) -> list[tuple[str, str]]:
    raw = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    if digest(raw.encode("utf-8")) != expected_hash:
        raise ValueError(f"FASTA hash mismatch: {path}")
    lines = raw.splitlines()
    if not lines or len(lines) % 2:
        raise ValueError(f"Expected two-line FASTA: {path}")
    rows = []
    for at in range(0, len(lines), 2):
        if not lines[at].startswith(">") or not lines[at + 1]:
            raise ValueError(f"Malformed FASTA: {path}")
        identifier = lines[at][1:].split()[0]
        sequence = lines[at + 1].upper()
        if not identifier or set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
            raise ValueError(f"Invalid FASTA record: {path}")
        rows.append((identifier, sequence))
    if len({key for key, _ in rows}) != len(rows):
        raise ValueError(f"Duplicate FASTA identifier: {path}")
    return rows


def v2_alignment_metrics(left: str, right: str) -> dict[str, float | bool]:
    """Compute true ungapped residue coverage and identity under RULE_V2."""
    left, right = left.upper(), right.upper()
    if not left or not right:
        raise ValueError("Sequences must be non-empty")
    if left == right or left in right or right in left:
        return {
            "exact_or_contained": True,
            "aligned_columns": float(min(len(left), len(right))),
            "identity": 1.0,
            "left_residue_coverage": 1.0 if left in right else len(right) / len(left),
            "right_residue_coverage": 1.0 if right in left else len(left) / len(right),
            "hit": True,
        }
    try:
        import parasail
    except ImportError as error:
        raise ValueError("parasail is required for the V2 sequence firewall") from error
    result = parasail.sw_trace_striped_16(left, right, 10, 1, parasail.blosum62)
    query = result.traceback.query
    reference = result.traceback.ref
    aligned_columns = len(query)
    if aligned_columns < RULE_V2["minimum_aligned_residues"]:
        return {
            "exact_or_contained": False,
            "aligned_columns": float(aligned_columns),
            "identity": 0.0,
            "left_residue_coverage": 0.0,
            "right_residue_coverage": 0.0,
            "hit": False,
        }
    matches = sum(a == b for a, b in zip(query, reference))
    left_ungapped = sum(a != "-" for a in query)
    right_ungapped = sum(b != "-" for b in reference)
    identity = matches / aligned_columns
    left_coverage = left_ungapped / len(left)
    right_coverage = right_ungapped / len(right)
    if not (0.0 <= left_coverage <= 1.0 and 0.0 <= right_coverage <= 1.0):
        raise ValueError("Ungapped residue coverage out of [0, 1] bounds")
    hit = (
        identity >= RULE_V2["minimum_identity"]
        and min(left_coverage, right_coverage) >= RULE_V2["minimum_each_sequence_residue_coverage"]
    )
    return {
        "exact_or_contained": False,
        "aligned_columns": float(aligned_columns),
        "identity": identity,
        "left_residue_coverage": left_coverage,
        "right_residue_coverage": right_coverage,
        "hit": hit,
    }


def v2_homolog_hit(left: str, right: str) -> bool:
    return bool(v2_alignment_metrics(left, right)["hit"])


def verify(manifest_path: Path) -> dict:
    manifest_path = manifest_path.resolve()
    base = manifest_path.parent
    raw_manifest = manifest_path.read_bytes()
    split = json.loads(raw_manifest)
    schema_path = Path(__file__).resolve().parents[1] / "schemas/dfl_family_split_v2.schema.json"
    if schema_path.is_file():
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
        Draft202012Validator.check_schema(schema)
        Draft202012Validator(schema).validate(split)
    if (split.get("schema_version") != "2.0.0"
            or split.get("status") != "preregistered_v2_family_disjoint_split"):
        raise ValueError("Unexpected V2 split schema/status")
    if split.get("similarity_rule") != RULE_V2:
        raise ValueError("V2 similarity rule mismatch")
    null_check = split.get("shuffle_null_verification") or {}
    max_null = RULE_V2["maximum_shuffle_null_rate"]
    if (null_check.get("total_pairs", 0) < 300
            or null_check.get("residue_shuffle_rate", 1.0) > max_null
            or null_check.get("block12_shuffle_rate", 1.0) > max_null):
        raise ValueError("Shuffle-null specificity gate failed or missing")

    policies = split.get("source_policies")
    if not isinstance(policies, dict) or not policies:
        raise ValueError("DFL source policies missing")
    for source, policy in policies.items():
        if (not source
                or policy.get("task") not in {
                    "DFL_curated_positive_and_experimental_order",
                    "DFL_experimental_order_negative_panel",
                    "DFL_positive_unlabeled",
                }
                or policy.get("zero_semantics") not in {"explicit_negative_curated", "none"}):
            raise ValueError(f"Incomplete V2 DFL source policy: {source}")
        if policy["task"] == "DFL_positive_unlabeled" and policy["zero_semantics"] != "none":
            raise ValueError(f"Positive-only source asserts negatives: {source}")

    entries = split.get("benchmark_registries")
    if not isinstance(entries, list):
        raise ValueError("Benchmark registries missing")
    names = [entry["name"] for entry in entries]
    if len(set(names)) != len(names) or not REQUIRED_REGISTRIES.issubset(names):
        raise ValueError("Required sealed benchmark registry omitted or duplicated")
    references = []
    registry_hashes = {}
    for entry in entries:
        path, registry = pinned_json(base, entry)
        expected_status, expected_seal = REGISTRY_STATES.get(
            entry["name"], ("sequence_only_labels_sealed", "unopened")
        )
        if (registry.get("schema_version") != "1.0.0"
                or registry.get("status") != expected_status
                or registry.get("label_seal_state") != expected_seal
                or not registry.get("source_release")):
            raise ValueError(f"Registry seal/status mismatch or incomplete registry: {path}")
        records = fasta(anchored(path.parent, registry["fasta_path"]), registry["fasta_sha256_canonical_lf"])
        if [key for key, _ in records] != registry["ids"]:
            raise ValueError(f"Registry ID roster changed: {path}")
        references.extend((entry["name"], key, sequence) for key, sequence in records)
        registry_hashes[entry["name"]] = entry["sha256"]
    if len({key for _, key, _ in references}) != len(references):
        raise ValueError("Ambiguous benchmark identifier across registries")

    _, family = pinned_json(base, split["family_manifest"])
    if (family.get("schema_version") != "2.0.0"
            or family.get("status") != "family_clustered_v2"
            or not family.get("method")
            or not family.get("version")
            or not family.get("source_sha256")
            or not isinstance(family.get("membership"), dict)):
        raise ValueError("V2 family manifest lacks versioned method or provenance")
    membership = family["membership"]
    if any(not isinstance(value, str) or not value for value in membership.values()):
        raise ValueError("Empty family assignment")
    if any(key not in membership for _, key, _ in references):
        raise ValueError("Benchmark lacks family assignment")
    benchmark_families = {membership[key] for _, key, _ in references}

    ids_by_role, sequences_by_role, families_by_role, counts, mask_counts = {}, {}, {}, {}, {}
    for role in ROLES:
        info = split["files"][role]
        raw = anchored(base, info["path"]).read_bytes()
        if digest(raw) != info["sha256"]:
            raise ValueError(f"Split file hash mismatch: {role}")
        rows = [json.loads(line) for line in raw.decode("utf-8").splitlines() if line]
        if not rows or len(rows) != info["records"]:
            raise ValueError(f"Empty or mismatched split: {role}")
        identifiers, sequences, families = set(), set(), set()
        masks = {"1": 0, "0": 0, "x": 0}
        for row in rows:
            key, sequence, labels = row["protein_id"], row["sequence"].upper(), row["labels"]
            if not key or key in identifiers or key not in membership:
                raise ValueError(f"Duplicate or unclustered identifier: {role}/{key}")
            if not sequence or len(sequence) != len(labels) or set(labels) - set("01x"):
                raise ValueError(f"Invalid 1/0/x mask: {role}/{key}")
            if "1" not in labels and "0" not in labels:
                raise ValueError(f"All-unknown record in split: {role}/{key}")
            for label in masks:
                masks[label] += labels.count(label)
            policy = policies.get(row.get("source"))
            if not policy or row.get("task") != policy["task"]:
                raise ValueError(f"Missing or mismatched source policy: {role}/{key}")
            if "0" in labels and policy["zero_semantics"] != "explicit_negative_curated":
                raise ValueError(f"Unverified negative mask: {role}/{key}")
            provenance = row.get("provenance")
            if not isinstance(provenance, dict) or not provenance.get("source_id") or not provenance.get("input_sha256"):
                raise ValueError(f"Missing record provenance: {role}/{key}")
            if row.get("family_component_id") != membership[key]:
                raise ValueError(f"Row family_component_id mismatch: {role}/{key}")
            if membership[key] in benchmark_families:
                raise ValueError(f"Benchmark-family collision: {role}/{key}")
            for name, ref_key, ref_sequence in references:
                if key == ref_key or v2_homolog_hit(sequence, ref_sequence):
                    raise ValueError(f"Benchmark V2 sequence collision: {role}/{key} vs {name}/{ref_key}")
            identifiers.add(key)
            sequences.add(sequence)
            families.add(membership[key])
        ids_by_role[role], sequences_by_role[role], families_by_role[role] = identifiers, sequences, families
        counts[role] = len(rows)
        mask_counts[role] = masks

    if not mask_counts["train"]["1"] or not mask_counts["development"]["1"]:
        raise ValueError("Train and development each require observed DFL positives ('1')")
    if not mask_counts["train"]["0"] or not mask_counts["development"]["0"]:
        raise ValueError("Train and development each require verified experimental-order negatives ('0')")
    if not mask_counts["auxiliary"]["0"] or mask_counts["auxiliary"]["1"] != 0:
        raise ValueError("Auxiliary order_calibration_panel must contain verified '0' residues and zero '1' residues")

    for at, left in enumerate(ROLES):
        for right in ROLES[at + 1:]:
            if ids_by_role[left] & ids_by_role[right] or sequences_by_role[left] & sequences_by_role[right]:
                raise ValueError(f"Cross-split exact overlap: {left}/{right}")
            if families_by_role[left] & families_by_role[right]:
                raise ValueError(f"Cross-split family overlap: {left}/{right}")
            for left_seq in sequences_by_role[left]:
                for right_seq in sequences_by_role[right]:
                    if v2_homolog_hit(left_seq, right_seq):
                        raise ValueError(f"Cross-split V2 sequence homolog overlap: {left}/{right}")

    return {
        "status": "pass_v2_preregistered_sequence_and_family_firewall",
        "split_manifest_sha256": digest(raw_manifest),
        "family_manifest_sha256": split["family_manifest"]["sha256"],
        "registry_sha256": registry_hashes,
        "registry_states": {
            name: REGISTRY_STATES.get(name, ("sequence_only_labels_sealed", "unopened"))[1]
            for name in names
        },
        "shuffle_null_verification": null_check,
        "counts": counts,
        "mask_counts": mask_counts,
        "gates": {
            "verified_negative_residues": mask_counts["train"]["0"] > 0 and mask_counts["development"]["0"] > 0,
            "frozen_family_disjoint_train_development_split": True,
            "independently_adjudicated_order_negative_panel": mask_counts["auxiliary"]["0"] > 0,
            "benchmark_firewall_preflight_for_split": True,
        },
        "limitations": split.get("limitations", []),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(verify(args.manifest), indent=2, sort_keys=True))
    except Exception as error:
        print(json.dumps({"status": "blocked", "reason": str(error)}, indent=2, sort_keys=True))
        raise SystemExit(2) from error
