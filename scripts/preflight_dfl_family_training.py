"""Sequence-only, fail-closed gate for a future family-disjoint DFL split.

This is a prerequisite, not a certificate of family quality or permission to
train. Benchmark labels are never opened. The launcher must invoke this before
allocating GPU and must also run the repo's detector status check.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROLES = ("train", "development", "auxiliary")
REQUIRED_REGISTRIES = frozenset({"caid3", "dp94"})
REGISTRY_STATES = {
    "caid3": ("sequence_only_protected", "retrospective_spent_exclusion_only"),
    "dp94": ("sequence_only_labels_sealed", "unopened"),
}
RULE = {"matrix": "BLOSUM62", "gap_open": 10, "gap_extend": 1,
        "minimum_identity": 0.25, "minimum_aligned_residues": 50,
        "minimum_shorter_sequence_coverage": 0.5}


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


def local_hit(left: str, right: str) -> bool:
    if left == right or left in right or right in left:
        return True
    try:
        import parasail
    except ImportError as error:
        raise ValueError("parasail is required for the local-similarity firewall") from error
    hit = parasail.sw_stats_striped_16(left, right, 10, 1, parasail.blosum62)
    return (hit.length >= 50 and hit.matches / hit.length >= 0.25
            and hit.length / min(len(left), len(right)) >= 0.5)


def verify(manifest_path: Path) -> dict:
    manifest_path = manifest_path.resolve()
    base = manifest_path.parent
    split = json.loads(manifest_path.read_bytes())
    if split.get("schema_version") != "1.0.0" or split.get("status") != "candidate_family_split_not_certified":
        raise ValueError("Unexpected split schema/status; old split is ineligible")
    if split.get("local_similarity_rule") != RULE:
        raise ValueError("Similarity rule changed")
    policies = split.get("source_policies")
    if not isinstance(policies, dict) or not policies:
        raise ValueError("DFL source policies missing")
    for source, policy in policies.items():
        if (not source or policy.get("task") not in {"DFL", "DFL_positive_unlabeled"}
                or policy.get("zero_semantics") not in {"explicit_negative_curated", "none", "published_negative_unverified"}):
            raise ValueError(f"Incomplete DFL source policy: {source}")
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
            entry["name"], ("sequence_only_labels_sealed", "unopened"))
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
    if (family.get("schema_version") != "1.0.0" or family.get("status") != "family_clustered"
            or not family.get("method") or not family.get("version")
            or not family.get("source_sha256") or not isinstance(family.get("membership"), dict)):
        raise ValueError("Family manifest lacks versioned method or provenance")
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
            if membership[key] in benchmark_families:
                raise ValueError(f"Benchmark-family collision: {role}/{key}")
            for name, ref_key, ref_sequence in references:
                if key == ref_key or local_hit(sequence, ref_sequence):
                    raise ValueError(f"Benchmark sequence collision: {role}/{key} vs {name}/{ref_key}")
            identifiers.add(key)
            sequences.add(sequence)
            families.add(membership[key])
        ids_by_role[role], sequences_by_role[role], families_by_role[role] = identifiers, sequences, families
        counts[role] = len(rows)
        mask_counts[role] = masks
    if not mask_counts["train"]["1"] or not mask_counts["development"]["1"]:
        raise ValueError("Train and development each require observed DFL positives")
    for at, left in enumerate(ROLES):
        for right in ROLES[at + 1:]:
            if ids_by_role[left] & ids_by_role[right] or sequences_by_role[left] & sequences_by_role[right]:
                raise ValueError(f"Cross-split exact overlap: {left}/{right}")
            if families_by_role[left] & families_by_role[right]:
                raise ValueError(f"Cross-split family overlap: {left}/{right}")
    return {"status": "pass_declared_sequence_and_family_manifest_only",
            "split_manifest_sha256": digest(manifest_path.read_bytes()),
            "family_manifest_sha256": split["family_manifest"]["sha256"],
            "registry_sha256": registry_hashes, "registry_states": {
                name: REGISTRY_STATES.get(name, ("sequence_only_labels_sealed", "unopened"))[1]
                for name in names}, "counts": counts, "mask_counts": mask_counts,
            "limitations": "Not a biological homology certificate, DFL benchmark approval, or training authorization by itself."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(verify(args.manifest), indent=2))
    except (KeyError, ValueError, OSError, TypeError, json.JSONDecodeError) as error:
        print(json.dumps({"status": "blocked", "reason": str(error)}, indent=2))
        raise SystemExit(2) from error
