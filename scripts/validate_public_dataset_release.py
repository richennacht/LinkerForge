"""Validate the two historical public dataset packages without training a model."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import zipfile


ROLES = ("train", "development", "auxiliary")
MINED_PARTS = (
    "mined_candidates_jsonl",
    "natural_authorities_jsonl",
    "template_authorities_jsonl",
)


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def rows(payload: bytes) -> list[dict]:
    return [json.loads(line) for line in io.BytesIO(payload) if line.strip()]


def validate_b1(directory: Path) -> dict:
    manifest_bytes = (directory / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest["status"] != "preregistered_v2_family_disjoint_split":
        raise ValueError("Unexpected B1 split status")
    counts = {}
    masks = {}
    for role in ROLES:
        entry = manifest["files"][role]
        payload = (directory / entry["path"]).read_bytes()
        if digest(payload) != entry["sha256"]:
            raise ValueError(f"B1 {role} hash mismatch")
        items = rows(payload)
        if len(items) != entry["records"]:
            raise ValueError(f"B1 {role} row count mismatch")
        observed = Counter()
        for item in items:
            if len(item["sequence"]) != len(item["labels"]):
                raise ValueError(f"B1 {role} sequence/mask mismatch")
            if set(item["labels"]) - set("01x"):
                raise ValueError(f"B1 {role} invalid mask")
            observed.update(item["labels"])
        if dict(observed) != {k: v for k, v in manifest["mask_counts"][role].items() if v}:
            raise ValueError(f"B1 {role} mask counts mismatch")
        counts[role] = len(items)
        masks[role] = dict(observed)
    for name in ("family_manifest",):
        entry = manifest[name]
        if digest((directory / entry["path"]).read_bytes()) != entry["sha256"]:
            raise ValueError(f"B1 {name} hash mismatch")
    for entry in manifest["benchmark_registries"]:
        if digest((directory / entry["path"]).read_bytes()) != entry["sha256"]:
            raise ValueError(f"B1 {entry['name']} registry hash mismatch")
        registry = json.loads((directory / entry["path"]).read_bytes())
        fasta = (directory / "registries" / registry["fasta_path"]).read_bytes()
        if digest(fasta.replace(b"\r\n", b"\n")) != registry["fasta_sha256_canonical_lf"]:
            raise ValueError(f"B1 {entry['name']} FASTA hash mismatch")
    return {"manifest_sha256": digest(manifest_bytes), "records": counts, "mask_counts": masks}


def validate_mined(archive: Path, manifest_path: Path) -> dict:
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest["dataset_id"] != "mined_and_natural_linkers_v1":
        raise ValueError("Unexpected mined dataset ID")
    counts = {}
    tiers = Counter()
    with zipfile.ZipFile(archive) as bundle:
        if set(bundle.namelist()) != {"manifest.json"} | {
            entry["path"] for entry in manifest["files"].values()
        }:
            raise ValueError("Mined archive member set mismatch")
        if bundle.read("manifest.json") != manifest_bytes:
            raise ValueError("Mined manifest differs from archived manifest")
        for key, entry in manifest["files"].items():
            payload = bundle.read(entry["path"])
            if digest(payload) != entry["sha256"]:
                raise ValueError(f"Mined {key} hash mismatch")
            items = rows(payload)
            if len(items) != entry["records"]:
                raise ValueError(f"Mined {key} row count mismatch")
            counts[key] = len(items)
            if key == "all_records_jsonl":
                tiers.update(item["authority_tier"] for item in items)
                if len({item["record_id"] for item in items}) != len(items):
                    raise ValueError("Mined record IDs are not unique")
    if counts["all_records_jsonl"] != sum(counts[k] for k in MINED_PARTS):
        raise ValueError("Mined partition counts do not sum to all records")
    for tier, count in tiers.items():
        if count != manifest["counts_by_tier"][tier]:
            raise ValueError(f"Mined {tier} count mismatch")
    return {"manifest_sha256": digest(manifest_bytes), "records": counts, "tiers": dict(tiers)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", type=Path, default=Path(__file__).resolve().parents[1] / "datasets")
    args = parser.parse_args()
    print(json.dumps({
        "classifier_b1_v2": validate_b1(args.datasets / "dfl_b1_v2_family_split"),
        "mined_and_natural_v1": validate_mined(
            args.datasets / "mined_and_natural_linkers_v1.zip",
            args.datasets / "mined_and_natural_linkers_v1.manifest.json",
        ),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
