"""Fail closed before CAID3/DP94-protected DFL training; never reads benchmark labels.

Supports:
1. The preregistered V2 family-disjoint split (`schema_version == "2.0.0"`,
   `status == "preregistered_v2_family_disjoint_split"`), delegating to
   `preflight_dfl_family_training_v2.verify` and additionally cross-checking
   that the split's `caid3` and `dp94` sequence registries match the root
   `artifacts/caid3-linker.fasta` and
   `data/processed/prospective_dfl_holdout_dp94_v1/sequences.fasta`.
2. Legacy V1 sequence-screened manifests (`caid3_protected_dfl_v1`), which fail
   closed due to prospective DP94 overlap (`DP00919` vs `DP03956`).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

if __package__:
    from scripts.build_caid3_protected_dfl_splits import OUTPUT, REFERENCE, ROOT, RULE, match, reference_sequences, sha
    from scripts.preflight_dfl_family_training_v2 import fasta as v2_fasta, pinned_json, verify as verify_v2
else:
    from build_caid3_protected_dfl_splits import OUTPUT, REFERENCE, ROOT, RULE, match, reference_sequences, sha
    from preflight_dfl_family_training_v2 import fasta as v2_fasta, pinned_json, verify as verify_v2


PROSPECTIVE = ROOT / "data/processed/prospective_dfl_holdout_dp94_v1"
V2_DEFAULT_OUTPUT = ROOT / "data/processed/dfl_b1_v2_family_split"


def prospective_references() -> list[tuple[str, str]]:
    """Read benchmark sequences only; never read its labels during training."""
    manifest_path = PROSPECTIVE / "manifest.json"
    fasta_path = PROSPECTIVE / "sequences.fasta"
    if not manifest_path.exists() or not fasta_path.exists():
        raise ValueError("Prospective benchmark sequence firewall is missing")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "candidate_not_certified_independent":
        raise ValueError("Unreviewed prospective benchmark manifest state")
    canonical = fasta_path.read_text(encoding="utf-8").replace("\r\n", "\n").encode("utf-8")
    if hashlib.sha256(canonical).hexdigest() != manifest["sequence_fasta_sha256"]:
        raise ValueError("Prospective benchmark sequence hash changed")
    lines = fasta_path.read_text(encoding="utf-8").splitlines()
    if len(lines) % 2:
        raise ValueError("Malformed prospective benchmark FASTA")
    references = []
    for at in range(0, len(lines), 2):
        if not lines[at].startswith(">") or not lines[at + 1]:
            raise ValueError("Malformed prospective benchmark FASTA record")
        references.append((lines[at][1:].split()[0], lines[at + 1].upper()))
    if len(references) != manifest["holdout_proteins"] or [r[0] for r in references] != manifest["holdout_ids"]:
        raise ValueError("Prospective benchmark IDs/count changed")
    return references


def _verify_v2_with_root_registries(manifest_path: Path) -> dict:
    manifest_path = manifest_path.resolve()
    if manifest_path == (OUTPUT / "manifest.json").resolve():
        raise ValueError("Deprecated split caid3_protected_dfl_v1 is prohibited")
    v2_result = verify_v2(manifest_path)
    refs = reference_sequences(REFERENCE)
    future_refs = prospective_references()
    base = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    by_name = {entry["name"]: entry for entry in manifest["benchmark_registries"]}
    caid_reg_path, caid_reg = pinned_json(base, by_name["caid3"])
    dp94_reg_path, dp94_reg = pinned_json(base, by_name["dp94"])
    caid_records = v2_fasta(caid_reg_path.parent / caid_reg["fasta_path"], caid_reg["fasta_sha256_canonical_lf"])
    dp94_records = v2_fasta(dp94_reg_path.parent / dp94_reg["fasta_path"], dp94_reg["fasta_sha256_canonical_lf"])
    if caid_records != refs:
        raise ValueError("V2 split caid3 registry does not match root artifacts/caid3-linker.fasta")
    if dp94_records != future_refs:
        raise ValueError("V2 split dp94 registry does not match root data/processed/prospective_dfl_holdout_dp94_v1")
    return {
        **v2_result,
        "reference_sha256": sha(REFERENCE),
        "prospective_reference_sha256": sha(PROSPECTIVE / "sequences.fasta"),
    }


def verify(manifest_path: Path) -> dict:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("schema_version") == "2.0.0"
            and manifest.get("status") == "preregistered_v2_family_disjoint_split"):
        return _verify_v2_with_root_registries(manifest_path)
    if manifest["status"] != "sequence_screened_not_full_homology_certification":
        raise ValueError("Unexpected data status")
    if manifest["rule"] != RULE or manifest["reference_sha256"] != sha(REFERENCE):
        raise ValueError("Blacklist rule or CAID3 reference changed")
    refs = reference_sequences(REFERENCE)
    future_refs = prospective_references()
    counts = {}
    ids = {}
    sequences = {}
    for name in ("train", "development", "auxiliary"):
        info = manifest["files"][name]
        path = ROOT / info["path"]
        if sha(path) != info["sha256"]:
            raise ValueError(f"Data hash changed: {name}")
        rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
        if len(rows) != info["retained_proteins"] or not rows:
            raise ValueError(f"Unexpected row count: {name}")
        for row in rows:
            if len(row["sequence"]) != len(row["labels"]):
                raise ValueError(f"Sequence/label mismatch: {name} {row['protein_id']}")
            if set(row["labels"]) - {"0", "1", "x"}:
                raise ValueError(f"Unknown label code: {name} {row['protein_id']}")
            hit = match(row, refs)
            if hit:
                raise ValueError(f"CAID3 blacklist violation in {name}: {hit}")
            future_hit = match(row, future_refs)
            if future_hit:
                future_hit["benchmark_id"] = future_hit.pop("caid3_id")
                raise ValueError(f"Prospective DP94 benchmark sequence overlap in {name}: {future_hit}")
        ids[name] = {row["protein_id"] for row in rows}
        sequences[name] = {row["sequence"].upper() for row in rows}
        counts[name] = len(rows)
    if ids["train"] & ids["development"] or sequences["train"] & sequences["development"]:
        raise ValueError("Exact train/development overlap")
    blacklist = manifest_path.parent / "blacklist.json"
    if sha(blacklist) != manifest["blacklist_sha256"]:
        raise ValueError("Blacklist file hash changed")
    return {"status": "pass_under_declared_sequence_screen", "reference_sha256": sha(REFERENCE),
            "prospective_reference_sha256": sha(PROSPECTIVE / "sequences.fasta"),
            "counts": counts, "limitations": manifest["limitations"]}


def failure_report(manifest_path: Path, error: Exception) -> dict:
    """Return a machine-readable fail-closed result without weakening the gate."""
    protected_inputs = [
        REFERENCE,
        PROSPECTIVE / "manifest.json",
        PROSPECTIVE / "sequences.fasta",
    ]
    return {
        "status": "failed_closed",
        "training_authorized": False,
        "manifest": str(manifest_path),
        "reason": str(error),
        "missing_protected_inputs": [str(path) for path in protected_inputs if not path.is_file()],
        "claim_boundary": (
            "Missing or mismatched protected sequence inputs prohibit training. "
            "The gate must never be bypassed or reconstructed from benchmark labels."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=V2_DEFAULT_OUTPUT / "manifest.json")
    args = parser.parse_args()
    try:
        result = verify(args.manifest)
    except (FileNotFoundError, KeyError, OSError, ValueError) as error:
        print(json.dumps(failure_report(args.manifest, error), indent=2, sort_keys=True))
        raise SystemExit(2) from None
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
