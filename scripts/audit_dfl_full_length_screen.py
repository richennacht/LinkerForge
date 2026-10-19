"""Compare old local firewall with a bounded full-length candidate rule.

This is a sequence-only sensitivity audit, not a validated family map or a
replacement for the binding B1 preflight. Benchmark labels are never parsed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

if __package__:
    from scripts.audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
    from scripts.preflight_dfl_family_training import local_hit
else:
    from audit_dfl_b1_source_capacity import caid3_sequences, two_line_fasta
    from preflight_dfl_family_training import local_hit


PROPOSED_RULE = {"matrix": "BLOSUM62", "gap_open": 10, "gap_extend": 1,
                 "minimum_identity": 0.30, "minimum_aligned_residues": 50,
                 "minimum_each_sequence_coverage": 0.80}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def full_length_hit(left: str, right: str) -> bool:
    if left == right or left in right or right in left:
        return True
    import parasail

    result = parasail.sw_trace_striped_16(left, right, 10, 1, parasail.blosum62)
    query, reference = result.traceback.query, result.traceback.ref
    aligned = len(query)
    if aligned < 50:
        return False
    identity = sum(a == b for a, b in zip(query, reference)) / aligned
    query_coverage = sum(a != "-" for a in query) / len(left)
    reference_coverage = sum(b != "-" for b in reference) / len(right)
    return identity >= 0.30 and min(query_coverage, reference_coverage) >= 0.80


def audit(root: Path) -> dict:
    root = root.resolve()
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
            raise ValueError(f"Source hash changed: {name}")
        sources.append((name, path, pool_manifest["files"][name]["records"]))
    if sha(catalog_path) != catalog_manifest["output_sha256"]:
        raise ValueError("Candidate catalog hash changed")
    sources.append(("disprot_candidate_catalog", catalog_path, catalog_manifest["records"]))
    caid_path = root / "artifacts/caid3-linker.fasta"
    if sha(caid_path) != pool_manifest["source_sha256"]["caid3"]:
        raise ValueError("CAID3 reference hash changed")
    dp_root = root / "data/processed/prospective_dfl_holdout_dp94_v1"
    dp_manifest_path = dp_root / "manifest.json"
    dp_manifest = json.loads(dp_manifest_path.read_text(encoding="utf-8"))
    if dp_manifest["status"] != "candidate_not_certified_independent":
        raise ValueError("DP94 registry status changed")
    dp_path = dp_root / "sequences.fasta"
    dp94 = two_line_fasta(dp_path)
    canonical = dp_path.read_text(encoding="utf-8").replace("\r\n", "\n")
    if (hashlib.sha256(canonical.encode()).hexdigest() != dp_manifest["sequence_fasta_sha256"]
            or [key for key, _ in dp94] != dp_manifest["holdout_ids"]):
        raise ValueError("DP94 registry changed")
    registries = {"CAID3": caid3_sequences(caid_path), "DP94": dp94}
    candidates = {}
    for name, path, expected_count in sources:
        rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
        if len(rows) != expected_count:
            raise ValueError(f"Source count changed: {name}")
        for row in rows:
            identifier, sequence, labels = row["protein_id"], row["sequence"], row["labels"]
            if len(sequence) != len(labels) or set(labels) - set("01x"):
                raise ValueError(f"Invalid candidate mask: {name}/{identifier}")
            if "1" in labels:
                candidates.setdefault(sequence, set()).add(identifier)
    counts = Counter()
    per_registry = {name: Counter() for name in registries}
    for sequence, ids in candidates.items():
        counts["unique_positive_sequences"] += 1
        old_hit = False
        new_hit = False
        for name, registry in registries.items():
            registry_old = False
            registry_new = False
            for reference_id, reference_sequence in registry:
                exact = reference_id in ids or sequence == reference_sequence or sequence in reference_sequence or reference_sequence in sequence
                registry_old |= exact or local_hit(sequence, reference_sequence)
                registry_new |= exact or full_length_hit(sequence, reference_sequence)
                if registry_old and registry_new:
                    break
            old_hit |= registry_old
            new_hit |= registry_new
            per_registry[name]["old_local_flagged"] += registry_old
            per_registry[name]["proposed_full_length_flagged"] += registry_new
        counts["old_local_flagged_either"] += old_hit
        counts["proposed_full_length_flagged_either"] += new_hit
        counts["old_flagged_new_not_flagged"] += old_hit and not new_hit
        counts["new_flagged_old_not_flagged"] += new_hit and not old_hit
    return {"status": "sequence_only_sensitivity_not_B1_or_training_authorization",
            "candidate_unique_sequences": len(candidates), "counts": dict(counts),
            "by_registry": {name: dict(values) for name, values in per_registry.items()},
            "proposed_rule": PROPOSED_RULE,
            "source_hashes": {name: sha(path) for name, path, _ in sources},
            "caid3_reference_sha256": sha(caid_path),
            "dp94_sequence_manifest_sha256": sha(dp_manifest_path),
            "limitations": ["Proposed full-length rule is unvalidated and does not replace the binding local exclusion.",
                            "Domain/Pfam families and ambiguous low-complexity matches remain unaudited.",
                            "No benchmark labels, scores or per-protein errors were read."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.source_root), indent=2, sort_keys=True))
