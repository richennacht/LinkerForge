#!/usr/bin/env python3
"""Build a deterministic, dual-curator review packet from the v1 linker corpus.

The packet reviews linker identity and source evidence.  It deliberately keeps
fusion compatibility as ``not_assessed`` and never promotes classifier-mined
records to experimental or curator-validated labels.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
from typing import Any


EVIDENCE_ROLES = {
    "mined_high_confidence_dfl_call": "classifier_prediction",
    "mined_provisional_dfl_candidate": "classifier_prediction",
    "natural_cath_pdb_verified": "natural_structural_observation",
    "natural_disprot_curated_dfl": "curated_dfl_annotation",
    "literature_engineered_template": "published_engineered_use",
}


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_rank(seed: str, record_id: str) -> str:
    return hashlib.sha256(f"{seed}\0{record_id}".encode("utf-8")).hexdigest()


def length_bin(length: int, boundaries: list[int]) -> str:
    lower = 1
    for upper in boundaries:
        if length <= upper:
            return f"{lower:02d}-{upper:02d}"
        lower = upper + 1
    return f"{lower:02d}-plus"


def score_band(row: dict[str, Any]) -> str:
    score = float(row.get("classifier_scores", {}).get("cnn_mean_score", 0.0))
    if score >= 0.7351981997489929:
        return "frozen_high_threshold_or_above"
    if score >= 0.6391201615333557:
        return "conformal_to_high_threshold"
    if score >= 0.5:
        return "provisional_to_conformal_threshold"
    return "below_provisional_or_not_scored"


def stratum(row: dict[str, Any], boundaries: list[int]) -> str:
    return "|".join(
        [
            str(row["authority_tier"]),
            str(row.get("mining_source_corpus") or "not_applicable"),
            length_bin(int(row["linker_length"]), boundaries),
            score_band(row),
        ]
    )


def stratified_select(
    records: list[dict[str, Any]],
    quota: int,
    *,
    seed: str,
    boundaries: list[int],
) -> list[dict[str, Any]]:
    if quota > len(records):
        raise ValueError(f"quota {quota} exceeds available records {len(records)}")
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in records:
        groups[stratum(row, boundaries)].append(row)
    for rows in groups.values():
        rows.sort(key=lambda row: stable_rank(seed, str(row["record_id"])))
    selected: list[dict[str, Any]] = []
    keys = sorted(groups)
    while len(selected) < quota:
        progressed = False
        for key in keys:
            if groups[key] and len(selected) < quota:
                selected.append(groups[key].pop(0))
                progressed = True
        if not progressed:
            raise RuntimeError("stratified selection exhausted before reaching quota")
    return selected


def empty_review() -> dict[str, Any]:
    return {
        "status": "pending",
        "curator_id": None,
        "reviewed_at": None,
        "source_checked": False,
        "coordinates_correct": "not_reviewed",
        "interdomain_connection": "not_reviewed",
        "linker_behavior_evidence": "not_reviewed",
        "structural_state": "not_reviewed",
        "source_support": "not_reviewed",
        "decision": "not_reviewed",
        "evidence_ids": [],
        "notes": "",
    }


def review_record(row: dict[str, Any], curator_slot: str) -> dict[str, Any]:
    tier = str(row["authority_tier"])
    return {
        "schema_version": "linkerforge.linker-curation-review-record.v1",
        "review_id": f"{curator_slot}:{row['record_id']}",
        "record_id": str(row["record_id"]),
        "source_dataset_id": "mined_and_natural_linkers_v1",
        "evidence_role": EVIDENCE_ROLES[tier],
        "authority_tier": tier,
        "compatibility_outcome": "not_assessed",
        "source_record": {
            "protein_id": str(row["protein_id"]),
            "linker_sequence": str(row["linker_sequence"]),
            "linker_start_1based": int(row["linker_start_1based"]),
            "linker_end_1based": int(row["linker_end_1based"]),
            "linker_length": int(row["linker_length"]),
            "homology_cluster_id": str(row["homology_cluster_id"]),
            "evidence_grade": str(row["evidence_grade"]),
            "mining_source_corpus": row.get("mining_source_corpus"),
            "benchmark_firewall_verified": bool(row["benchmark_firewall_verified"]),
            "classifier_scores": dict(row.get("classifier_scores", {})),
            "flanking_context": dict(row.get("flanking_context", {})),
            "geometry": row.get("geometry"),
            "biophysical_properties": dict(row.get("biophysical_properties", {})),
        },
        "review": empty_review(),
    }


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def git_provenance(root: Path) -> dict[str, Any]:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, check=True,
            capture_output=True, text=True,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"], cwd=root, check=True,
                capture_output=True, text=True,
            ).stdout.strip()
        )
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = "unavailable", True
    return {"commit_sha": commit, "worktree_modified": dirty}


def build_packet(dataset_dir: Path, config_path: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing review directory: {output_dir}")
    manifest_path = dataset_dir / "manifest.json"
    source_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if source_manifest.get("dataset_id") != config.get("source_dataset_id"):
        raise ValueError("source dataset does not match selection config")

    integrity: dict[str, Any] = {"manifest_sha256": sha256_path(manifest_path), "verified_files": {}}
    for name, info in source_manifest["files"].items():
        path = dataset_dir / str(info["path"])
        actual_hash = sha256_path(path)
        actual_count = sum(1 for line in path.open("r", encoding="utf-8") if line.strip())
        if actual_hash != info["sha256"] or actual_count != int(info["records"]):
            raise ValueError(f"source integrity mismatch: {name}")
        integrity["verified_files"][name] = {
            "path": str(info["path"]),
            "records": actual_count,
            "sha256": actual_hash,
        }

    all_path = dataset_dir / source_manifest["files"]["all_records_jsonl"]["path"]
    records = [json.loads(line) for line in all_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    by_tier: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in records:
        by_tier[str(row["authority_tier"])].append(row)

    quotas = {str(key): int(value) for key, value in config["tier_quotas"].items()}
    target = sum(quotas.values())
    if not int(config["minimum_total_records"]) <= target <= int(config["maximum_total_records"]):
        raise ValueError("configured review total is outside the declared review range")
    unknown = sorted(set(quotas) - EVIDENCE_ROLES.keys())
    if unknown:
        raise ValueError(f"unsupported authority tiers in config: {', '.join(unknown)}")

    seed = str(config["selection_seed"])
    boundaries = [int(value) for value in config["length_bins"]]
    selected: list[dict[str, Any]] = []
    for tier in sorted(quotas):
        selected.extend(
            stratified_select(by_tier[tier], quotas[tier], seed=seed, boundaries=boundaries)
        )
    if len({row["record_id"] for row in selected}) != len(selected):
        raise ValueError("selection contains duplicate record IDs")
    selected.sort(key=lambda row: (str(row["authority_tier"]), stable_rank(seed, str(row["record_id"]))))

    output_dir.mkdir(parents=True)
    curator_a = [review_record(row, "curator_a") for row in selected]
    curator_b = [review_record(row, "curator_b") for row in selected]
    index_rows = [
        {
            "record_id": str(row["record_id"]),
            "authority_tier": str(row["authority_tier"]),
            "evidence_role": EVIDENCE_ROLES[str(row["authority_tier"])],
            "homology_cluster_id": str(row["homology_cluster_id"]),
            "selection_stratum": stratum(row, boundaries),
            "compatibility_outcome": "not_assessed",
        }
        for row in selected
    ]
    packet_paths = {
        "curator_a": output_dir / "curator_a.jsonl",
        "curator_b": output_dir / "curator_b.jsonl",
        "selection_index": output_dir / "selection_index.jsonl",
    }
    write_jsonl(packet_paths["curator_a"], curator_a)
    write_jsonl(packet_paths["curator_b"], curator_b)
    write_jsonl(packet_paths["selection_index"], index_rows)

    root = Path(__file__).resolve().parents[1]
    selected_by_tier = dict(Counter(str(row["authority_tier"]) for row in selected))
    manifest = {
        "schema_version": "linkerforge.linker-curation-review-manifest.v1",
        "status": "awaiting_independent_dual_review",
        "source_dataset_id": "mined_and_natural_linkers_v1",
        "review_target": "linker_identity_and_evidence_not_fusion_compatibility",
        "selection": {
            "seed": seed,
            "config_sha256": sha256_path(config_path),
            "selected_total": len(selected),
            "selected_by_tier": selected_by_tier,
            "selected_family_clusters": len({str(row["homology_cluster_id"]) for row in selected}),
            "selected_strata": len({stratum(row, boundaries) for row in selected}),
        },
        "source_integrity": integrity,
        "files": {
            name: {"path": path.name, "sha256": sha256_path(path), "records": len(selected)}
            for name, path in packet_paths.items()
        },
        "provenance": {
            **git_provenance(root),
            "created_at": datetime.now(timezone.utc).isoformat(),
            "versions": {"python": platform.python_version(), "packet_builder": "1.0.0"},
        },
        "claim_boundary": str(config["claim_boundary"]),
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    (output_dir / "README.md").write_text(
        "# Linker curation review packet\n\n"
        "This packet supports independent review of linker identity, coordinates and source evidence. "
        "It does **not** assign fusion-compatibility or biological-success labels.\n\n"
        "- Give `curator_a.jsonl` and `curator_b.jsonl` to independent reviewers.\n"
        "- Each reviewer fills only the `review` object and preserves every source field.\n"
        "- Do not reveal one review to the other reviewer before both are locked.\n"
        "- Adjudicate disagreements in a separately versioned artifact.\n"
        "- Never relabel classifier predictions as experimental outcomes.\n"
        "- Keep homology clusters intact when a future benchmark is split.\n",
        encoding="utf-8",
        newline="\n",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = build_packet(args.dataset_dir, args.config, args.output_dir)
    print(json.dumps(manifest["selection"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
