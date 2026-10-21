#!/usr/bin/env python3
"""Summarize the matched 1UBQ/2GB1 ColabFold baseline panel.

No composite winner is selected because acceptance thresholds and weights have
not been frozen on development data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_fusion_structure_evidence import aligned_rmsd, coordinate_audit, read_ca  # noqa: E402


ALIASES = {"ggggs_15_fixed": "ggggs_15", "eaaak_15_rigid": "eaaak_15"}
MATCHED_CONFIG_KEYS = [
    "use_templates",
    "msa_mode",
    "model_type",
    "num_models",
    "rank_by",
    "random_seed",
    "num_seeds",
    "use_fast_kernels",
    "version",
    "commit",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def matching_file(directory: Path, parts: list[str], suffix: str) -> Path:
    matches = [path for path in directory.glob(f"*{suffix}") if all(part in path.name for part in parts)]
    if len(matches) != 1:
        raise ValueError(f"Expected one {parts} {suffix} file in {directory}; found {len(matches)}")
    return matches[0]


def metrics(candidate_id: str, linker_length: int, prediction_dir: Path, domain_a_reference: Path, domain_b_reference: Path, prefix_override: str | None = None) -> dict:
    prefix = prefix_override or f"LinkerForge_1UBQ_{candidate_id}_2GB1"
    pdb_path = matching_file(prediction_dir, [prefix, "unrelaxed_rank_001"], ".pdb")
    score_path = matching_file(prediction_dir, [prefix, "scores_rank_001"], ".json")
    pdb = read_ca(pdb_path)
    score = json.loads(score_path.read_text(encoding="utf-8"))
    plddt = [float(value) for value in score["plddt"]]
    domain_a_length = 76
    domain_b_start = domain_a_length + linker_length
    pae = score["pae"]
    cross = [
        float(pae[left][right])
        for left in range(domain_a_length)
        for right in range(domain_b_start, len(pdb))
    ] + [
        float(pae[left][right])
        for left in range(domain_b_start, len(pdb))
        for right in range(domain_a_length)
    ]
    reference_a = read_ca(domain_a_reference, "A")
    reference_b = read_ca(domain_b_reference, "A")
    return {
        "candidate_id": candidate_id,
        "linker_length": linker_length,
        "mean_plddt": round(mean(plddt) or 0, 4),
        "domain_a_mean_plddt": round(mean(plddt[:domain_a_length]) or 0, 4),
        "linker_mean_plddt": round(mean(plddt[domain_a_length:domain_b_start]) or 0, 4) if linker_length else None,
        "domain_b_mean_plddt": round(mean(plddt[domain_b_start:]) or 0, 4),
        "ptm": score.get("ptm"),
        "mean_cross_domain_pae_angstrom": round(mean(cross) or 0, 4),
        "domain_a_rmsd_angstrom": round(aligned_rmsd(pdb[:domain_a_length], reference_a), 4),
        "domain_b_rmsd_angstrom": round(aligned_rmsd(pdb[domain_b_start:], reference_b), 4),
        "coordinate_pre_screen": coordinate_audit(pdb),
        "pdb": str(pdb_path),
        "pdb_sha256": sha256(pdb_path),
        "score_json": str(score_path),
        "score_json_sha256": sha256(score_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", type=Path, required=True)
    parser.add_argument("--existing-run-dir", type=Path, required=True)
    parser.add_argument("--dflpred-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.panel_dir / "panel_manifest.json").read_text(encoding="utf-8"))
    panel_config_path = args.panel_dir / "prediction" / "config.json"
    existing_config_path = args.existing_run_dir / "prediction" / "config.json"
    panel_config = json.loads(panel_config_path.read_text(encoding="utf-8"))
    existing_config = json.loads(existing_config_path.read_text(encoding="utf-8"))
    mismatches = {
        key: {"panel": panel_config.get(key), "existing": existing_config.get(key)}
        for key in MATCHED_CONFIG_KEYS
        if panel_config.get(key) != existing_config.get(key)
    }
    if mismatches:
        raise ValueError(f"Panel and existing GS10 run do not use matched predictor settings: {mismatches}")
    reference_a = args.existing_run_dir / "references" / "1UBQ.pdb"
    reference_b = args.existing_run_dir / "references" / "2GB1.pdb"
    rows = [
        metrics(record["candidate_id"], record["linker_length"], args.panel_dir / "prediction", reference_a, reference_b)
        for record in manifest["records"]
    ]
    existing = metrics(
        "ggggs_10", 10, args.existing_run_dir / "prediction", reference_a, reference_b,
        prefix_override="LinkerForge_acceptance_1UBQ_GGGGSGGGGS_2GB1",
    )
    rows.append(existing)
    detector = json.loads(args.dflpred_summary.read_text(encoding="utf-8"))
    detector_rows = {(row["domain_pair"], row["linker_id"]): row for row in detector["results"]}
    for row in rows:
        detector_id = ALIASES.get(row["candidate_id"], row["candidate_id"])
        dfl = detector_rows.get(("1UBQ_2GB1", detector_id))
        row["dflpred"] = None if dfl is None else {
            "linker_residues_called": dfl["linker_residues_called"],
            "linker_recall": dfl["linker_recall"],
            "domain_call_fraction": dfl["domain_call_fraction"],
        }
    analysis_config = {
        "domain_a_length": 76,
        "domain_b_length": 56,
        "matched_predictor_keys": MATCHED_CONFIG_KEYS,
        "aliases": ALIASES,
        "selection_policy": "no_winner_selected_thresholds_not_frozen",
    }
    analysis_config_sha = hashlib.sha256(
        json.dumps(analysis_config, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    commit_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    worktree_modified = bool(subprocess.run(
        ["git", "status", "--porcelain"], cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip())
    result = {
        "schema_version": "linkerforge.matched-structure-comparison.v1",
        "status": "complete",
        "evidence_quality": "matched_single_public_domain_pair_computational_sensitivity",
        "source_ids": [
            f"panel-manifest-sha256:{sha256(args.panel_dir / 'panel_manifest.json')}",
            f"existing-run-manifest-sha256:{sha256(args.existing_run_dir / 'run_manifest.json')}",
            f"dflpred-summary-sha256:{sha256(args.dflpred_summary)}",
        ],
        "provenance": {
            "commit_sha": commit_sha,
            "worktree_modified": worktree_modified,
            "config_sha256": analysis_config_sha,
            "versions": {
                "python": platform.python_version(),
                "colabfold": panel_config.get("version", "unavailable"),
                "colabfold_commit": panel_config.get("commit", "unavailable"),
            },
        },
        "analysis_config": analysis_config,
        "panel_manifest_sha256": sha256(args.panel_dir / "panel_manifest.json"),
        "existing_run_manifest_sha256": sha256(args.existing_run_dir / "run_manifest.json"),
        "panel_config_sha256": sha256(panel_config_path),
        "existing_config_sha256": sha256(existing_config_path),
        "dflpred_summary_sha256": sha256(args.dflpred_summary),
        "predictor_settings": {key: panel_config.get(key) for key in MATCHED_CONFIG_KEYS},
        "predictor_settings_match": True,
        "records": sorted(rows, key=lambda row: row["candidate_id"]),
        "selection_status": "no_winner_selected_thresholds_not_frozen",
        "limitations": [
            "One public domain pair cannot estimate generalization accuracy.",
            "AlphaFold-family confidence is not experimental fusion success.",
            "Each candidate uses one seed; stochastic sensitivity remains unmeasured.",
            "DFLpred detects flexible-linker propensity and is not a compatibility score.",
        ],
        "claim_boundary": "Matched computational sensitivity comparison only; no optimal-linker claim.",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
