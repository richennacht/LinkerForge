#!/usr/bin/env python3
"""Prepare matched public 1UBQ/2GB1 fusion controls for ColabFold.

The previously completed 10-aa GGGGS repeat is intentionally excluded from the
batch FASTA and referenced by hash in the manifest.  All other candidates use
the same public domains and are predictions, not biological success labels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


DOMAIN_A = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
DOMAIN_B = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
LINKERS = {
    "direct_0": "",
    "ggggs_5": "GGGGS",
    "ggggs_15_fixed": "GGGGSGGGGSGGGGS",
    "waldo_12": "GSAGSAAGSGEF",
    "bird_18": "KESGSVSSEQLAQFRSLD",
    "eaaak_15_rigid": "EAAAKEAAAKEAAAK",
}
EXISTING = {
    "candidate_id": "ggggs_10",
    "linker": "GGGGSGGGGS",
    "run_id": "acceptance-20261020-ubiquitin-gb1-target10",
    "sequence_sha256": "c690ceebc0a93d73a544b8f24de10572d604059e8c56badd4f28a33b01a1a260",
}


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def git_value(repo: Path, *args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "not-available"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    records = []
    fasta_parts = []
    for candidate_id, linker in LINKERS.items():
        sequence = DOMAIN_A + linker + DOMAIN_B
        record = {
            "candidate_id": candidate_id,
            "linker_sequence": linker,
            "linker_length": len(linker),
            "domain_a_length": len(DOMAIN_A),
            "domain_b_length": len(DOMAIN_B),
            "full_sequence_length": len(sequence),
            "sequence_sha256": sha256_text(sequence),
        }
        records.append(record)
        fasta_parts.append(f">LinkerForge_1UBQ_{candidate_id}_2GB1\n{sequence}\n")
    fasta = "".join(fasta_parts)
    (args.output_dir / "baselines.fasta").write_text(fasta, encoding="utf-8")
    repo = Path(__file__).resolve().parents[1]
    generation_config_sha = sha256_text(json.dumps({
        "domain_a": DOMAIN_A,
        "domain_b": DOMAIN_B,
        "linkers": LINKERS,
        "existing": EXISTING,
    }, sort_keys=True, separators=(",", ":")))
    working_tree_modified = bool(git_value(repo, "status", "--porcelain"))
    manifest = {
        "schema_version": "linkerforge.matched-colabfold-panel.v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "Matched structural sensitivity panel; not experimental validation or ranker training data",
        "evidence_quality": "matched_public_sequence_computational_sensitivity_panel",
        "source_ids": ["PDB:1UBQ:A", "PDB:2GB1:A", f"existing-run:{EXISTING['run_id']}"],
        "public_sources": ["PDB 1UBQ chain A", "PDB 2GB1 chain A"],
        "records": records,
        "existing_matched_candidate": EXISTING,
        "fasta_sha256": sha256_text(fasta),
        "planned_settings": {
            "predictor": "ColabFold/AlphaFold2-ptm",
            "models": 5,
            "seed": 0,
            "templates": False,
            "msa_mode": "mmseqs2_uniref_env",
            "rank_by": "plddt",
            "fast_kernels": True,
        },
        "provenance": {
            "commit_sha": git_value(repo, "rev-parse", "HEAD"),
            "worktree_modified": working_tree_modified,
            "config_sha256": generation_config_sha,
            "versions": {"python": platform.python_version(), "panel_builder": "1.0.0"},
        },
        "claim_boundary": "Predictions are compared under common settings; they do not establish expression, activity, or an optimal linker.",
    }
    (args.output_dir / "panel_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
