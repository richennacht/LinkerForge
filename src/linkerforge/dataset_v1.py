"""Loader and query interface for the Part 1 -> Part 2 Mined and Natural Linker Dataset (`v1`).

Provides provenance-stratified access to:
  1. `natural_cath_pdb_verified`: 91 coordinate-complete CATH v4.4.0 natural two-domain
     inter-domain linkers with 3D attachment geometry (`inter_port_ca_distance_A`,
     `physical_3p8A_feasible`), flanking domain sequences/CATH families, Seed-31 CNN
     scores, and `linkerforge.properties` biophysical calls.
  2. `natural_disprot_curated_dfl`: Benchmark-disjoint (`CAID3`/`DP94`-screened)
     curator-validated `DisProt_2026_06` `IDPO:0000033` natural flexible linker spans
     with N/C-terminal flanking sequences and domain architecture.
  3. `mined_domain_gated_interdomain`: Benchmark-disjoint hetero-domain (`left != right`)
     inter-domain spans (`5..60` aa) extracted via `domain_gap_candidates` and scored by
     the frozen Seed-31 `ThreeLayerDFLCNN` (`c0d356c55ed83807353e66b875501c61b17f0b9b83f2910c0f1147f6409fdf6f`).
  4. `literature_engineered_template`: Cited Chen et al. (2013) engineered linker
     templates with rigidity/cleavability classifications.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable

DEFAULT_DATASET_DIR = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "processed"
    / "mined_and_natural_linkers_v1"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_dataset_manifest(dataset_dir: Path = DEFAULT_DATASET_DIR) -> dict:
    manifest_path = Path(dataset_dir) / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for key, info in manifest["files"].items():
        file_path = Path(dataset_dir) / info["path"]
        if _sha256(file_path) != info["sha256"]:
            raise ValueError(f"Dataset file hash mismatch for {key}: {file_path}")
    return manifest


def load_linker_records(
    dataset_dir: Path = DEFAULT_DATASET_DIR,
    *,
    tiers: Iterable[str] | None = None,
    only_feasible_3p8A: bool = False,
    min_length: int = 5,
    max_length: int = 60,
) -> list[dict]:
    """Load hash-verified linker records filtered by authority tier and length."""
    dataset_dir = Path(dataset_dir)
    manifest = load_dataset_manifest(dataset_dir)
    all_path = dataset_dir / manifest["files"]["all_records_jsonl"]["path"]
    allowed_tiers = set(tiers) if tiers is not None else None
    records: list[dict] = []
    for line in all_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if allowed_tiers is not None and row["authority_tier"] not in allowed_tiers:
            continue
        if not (min_length <= int(row["linker_length"]) <= max_length):
            continue
        if only_feasible_3p8A:
            geom = row.get("geometry")
            if geom and geom.get("physical_3p8A_feasible") is False:
                continue
        records.append(row)
    return records
