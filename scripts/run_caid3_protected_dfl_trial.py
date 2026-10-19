"""Run a DFL baseline only after the V2 family-disjoint CAID3/DP94 firewall passes.

Reuses the previously frozen SAM-like/YOLO-1D training function with explicit
versioned V2 family-disjoint data paths (`data/processed/dfl_b1_v2_family_split`).
The deprecated `caid3_protected_dfl_v1` split is explicitly prohibited.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

if __package__:
    from scripts.check_detector_autoresearch import DEPRECATED_V1_SPLIT, inspect as inspect_autoresearch
    from scripts.preflight_caid3_protected_training import verify
else:
    from check_detector_autoresearch import DEPRECATED_V1_SPLIT, inspect as inspect_autoresearch
    from preflight_caid3_protected_training import verify

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/dfl_b1_v2_family_split"
CONFIG = ROOT / "configs/dfl_b1_v2_family_round_2026-10-18.json"
ACTIVE_POINTER = ROOT / "configs/detector_autoresearch_active.json"
ALIAS = DATA / "auxiliary_alias_manifest.json"
FROZEN_TRAINING_SOURCE = ROOT / "scripts/run_dfl_disorder_transfer.py"
FROZEN_TRAINING_SOURCE_SHA = "43c529dd72cc42f4ff1a5ade38a5d5601a6e4c7ecdce6501c162219a8128a4d4"

VALID_PREFLIGHT_STATUSES = frozenset({
    "pass_v2_preregistered_sequence_and_family_firewall",
    "pass_under_declared_sequence_screen",
})


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha_canonical_lf(path: Path) -> str:
    canonical = path.read_text(encoding="utf-8").replace("\r\n", "\n").encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def preflight() -> tuple[list[dict], list[dict]]:
    if FROZEN_TRAINING_SOURCE_SHA not in {sha(FROZEN_TRAINING_SOURCE), sha_canonical_lf(FROZEN_TRAINING_SOURCE)}:
        raise ValueError("Original training implementation changed; review before running")
    manifest_rel = str((DATA / "manifest.json").relative_to(ROOT)).replace("\\", "/")
    if manifest_rel == DEPRECATED_V1_SPLIT:
        raise ValueError(f"Deprecated split {DEPRECATED_V1_SPLIT} is prohibited")
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    if config["data_manifest"] != manifest_rel:
        raise ValueError("Training config does not name the protected data manifest")
    verified = verify(DATA / "manifest.json")
    if verified["status"] not in VALID_PREFLIGHT_STATUSES:
        raise ValueError("CAID3/DP94 family firewall did not pass")
    if ACTIVE_POINTER.exists():
        pointer = json.loads(ACTIVE_POINTER.read_text(encoding="utf-8"))
        if pointer.get("split_manifest") != manifest_rel:
            raise ValueError("Active autoresearch pointer does not match V2 split manifest")
    alias = json.loads(ALIAS.read_text(encoding="utf-8"))
    if alias["files"]["pretrain"]["sha256"] != sha(DATA / "auxiliary.jsonl"):
        raise ValueError("Auxiliary alias hash mismatch")
    train = [json.loads(line) for line in (DATA / "train.jsonl").read_text(encoding="utf-8").splitlines()]
    development = [json.loads(line) for line in (DATA / "development.jsonl").read_text(encoding="utf-8").splitlines()]
    return train, development


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arch", choices=("samlike", "yolo1d"), required=True)
    parser.add_argument("--condition", choices=("direct", "disorder_transfer"), required=True)
    parser.add_argument("--seed", type=int, choices=(7, 17, 31), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-seconds", type=float, default=1800)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    preflight()
    if args.preflight_only:
        print(json.dumps({
            "status": "preflight_passed",
            "data_manifest": str((DATA / "manifest.json").relative_to(ROOT)).replace("\\", "/"),
            "data_manifest_sha256": sha(DATA / "manifest.json"),
        }, indent=2, sort_keys=True))
        return
    import run_dfl_disorder_transfer as trial

    trial.CONFIG_PATH = CONFIG
    trial.TRAIN = DATA / "train.jsonl"
    trial.DEV = DATA / "development.jsonl"
    trial.AUX = DATA / "auxiliary.jsonl"
    trial.AUX_MANIFEST = ALIAS
    trial.DATA_MANIFEST = DATA / "manifest.json"
    trial.preflight = preflight
    try:
        trial.run_trial(args)
    finally:
        manifest_path = args.output_dir / "manifest.json"
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["launcher_sha256"] = sha(Path(__file__))
            manifest["blacklist_manifest_sha256"] = sha(DATA / "manifest.json")
            manifest["blacklist_status"] = "passed_v2_family_firewall_before_training"
            manifest["claim_status"] = "V2_family_disjoint_development_only"
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
