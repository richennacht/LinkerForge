"""Read-only status check for the active detector autoresearch round."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import psutil

if __package__:
    from scripts.preflight_caid3_protected_training import verify
else:
    from preflight_caid3_protected_training import verify


DEPRECATED_V1_SPLIT = "data/processed/caid3_protected_dfl_v1/manifest.json"
DEFAULT_V2_SPLIT = "data/processed/dfl_b1_v2_family_split/manifest.json"


def inspect(ledger_path: Path) -> dict:
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    root = ledger_path.resolve().parents[1]
    active = [t for t in ledger.get("trials", []) if t.get("status") == "running"]
    if len(active) > 1:
        return {"action": "integrity_stop", "reason": "multiple running trials in ledger"}
    if not active:
        split_rel = str(ledger.get("split_manifest", DEFAULT_V2_SPLIT)).replace("\\", "/")
        if split_rel == DEPRECATED_V1_SPLIT:
            return {
                "action": "data_integrity_audit_only",
                "round_status": ledger.get("status"),
                "reason": f"Deprecated split {DEPRECATED_V1_SPLIT} is prohibited",
                "next_action": "Use data/processed/dfl_b1_v2_family_split/manifest.json and rerun preflight.",
            }
        protected_manifest = root / split_rel
        try:
            verified = verify(protected_manifest)
        except (ValueError, KeyError, OSError, ImportError) as exc:
            return {
                "action": "data_integrity_audit_only",
                "round_status": ledger.get("status"),
                "split_manifest": split_rel,
                "reason": str(exc),
                "next_action": "Rebuild a versioned benchmark-excluded family split and rerun preflight before any new training job.",
            }
        return {
            "action": "launch_or_audit",
            "round_status": ledger.get("status"),
            "split_manifest": split_rel,
            "preflight_status": verified.get("status"),
            "split_manifest_sha256": verified.get("split_manifest_sha256"),
            "counts": verified.get("counts"),
            "next_action": ledger.get("next_action"),
        }

    trial = active[0]
    pid = trial.get("launcher_pid")
    output_dir = root / trial["output_dir"]
    manifest_path = output_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else None
    process_alive = bool(pid and psutil.pid_exists(pid))
    if process_alive:
        try:
            process = psutil.Process(pid)
            command = " ".join(process.cmdline()).lower()
            expected_script = Path(trial.get("script", "train_lora_linker_detector.py")).name.lower()
            process_alive = process.is_running() and expected_script in command
        except (psutil.Error, OSError):
            process_alive = False
    state = manifest.get("status") if manifest else "missing"
    if process_alive and state == "running":
        action = "monitor_running"
    elif process_alive:
        action = "inspect_process_manifest_mismatch"
    else:
        action = "reconcile_finished"
    return {
        "action": action,
        "trial_id": trial["trial_id"],
        "pid": pid,
        "process_alive": process_alive,
        "manifest_status": state,
        "manifest_path": str(manifest_path),
        "next_action": ledger.get("next_action"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, help="Override the active-round ledger")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.ledger is None:
        pointer = json.loads((root / "configs/detector_autoresearch_active.json").read_text(encoding="utf-8"))
        ledger_path = root / pointer["ledger"]
    else:
        ledger_path = args.ledger
    print(json.dumps(inspect(ledger_path), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
