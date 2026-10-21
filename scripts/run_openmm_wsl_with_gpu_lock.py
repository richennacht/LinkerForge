#!/usr/bin/env python3
"""Launch OpenMM in WSL under the shared LinkerForge GPU lock."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CREATE_NO_WINDOW = 0x08000000 if os.name == "nt" else 0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        raise ValueError(f"A Windows drive path is required: {resolved}")
    return f"/mnt/{drive}/{resolved.as_posix().split(':', 1)[1].lstrip('/')}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--owner", required=True)
    parser.add_argument("--wsl-python", default="/home/linkerforge/localcolabfold/.pixi/envs/default/bin/python")
    args = parser.parse_args()

    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite production output: {args.output_dir}")
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    if preflight.get("status") != "pass" or preflight.get("run_authorized") is not True:
        raise ValueError("Production preflight did not authorize this run")
    if preflight["inputs"]["production_config_sha256"] != sha256(args.config):
        raise ValueError("Preflight refers to a different production configuration")
    smoke_manifest_path = Path(preflight["inputs"]["smoke_manifest"])
    smoke = json.loads(smoke_manifest_path.read_text(encoding="utf-8"))
    if smoke.get("input_pdb_sha256") != sha256(args.pdb):
        raise ValueError("Production PDB differs from the same-construct smoke/preflight input")

    args.lock.parent.mkdir(parents=True, exist_ok=True)
    lock_contents = (
        f"pid={os.getpid()}\nowner={args.owner}\noutput={args.output_dir.resolve()}\n"
        f"created_utc={utc_now()}\n"
    )
    try:
        with args.lock.open("x", encoding="utf-8") as handle:
            handle.write(lock_contents)
    except FileExistsError as exc:
        raise RuntimeError(
            f"Another LinkerForge GPU job may be running. Inspect {args.lock} and remove it only after confirming no prediction or simulation process is active."
        ) from exc

    command = [
        "wsl.exe", "--exec", args.wsl_python,
        wsl_path(ROOT / "scripts/run_openmm_validation.py"),
        "--pdb", wsl_path(args.pdb),
        "--config", wsl_path(args.config),
        "--output-dir", wsl_path(args.output_dir),
    ]
    started = utc_now()
    exit_code: int | None = None
    error: str | None = None
    try:
        completed = subprocess.run(command, cwd=ROOT, check=False, creationflags=CREATE_NO_WINDOW)
        exit_code = completed.returncode
        if exit_code != 0:
            raise RuntimeError(f"OpenMM WSL process exited {exit_code}")
        manifest_path = args.output_dir / "run_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("status") != "complete" or manifest.get("input_pdb_sha256") != sha256(args.pdb):
            raise RuntimeError("OpenMM process returned without a matching complete manifest")
    except Exception as exc:
        error = str(exc)
        raise
    finally:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        launcher = {
            "schema_version": "linkerforge.openmm-wsl-launch.v1",
            "status": "complete" if error is None and exit_code == 0 else "failed",
            "started_utc": started,
            "completed_utc": utc_now(),
            "command": command,
            "exit_code": exit_code,
            "error": error,
            "inputs": {
                "pdb": str(args.pdb.resolve()), "pdb_sha256": sha256(args.pdb),
                "config": str(args.config.resolve()), "config_sha256": sha256(args.config),
                "preflight": str(args.preflight.resolve()), "preflight_sha256": sha256(args.preflight),
            },
            "provenance": {
                "commit_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "script_sha256": sha256(Path(__file__)),
                "python": platform.python_version(),
            },
        }
        (args.output_dir / "wsl_launcher.json").write_text(json.dumps(launcher, indent=2) + "\n", encoding="utf-8")
        if args.lock.exists() and args.lock.read_text(encoding="utf-8") == lock_contents:
            args.lock.unlink()


if __name__ == "__main__":
    main()
