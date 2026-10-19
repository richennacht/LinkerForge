"""Run a bounded ColabFold prediction for one complete fusion FASTA.

This is an optional compute adapter, not an installation of ColabFold. It never
models the linker alone: the submitted sequence is Domain A + linker + Domain B.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

AA3 = dict(zip(
    "ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCEQGHILKMFPSTWYV",
))
CANONICAL = set("ACDEFGHIKLMNPQRSTVWY")


def read_one_fasta(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    headers = [line for line in lines if line.startswith(">")]
    if len(headers) != 1:
        raise ValueError("Provide exactly one full-fusion FASTA entry")
    sequence = "".join(line.strip().upper() for line in lines if line and not line.startswith(">"))
    if not sequence or len(sequence) > 2100 or set(sequence) - CANONICAL:
        raise ValueError("Full fusion must have 1-2100 standard amino-acid residues")
    return sequence


def pdb_ca_sequence(path: Path) -> str:
    residues: dict[tuple[str, str], str] = {}
    chains: set[str] = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM  ") or line[12:16].strip() != "CA":
            continue
        chain = line[21:22]
        residue = AA3.get(line[17:20].strip())
        if residue is None:
            raise ValueError(f"Unknown residue in {path.name}")
        chains.add(chain)
        residues[(chain, line[22:27])] = residue
    if len(chains) != 1 or not residues:
        raise ValueError("Prediction must contain one full-fusion chain")
    return "".join(residues.values())


def run(fasta: Path, output: Path, executable: str, timeout_minutes: int) -> Path:
    sequence = read_one_fasta(fasta)
    resolved = shutil.which(executable)
    if not resolved:
        raise FileNotFoundError(f"ColabFold executable unavailable: {executable}")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output}")
    output.mkdir(parents=True, exist_ok=True)
    command = [resolved, str(fasta.resolve()), str(output.resolve())]
    started = datetime.now(timezone.utc).isoformat()
    configuration = {"executable": str(Path(resolved).resolve()), "timeout_minutes": timeout_minutes}
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unavailable"
    try:
        predictor_version = importlib.metadata.version("colabfold")
    except importlib.metadata.PackageNotFoundError:
        predictor_version = "not-resolved-from-current-python"
    log_path = output / "colabfold.log"
    manifest_path = output / "run_manifest.json"
    manifest = {
        "status": "running",
        "started_utc": started,
        "input_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "sequence_length": len(sequence),
        "source_commit": commit,
        "config_sha256": hashlib.sha256(json.dumps(configuration, sort_keys=True).encode()).hexdigest(),
        "adapter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "command": command,
        "timeout_minutes": timeout_minutes,
        "predictor": "ColabFold full-fusion sequence",
        "predictor_package_version": predictor_version,
        "evidence_quality": "predicted_unvalidated",
        "claim": "predicted structure, not experimental validation",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    try:
        with log_path.open("w", encoding="utf-8") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=timeout_minutes * 60, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"ColabFold exited {result.returncode}; see {log_path}")
        matches = [path for path in output.rglob("*.pdb") if pdb_ca_sequence(path) == sequence]
        if not matches:
            raise RuntimeError("No single-chain PDB matched the complete submitted fusion sequence")
        chosen = sorted(matches, key=lambda path: ("rank_001" not in path.name, path.name))[0]
        manifest.update(status="complete", selected_pdb=str(chosen.resolve()), selected_pdb_sha256=hashlib.sha256(chosen.read_bytes()).hexdigest())
        return chosen
    except Exception as exc:
        manifest.update(status="failed", failure_reason=str(exc))
        raise
    finally:
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--colabfold", default="colabfold_batch")
    parser.add_argument("--timeout-minutes", type=int, default=120)
    args = parser.parse_args()
    if not 1 <= args.timeout_minutes <= 360:
        parser.error("--timeout-minutes must be between 1 and 360")
    try:
        print(run(args.fasta, args.output, args.colabfold, args.timeout_minutes))
    except Exception as exc:
        print(f"Full-fusion prediction unavailable: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
