"""Post-hoc CAID3 common-27 comparison of a development-frozen B1 V2 CNN.

CAID3 has already been viewed in this project. Passing this lock does not make
the test untouched, prospective, or suitable for model selection.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import torch
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score
from transformers import AutoModel, AutoTokenizer

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from audit_caid3_dld_overlap import caid_records
from compare_dfl_caid3_common import competitor_scores
from preflight_dfl_family_training_v2 import verify as verify_v2
from linkerforge.detector.cnn import ThreeLayerDFLCNN

REFERENCE = ROOT / "artifacts/caid3-common-coverage-001/linker-common.fasta"
BASELINES = ROOT / "artifacts/caid3-common-coverage-001/predictions"
EVALUATOR = ROOT / "artifacts/caid-evaluator/caid.py"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def known_metrics(reference, by_id: dict[str, list[float]]) -> dict:
    """Local cross-check; the official CAID evaluator remains authoritative."""
    pairs = [(label, value)
             for identifier, sequence, labels in reference
             for label, value in zip(labels, by_id[identifier], strict=True)
             if label in "01"]
    if not pairs:
        raise ValueError("No known CAID3 residues")
    truth = np.asarray([label == "1" for label, _ in pairs], dtype=np.uint8)
    scores = np.asarray([value for _, value in pairs], dtype=np.float64)
    if not np.isfinite(scores).all() or set(truth.tolist()) != {0, 1}:
        raise ValueError("Invalid known-label score vector")
    precision, recall, _ = precision_recall_curve(truth, scores)
    f1 = np.divide(2 * precision * recall, precision + recall,
                   out=np.zeros_like(precision), where=precision + recall > 0)
    return {"proteins": len(reference), "known_residues": len(pairs),
            "positive_residues": int(truth.sum()),
            "aps": float(average_precision_score(truth, scores)),
            "auroc": float(roc_auc_score(truth, scores)), "fmax": float(f1.max())}


def validate_selection(selection_path: Path) -> tuple[dict, dict, dict, Path]:
    """Validate frozen choices before opening the CAID3 reference."""
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    if (selection.get("status") != "frozen_before_caid3_b1_v2_retrospective"
            or selection.get("selection_basis") != "b1_v2_development_only"
            or selection.get("claim_status") != "post_hoc_retrospective_not_untouched"):
        raise ValueError("This comparison requires a development-only retrospective selection lock")
    if selection.get("scorer_sha256") != sha(Path(__file__)):
        raise ValueError("Scorer changed after selection freeze")
    checkpoint_path = (ROOT / selection["checkpoint"]).resolve()
    if not checkpoint_path.is_relative_to(ROOT / "artifacts") or sha(checkpoint_path) != selection["checkpoint_sha256"]:
        raise ValueError("Frozen checkpoint missing or changed")
    manifest_path = checkpoint_path.parent / "manifest.json"
    if sha(manifest_path) != selection["run_manifest_sha256"]:
        raise ValueError("Run manifest changed after freeze")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "completed" or manifest.get("seed") != selection.get("seed"):
        raise ValueError("Frozen run incomplete or seed mismatch")
    config_path = (ROOT / selection["config"]).resolve()
    if sha(config_path) != selection["config_sha256"]:
        raise ValueError("Frozen configuration changed")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if (manifest.get("config_sha256") != selection["config_sha256"]
            or manifest.get("split_manifest_sha256") != selection["split_manifest_sha256"]
            or config.get("split_manifest_sha256") != selection["split_manifest_sha256"]):
        raise ValueError("Frozen training lineage mismatch")
    audit_path = (ROOT / selection["development_audit"]).resolve()
    if not audit_path.is_relative_to(ROOT / "artifacts") or sha(audit_path) != selection["development_audit_sha256"]:
        raise ValueError("Development-only audit changed after freeze")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (audit.get("role") != "development"
            or audit.get("checkpoint_sha256") != selection["checkpoint_sha256"]
            or audit.get("split_manifest_sha256") != selection["split_manifest_sha256"]):
        raise ValueError("Development-only audit lineage mismatch")
    split_path = ROOT / config["split_manifest"]
    if sha(split_path) != selection["split_manifest_sha256"]:
        raise ValueError("Protected split changed")
    verify_v2(split_path)
    if config.get("encoder") != selection.get("encoder") or config.get("max_residues_per_window") != selection.get("window"):
        raise ValueError("Frozen inference policy mismatch")
    if selection.get("reference_sha256") != sha(REFERENCE) or selection.get("official_evaluator_sha256") != sha(EVALUATOR):
        raise ValueError("Reference or official evaluator changed after freeze")
    baselines = sorted(BASELINES.glob("*.caid"))
    if len(baselines) != 7 or {p.name: sha(p) for p in baselines} != selection.get("baseline_sha256"):
        raise ValueError("Common-target competitor files changed after freeze")
    return selection, manifest, config, checkpoint_path


@torch.no_grad()
def predict(backbone, tokenizer, head, sequence: str, window: int, device: torch.device) -> list[float]:
    values: list[float] = []
    for start in range(0, len(sequence), window):
        fragment = sequence[start:start + window]
        tokens = tokenizer([fragment], return_tensors="pt", add_special_tokens=True, truncation=False)
        hidden = backbone(**{key: value.to(device) for key, value in tokens.items()}).last_hidden_state[:, 1:1 + len(fragment)]
        part = head(hidden).sigmoid().squeeze(0).cpu().tolist()
        if len(part) != len(fragment):
            raise ValueError("Prediction length mismatch")
        values.extend(part)
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.selection = args.selection.resolve()
    args.output_dir = args.output_dir.resolve()
    if args.output_dir.exists():
        raise FileExistsError("Retrospective output already exists; no repeated overwrite")
    selection, manifest, config, checkpoint_path = validate_selection(args.selection)
    # All selection and input hashes are checked before labels enter memory.
    reference = caid_records(REFERENCE)
    if len(reference) != 27:
        raise ValueError("Unexpected common-target reference size")
    baseline_paths = sorted(BASELINES.glob("*.caid"))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.set_num_threads(min(4, torch.get_num_threads()))
    tokenizer = AutoTokenizer.from_pretrained(config["encoder"], local_files_only=True)
    backbone = AutoModel.from_pretrained(config["encoder"], local_files_only=True).to(device).eval()
    backbone.requires_grad_(False)
    head = ThreeLayerDFLCNN(backbone.config.hidden_size).to(device).eval()
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    if checkpoint.get("config_sha256") != manifest["config_sha256"]:
        raise ValueError("Checkpoint/config mismatch")
    head.load_state_dict(checkpoint["head"])
    by_id = {identifier: predict(backbone, tokenizer, head, sequence, selection["window"], device)
             for identifier, sequence, _ in reference}
    results = {path.stem: known_metrics(reference, competitor_scores(path, reference)) for path in baseline_paths}
    results["LinkerForge_B1_V2_frozen"] = known_metrics(reference, by_id)
    args.output_dir.mkdir(parents=True)
    predictions = args.output_dir / "predictions"
    predictions.mkdir()
    for path in baseline_paths:
        shutil.copyfile(path, predictions / path.name)
    ours = predictions / "LinkerForge_B1_V2_frozen.caid"
    with ours.open("w", encoding="utf-8") as stream:
        for identifier, sequence, _ in reference:
            stream.write(f">{identifier}\n")
            for position, (residue, score) in enumerate(zip(sequence, by_id[identifier]), 1):
                stream.write(f"{position}\t{residue}\t{score:.8f}\t{int(score >= selection['display_threshold'])}\n")
    official = args.output_dir / "official_evaluator"
    subprocess.run([sys.executable, str(EVALUATOR), str(REFERENCE), str(predictions), "-o", str(official)],
                   cwd=EVALUATOR.parent, check=True, timeout=3600)
    report = {
        "schema_version": "1.0.0", "status": "post_hoc_public_retrospective_not_untouched",
        "created_utc": datetime.now(timezone.utc).isoformat(), "device": str(device),
        "reference_sha256": sha(REFERENCE), "selection_sha256": sha(args.selection),
        "checkpoint_sha256": sha(checkpoint_path), "run_manifest_sha256": sha(checkpoint_path.parent / "manifest.json"),
        "split_manifest_sha256": selection["split_manifest_sha256"], "scorer_sha256": sha(Path(__file__)),
        "official_evaluator_sha256": sha(EVALUATOR), "prediction_sha256": sha(ours),
        "baseline_sha256": {p.name: sha(p) for p in baseline_paths},
        "results_local_crosscheck": results, "official_evaluator_dir": str(official),
        "limitations": ["CAID3 scores were viewed before this run; not a fresh independent test.",
                        "CAID3 common-27 lacks verified no-linker proteins for specificity.",
                        "The B1 V2 firewall excludes declared sequence/family matches; it cannot audit all ESM pretraining or unknown remote homology."],
    }
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
