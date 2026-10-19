"""Label-agnostic B1 V2 scoring; auxiliary scoring requires a frozen selection.

This is a development-equivalence and order-only specificity diagnostic, not
an independent 99%-precision certificate or a proteome mining entry point.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
if __package__:
    from scripts.audit_dfl_b1_v2_precision import precision_frontier
    from scripts.train_dfl_b1_v2_three_layer_cnn import CONFIG, preflight
else:
    from audit_dfl_b1_v2_precision import precision_frontier
    from train_dfl_b1_v2_three_layer_cnn import CONFIG, preflight


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def all_windows(rows: list[dict], width: int):
    """Window every residue; labels must not influence sequence inclusion/order."""
    return [(row["sequence"][start:start + width],
             row["labels"][start:start + width], row["protein_id"])
            for row in rows for start in range(0, len(row["sequence"]), width)]


def score(run_dir: Path, role: str, selection: Path | None = None) -> dict:
    import numpy as np
    import torch
    from sklearn.metrics import average_precision_score, roc_auc_score
    from transformers import AutoModel, AutoTokenizer
    from linkerforge.detector.cnn import ThreeLayerDFLCNN

    if role not in {"development", "auxiliary"}:
        raise ValueError("Only development or auxiliary is available")
    run_dir = run_dir.resolve()
    manifest_path = run_dir / "manifest.json"
    checkpoint_path = run_dir / "best.pt"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["status"] != "completed":
        raise ValueError("Only completed runs can be scored")
    config_path = ROOT / manifest.get("config_path", str(CONFIG.relative_to(ROOT)))
    config, _, _, dev = preflight(config_path)
    if sha(config_path) != manifest["config_sha256"]:
        raise ValueError("Config hash mismatch")
    if manifest["split_manifest_sha256"] != config["split_manifest_sha256"]:
        raise ValueError("Split hash mismatch")
    if role == "auxiliary":
        if selection is None:
            raise ValueError("Auxiliary labels require a frozen selection manifest")
        selected = json.loads(selection.read_text(encoding="utf-8"))
        if selected.get("status") != "frozen_before_auxiliary" or selected.get("evaluator_sha256") != sha(Path(__file__)):
            raise ValueError("Selection state or evaluator hash mismatch")
        if selected.get("checkpoint_sha256") != sha(checkpoint_path) or selected.get("run_manifest_sha256") != sha(manifest_path):
            raise ValueError("Frozen checkpoint or run manifest mismatch")
        if selected.get("config_sha256") != manifest["config_sha256"]:
            raise ValueError("Frozen config mismatch")
        if selected.get("split_manifest_sha256") != config["split_manifest_sha256"]:
            raise ValueError("Frozen split mismatch")
        if selected.get("auxiliary_sha256") != config["auxiliary_sha256"]:
            raise ValueError("Frozen auxiliary panel hash mismatch")
        dev_audit_path = ROOT / selected["development_audit_path"]
        if not dev_audit_path.resolve().is_relative_to(ROOT / "artifacts") or sha(dev_audit_path) != selected["development_audit_sha256"]:
            raise ValueError("Frozen development audit missing or changed")
        threshold = selected["threshold"]
        if not isinstance(threshold, (int, float)) or not 0 <= threshold <= 1:
            raise ValueError("Frozen threshold invalid")
        base = ROOT / config["split_manifest"]
        rows = [json.loads(line) for line in (base.parent / "auxiliary.jsonl").read_text(encoding="utf-8").splitlines()]
        if any("1" in row["labels"] for row in rows):
            raise ValueError("Expected a negative-only auxiliary panel")
    else:
        rows = dev
        threshold = None

    torch.set_num_threads(min(4, torch.get_num_threads()))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    if checkpoint["config_sha256"] != manifest["config_sha256"]:
        raise ValueError("Checkpoint config mismatch")
    tokenizer = AutoTokenizer.from_pretrained(config["encoder"], local_files_only=True)
    backbone = AutoModel.from_pretrained(config["encoder"], local_files_only=True).to(device).eval()
    backbone.requires_grad_(False)
    head = ThreeLayerDFLCNN(backbone.config.hidden_size).to(device).eval()
    head.load_state_dict(checkpoint["head"])
    labels, scores, unknown_scores = [], [], []
    protein_calls = {}
    for sequence, mask, protein_id in all_windows(rows, config["max_residues_per_window"]):
        tokens = tokenizer([sequence], return_tensors="pt", add_special_tokens=True)
        tokens = {key: value.to(device) for key, value in tokens.items()}
        with torch.no_grad():
            features = backbone(**tokens).last_hidden_state[:, 1:1 + len(sequence)]
            values = head(features).sigmoid().squeeze(0).cpu().numpy()
        if len(values) != len(mask):
            raise ValueError("Inference did not cover every residue")
        for label, value in zip(mask, values, strict=True):
            if label == "x":
                unknown_scores.append(float(value))
            elif label in "01":
                labels.append(int(label))
                scores.append(float(value))
                if threshold is not None and value >= threshold:
                    protein_calls[protein_id] = protein_calls.get(protein_id, 0) + 1
            else:
                raise ValueError("Unexpected label symbol")
    result = {
        "schema_version": "1.0.0", "role": role,
        "claim_status": "deployment_scoring_diagnostic_not_mining_certified",
        "scoring_policy": "all_windows_one_at_a_time_without_using_labels_for_inclusion",
        "device": str(device), "seed": manifest["seed"],
        "run_manifest_sha256": sha(manifest_path), "checkpoint_sha256": sha(checkpoint_path),
        "config_sha256": manifest["config_sha256"],
        "split_manifest_sha256": manifest["split_manifest_sha256"],
        "evaluator_sha256": sha(Path(__file__)),
        "proteins": len(rows), "known_positive_residues": int(sum(labels)),
        "known_negative_residues": len(labels) - int(sum(labels)),
        "unknown_residues": len(unknown_scores),
        "total_residues": len(labels) + len(unknown_scores),
    }
    if role == "development":
        if set(labels) != {0, 1}:
            raise ValueError("Development requires both labelled classes")
        result["aps"] = float(average_precision_score(labels, scores))
        result["auroc"] = float(roc_auc_score(labels, scores))
        result["precision_99_min50tp"] = precision_frontier(labels, scores)
        point = result["precision_99_min50tp"]
        if point["attained_on_known_development_labels"]:
            result["unknown_residues_called_at_development_threshold"] = int(
                np.sum(np.asarray(unknown_scores) >= point["threshold"]))
    else:
        result["frozen_threshold"] = float(threshold)
        result["false_positive_order_residues"] = int(np.sum(np.asarray(scores) >= threshold))
        result["known_order_residue_specificity"] = float(1 - result["false_positive_order_residues"] / len(scores))
        result["unknown_residues_called"] = int(np.sum(np.asarray(unknown_scores) >= threshold))
        result["proteins_with_known_order_calls"] = len(protein_calls)
        result["precision_estimable"] = False
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--role", choices=("development", "auxiliary"), required=True)
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = score(args.run_dir, args.role, args.selection)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
