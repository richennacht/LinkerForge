"""Audit tied-score high-precision operating points on known B1 V2 development labels.

This cannot certify proteome-mining precision: x residues are unknown, not
negative, and the small curated-order panel is not a representative sample.
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
    from scripts.train_dfl_b1_v2_three_layer_cnn import CONFIG, preflight, windows
else:
    from train_dfl_b1_v2_three_layer_cnn import CONFIG, preflight, windows


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def precision_frontier(truth, scores, *, target: float = 0.99, minimum_tp: int = 50):
    """Return the maximum-recall qualifying point; never split equal scores."""
    import numpy as np
    from scipy.stats import beta

    truth = np.asarray(truth, dtype=np.uint8)
    scores = np.asarray(scores, dtype=np.float64)
    if len(truth) != len(scores) or not len(truth) or set(truth.tolist()) != {0, 1}:
        raise ValueError("Expected aligned binary labels with both classes")
    if not np.isfinite(scores).all():
        raise ValueError("Non-finite scores")
    order = np.argsort(-scores, kind="stable")
    sorted_truth, sorted_scores = truth[order], scores[order]
    last_in_tie = np.r_[sorted_scores[1:] != sorted_scores[:-1], True]
    tp = np.cumsum(sorted_truth, dtype=np.int64)[last_in_tie]
    selected = np.flatnonzero(last_in_tie) + 1
    fp = selected - tp
    precision = tp / selected
    candidates = np.flatnonzero((precision >= target) & (tp >= minimum_tp))
    if len(candidates) == 0:
        return {"attained_on_known_development_labels": False,
                "target_precision": target, "minimum_true_positives": minimum_tp}
    i = candidates[-1]
    # One-sided 95% exact binomial lower bound, conditional on the selected
    # threshold. Because the threshold is chosen here, it is descriptive only.
    lower = float(beta.ppf(0.05, int(tp[i]), int(fp[i]) + 1))
    return {"attained_on_known_development_labels": True,
            "target_precision": target, "minimum_true_positives": minimum_tp,
            "threshold": float(sorted_scores[selected[i] - 1]),
            "tp": int(tp[i]), "fp": int(fp[i]),
            "fn": int(truth.sum() - tp[i]), "tn": int((truth == 0).sum() - fp[i]),
            "empirical_precision": float(precision[i]),
            "known_positive_recall": float(tp[i] / truth.sum()),
            "one_sided_95pct_binomial_lower_bound_descriptive": lower,
            "lower_bound_qualifies_99pct": bool(lower >= target)}


def audit(run_dir: Path) -> dict:
    import numpy as np
    import torch
    from transformers import AutoModel, AutoTokenizer
    from linkerforge.detector.cnn import ThreeLayerDFLCNN

    torch.set_num_threads(min(4, torch.get_num_threads()))

    run_dir = run_dir.resolve()
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest["status"] != "completed":
        raise ValueError("Training run is not complete")
    config_path = ROOT / manifest.get("config_path", str(CONFIG.relative_to(ROOT)))
    config, _, _, dev = preflight(config_path)
    if sha(config_path) != manifest["config_sha256"]:
        raise ValueError("Run config hash mismatch")
    if manifest["split_manifest_sha256"] != config["split_manifest_sha256"]:
        raise ValueError("Run split hash mismatch")
    checkpoint_path = run_dir / "best.pt"
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    if checkpoint["config_sha256"] != manifest["config_sha256"]:
        raise ValueError("Checkpoint config hash mismatch")
    tokenizer = AutoTokenizer.from_pretrained(config["encoder"], local_files_only=True)
    backbone = AutoModel.from_pretrained(config["encoder"], local_files_only=True).eval()
    backbone.requires_grad_(False)
    head = ThreeLayerDFLCNN(backbone.config.hidden_size).eval()
    head.load_state_dict(checkpoint["head"])
    truth, scores = [], []
    unknown_count = 0
    dev_windows = windows(dev, config["max_residues_per_window"])
    for offset in range(0, len(dev_windows), config["batch_size"]):
        group = dev_windows[offset:offset + config["batch_size"]]
        sequences = [item[0] for item in group]
        tokens = tokenizer(sequences, return_tensors="pt", padding=True,
                           add_special_tokens=True, truncation=False)
        with torch.no_grad():
            features = backbone(**tokens).last_hidden_state[:, 1:1 + max(map(len, sequences))]
            probabilities = head(features).sigmoid().cpu().numpy()
        for index, (_, labels, _) in enumerate(group):
            for label, score in zip(labels, probabilities[index, :len(labels)], strict=True):
                if label == "x":
                    unknown_count += 1
                else:
                    truth.append(int(label))
                    scores.append(float(score))
    expected = manifest["best_development"]
    if sum(truth) != expected["positive_residues"] or len(truth) - sum(truth) != expected["negative_residues"]:
        raise ValueError("Development label counts differ from training manifest")
    from sklearn.metrics import average_precision_score
    aps = float(average_precision_score(np.asarray(truth), np.asarray(scores)))
    # CPU recomputation differs slightly from CUDA around nearly tied scores.
    if abs(aps - expected["masked_residue_aps"]) > 1e-4:
        raise ValueError(f"Recomputed APS mismatch: {aps} versus {expected['masked_residue_aps']}")
    return {
        "schema_version": "1.0.0", "claim_status": "known_development_labels_only_not_mining_certified",
        "run_manifest_sha256": sha(run_dir / "manifest.json"),
        "checkpoint_sha256": sha(checkpoint_path),
        "config_sha256": manifest["config_sha256"],
        "split_manifest_sha256": manifest["split_manifest_sha256"],
        "seed": manifest["seed"], "known_positive_residues": sum(truth),
        "known_negative_residues": len(truth) - sum(truth),
        "unknown_residues_in_scored_windows": unknown_count,
        "recomputed_known_development_aps": aps,
        "precision_99_min50tp": precision_frontier(truth, scores),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.run_dir)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
