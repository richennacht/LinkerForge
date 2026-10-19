"""Bounded, fail-closed three-layer CNN DFL development trial on the B1 V2 split.

This is an exploratory classifier fit, not a high-precision mining certificate.
The CAID3/DP94 labels are never opened. The family-disjoint auxiliary order
panel remains untouched until a separate, frozen calibration protocol exists.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/dfl_b1_v2_three_layer_cnn_2026-10-18.json"
LOCK = ROOT / "artifacts/dfl-b1-v2-cnn-gpu.lock"
sys.path.insert(0, str(ROOT / "src"))

if __package__:
    from scripts.preflight_caid3_protected_training import verify
else:
    from preflight_caid3_protected_training import verify


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preflight(config_path: Path = CONFIG) -> tuple[dict, dict, list[dict], list[dict]]:
    config_path = config_path.resolve()
    if not config_path.is_relative_to(ROOT / "configs"):
        raise ValueError("Trial config must be in the repository configs directory")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    split_path = (ROOT / config["split_manifest"]).resolve()
    if not split_path.is_relative_to(ROOT) or sha(split_path) != config["split_manifest_sha256"]:
        raise ValueError("Pinned B1 V2 split manifest missing or changed")
    result = verify(split_path)
    if result["status"] != "pass_v2_preregistered_sequence_and_family_firewall":
        raise ValueError("B1 V2 benchmark/family preflight did not pass")
    base = split_path.parent
    for role in ("train", "development", "auxiliary"):
        if sha(base / f"{role}.jsonl") != config[f"{role}_sha256"]:
            raise ValueError(f"Pinned B1 V2 {role} data changed")
    train = [json.loads(line) for line in (base / "train.jsonl").read_text(encoding="utf-8").splitlines()]
    dev = [json.loads(line) for line in (base / "development.jsonl").read_text(encoding="utf-8").splitlines()]
    if config["architecture"] != "three_layer_cnn_1_5_1_50_5" or config["encoder"] != "facebook/esm2_t12_35M_UR50D":
        raise ValueError("Unreviewed architecture or encoder")
    if not 0 < config["max_seconds"] <= 1800 or not 1 <= config["epochs"] <= 12:
        raise ValueError("Trial bounds invalid")
    return config, result, train, dev


def windows(rows: list[dict], width: int) -> list[tuple[str, str, str]]:
    return [(row["sequence"][start:start + width], row["labels"][start:start + width], row["protein_id"])
            for row in rows for start in range(0, len(row["sequence"]), width)
            if any(code in "01" for code in row["labels"][start:start + width])]


def run(output_dir: Path, config_path: Path = CONFIG) -> None:
    config, verified, train, dev = preflight(config_path)  # Must finish before GPU allocation.
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite trial output: {output_dir}")
    if LOCK.exists():
        raise RuntimeError(f"Another B1 GPU job may be active: {LOCK}")
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    lock_fd = os.open(LOCK, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    try:
        os.write(lock_fd, json.dumps({"pid": os.getpid(), "output_dir": str(output_dir)}).encode())
    finally:
        os.close(lock_fd)
    try:
        import numpy as np
        import torch
        from sklearn.metrics import average_precision_score, roc_auc_score
        from transformers import AutoModel, AutoTokenizer
        from linkerforge.detector.cnn import ThreeLayerDFLCNN, masked_bce_with_logits

        if not torch.cuda.is_available():
            raise RuntimeError("CUDA unavailable; no CPU fallback for this bounded GPU trial")
        output_dir.mkdir(parents=True)
        torch.manual_seed(config["seed"])
        torch.cuda.manual_seed_all(config["seed"])
        random.seed(config["seed"])
        np.random.seed(config["seed"])
        device = torch.device("cuda")
        started = time.monotonic()
        deadline = started + config["max_seconds"]
        manifest_path = output_dir / "manifest.json"
        manifest = {
            "schema_version": "1.0.0", "status": "initializing", "pid": os.getpid(),
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "config_path": str(config_path.resolve().relative_to(ROOT)),
            "config_sha256": sha(config_path), "split_manifest_sha256": verified["split_manifest_sha256"],
            "train_sha256": config["train_sha256"], "development_sha256": config["development_sha256"],
            "auxiliary_sha256": config["auxiliary_sha256"], "source_sha256": sha(Path(__file__)),
            "pipeline_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "seed": config["seed"], "gpu": torch.cuda.get_device_name(0),
            "model_id": config["encoder"], "claim_status": config["claim_status"],
            "versions": {"python": sys.version.split()[0], "torch": torch.__version__,
                         "transformers": __import__("transformers").__version__},
        }
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        try:
            tokenizer = AutoTokenizer.from_pretrained(config["encoder"], local_files_only=True)
            backbone = AutoModel.from_pretrained(config["encoder"], local_files_only=True).to(device).eval()
            backbone.requires_grad_(False)
            head = ThreeLayerDFLCNN(backbone.config.hidden_size).to(device)
            optimizer = torch.optim.AdamW(head.parameters(), lr=config["learning_rate"])
            train_windows = windows(train, config["max_residues_per_window"])
            dev_windows = windows(dev, config["max_residues_per_window"])
            if not train_windows or not dev_windows:
                raise ValueError("No observed labels in training or development windows")
            manifest["status"] = "running"
            manifest["train_windows"] = len(train_windows)
            manifest["development_windows"] = len(dev_windows)
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

            def embed(group):
                sequences = [item[0] for item in group]
                tokens = tokenizer(sequences, return_tensors="pt", padding=True,
                                   add_special_tokens=True, truncation=False)
                tokens = {name: value.to(device) for name, value in tokens.items()}
                with torch.no_grad():
                    features = backbone(**tokens).last_hidden_state[:, 1:1 + max(map(len, sequences))]
                labels = torch.full(features.shape[:2], -1, dtype=torch.int8, device=device)
                for index, (_, mask, _) in enumerate(group):
                    labels[index, :len(mask)] = torch.tensor(
                        [{"1": 1, "0": 0, "x": -1}[code] for code in mask],
                        dtype=torch.int8, device=device)
                return features, labels

            best_aps = -1.0
            for epoch in range(config["epochs"]):
                if time.monotonic() >= deadline:
                    manifest["status"] = "time_limit"
                    break
                order = list(range(len(train_windows)))
                random.Random(config["seed"] + epoch).shuffle(order)
                head.train()
                losses = []
                for offset in range(0, len(order), config["batch_size"]):
                    if time.monotonic() >= deadline:
                        manifest["status"] = "time_limit"
                        break
                    group = [train_windows[index] for index in order[offset:offset + config["batch_size"]]]
                    features, labels = embed(group)
                    loss = masked_bce_with_logits(head(features), labels)
                    optimizer.zero_grad(set_to_none=True)
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(head.parameters(), 1.0)
                    optimizer.step()
                    losses.append(float(loss.item()))
                if manifest["status"] == "time_limit":
                    break
                head.eval()
                truth, scores = [], []
                with torch.no_grad():
                    for offset in range(0, len(dev_windows), config["batch_size"]):
                        if time.monotonic() >= deadline:
                            manifest["status"] = "time_limit"
                            break
                        group = dev_windows[offset:offset + config["batch_size"]]
                        features, labels = embed(group)
                        probabilities = head(features).sigmoid()
                        known = labels >= 0
                        truth.extend(labels[known].cpu().numpy().tolist())
                        scores.extend(probabilities[known].cpu().numpy().tolist())
                if manifest["status"] == "time_limit":
                    break
                if set(truth) != {0, 1}:
                    raise ValueError("Development scored residues lost a label class")
                truth_array = np.asarray(truth, dtype=np.uint8)
                score_array = np.asarray(scores, dtype=np.float32)
                aps = float(average_precision_score(truth_array, score_array))
                auroc = float(roc_auc_score(truth_array, score_array))
                fmax = 0.0
                for threshold in np.linspace(0, 1, 101):
                    predicted = score_array >= threshold
                    tp = int(np.sum(predicted & (truth_array == 1)))
                    fp = int(np.sum(predicted & (truth_array == 0)))
                    fn = int(np.sum(~predicted & (truth_array == 1)))
                    value = 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0
                    fmax = max(fmax, value)
                result = {"epoch": epoch + 1, "train_loss": sum(losses) / len(losses),
                          "development": {"masked_residue_aps": aps, "masked_residue_auroc": auroc,
                                          "fmax_grid_101": fmax, "positive_residues": int(truth_array.sum()),
                                          "negative_residues": int(len(truth_array) - truth_array.sum())}}
                with (output_dir / "epochs.jsonl").open("a", encoding="utf-8") as stream:
                    stream.write(json.dumps(result) + "\n")
                if aps > best_aps:
                    best_aps = aps
                    checkpoint = output_dir / "best.pt"
                    temporary = output_dir / "best.pt.tmp"
                    torch.save({"head": head.state_dict(), "epoch": epoch + 1,
                                "development": result["development"],
                                "config_sha256": manifest["config_sha256"]}, temporary)
                    temporary.replace(checkpoint)
                    manifest["best_epoch"] = epoch + 1
                    manifest["best_development"] = result["development"]
                manifest["completed_epochs"] = epoch + 1
                manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
                print(json.dumps(result), flush=True)
            if manifest["status"] == "running":
                manifest["status"] = "completed"
        except Exception as error:
            manifest["status"] = "failed"
            manifest["error"] = f"{type(error).__name__}: {error}"
            raise
        finally:
            manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
            manifest["peak_pytorch_vram_bytes"] = torch.cuda.max_memory_allocated()
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    finally:
        LOCK.unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--config", type=Path, default=CONFIG)
    args = parser.parse_args()
    if args.preflight_only:
        config, verified, train, dev = preflight(args.config)
        print(json.dumps({"status": "preflight_passed", "config_sha256": sha(args.config),
                          "split_manifest_sha256": verified["split_manifest_sha256"],
                          "train_proteins": len(train), "development_proteins": len(dev),
                          "max_seconds": config["max_seconds"]}, indent=2, sort_keys=True))
    else:
        run(args.output_dir, args.config)
