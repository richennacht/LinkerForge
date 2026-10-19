"""Score one FASTA with a previously frozen DFL checkpoint.

This is inference, not training, benchmark evaluation, or a certified linker
caller. Without an explicit threshold it abstains from calling spans.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.detector.regions import AMINO_ACIDS, intervals  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: Path) -> tuple[str, str]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()
             if line.strip()]
    if len(lines) < 2 or not lines[0].startswith(">") or any(
            line.startswith(">") for line in lines[1:]):
        raise ValueError("Expected exactly one nonempty FASTA entry")
    identifier = lines[0][1:].split()[0]
    sequence = "".join(lines[1:]).upper()
    if not identifier or not sequence or set(sequence) - AMINO_ACIDS:
        raise ValueError("FASTA needs an identifier and canonical amino acids")
    return identifier, sequence


def verified_checkpoint(selection_path: Path, artifact_root: Path) -> tuple[Path, dict, dict]:
    """Check the immutable selection and both file hashes before loading pickle."""
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    if (selection.get("frozen_before_caid3") is not True or
            selection.get("selection_basis") != "filtered_development_only" or
            selection.get("arch") != "samlike" or selection.get("variant") is not None):
        raise ValueError("Unsupported or non-frozen checkpoint selection")
    relative = Path(selection["checkpoint"].replace("\\", "/"))
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Checkpoint path must remain inside artifact root")
    root = artifact_root.resolve()
    checkpoint = (root / relative).resolve()
    if not checkpoint.is_relative_to(root):
        raise ValueError("Checkpoint escapes artifact root")
    manifest_path = checkpoint.parent / "manifest.json"
    if sha256(checkpoint) != selection["checkpoint_sha256"]:
        raise ValueError("Frozen checkpoint SHA-256 mismatch")
    if sha256(manifest_path) != selection["checkpoint_manifest_sha256"]:
        raise ValueError("Frozen checkpoint manifest SHA-256 mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "completed" or manifest.get("arch") != "samlike" or
            manifest.get("blacklist_status") != "passed_before_training" or
            manifest.get("model_id") != "facebook/esm2_t12_35M_UR50D"):
        raise ValueError("Checkpoint manifest is not the expected frozen DFL run")
    return checkpoint, selection, manifest


def score_sequence(sequence: str, checkpoint_path: Path, model_id: str,
                   revision: str, device: str) -> list[float]:
    """Use non-overlapping 510-residue windows, as in the frozen evaluator."""
    import torch
    from transformers import AutoModel, AutoTokenizer
    from linkerforge.models.promptable_linker_segmenter import PromptableLinkerSegmenter

    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA unavailable; select --device cpu")
    state = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    if state.get("arch") != "samlike" or state.get("condition") != "disorder_transfer":
        raise ValueError("Checkpoint architecture or condition mismatch")
    tokenizer = AutoTokenizer.from_pretrained(model_id, revision=revision,
                                              local_files_only=True)
    backbone = AutoModel.from_pretrained(model_id, revision=revision,
                                         local_files_only=True).to(device).eval()
    head = PromptableLinkerSegmenter(backbone.config.hidden_size).to(device).eval()
    head.load_state_dict(state["head"], strict=True)
    scores: list[float] = []
    with torch.inference_mode():
        for offset in range(0, len(sequence), 510):
            fragment = sequence[offset:offset + 510]
            tokens = tokenizer(fragment, return_tensors="pt", add_special_tokens=True,
                               truncation=False)
            hidden = backbone(**{key: value.to(device) for key, value in tokens.items()}
                              ).last_hidden_state[:, 1:1 + len(fragment)]
            valid = torch.ones((1, len(fragment)), dtype=torch.bool, device=device)
            concept = torch.zeros(1, dtype=torch.long, device=device)
            logits = head(hidden, valid, concept)["semantic_logits"]
            values = logits.sigmoid()[0].cpu().tolist()
            if len(values) != len(fragment):
                raise ValueError("Encoder length differs from sequence length")
            scores.extend(values)
    return scores


def make_report(identifier: str, sequence: str, scores: list[float], *,
                checkpoint_sha256: str, model_id: str, revision: str,
                threshold: float | None) -> dict:
    if len(scores) != len(sequence) or any(
            not math.isfinite(score) or not 0 <= score <= 1 for score in scores):
        raise ValueError("Invalid residue scores")
    candidates = [] if threshold is None else intervals(sequence, scores, threshold)
    return {
        "schema_version": "1.0.0", "task": "dfl_sequence_inference",
        "sequence_id": identifier, "sequence_length": len(sequence),
        "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        "coordinate_system": "1-based inclusive",
        "residue_scores": scores,
        "abstained_from_span_calls": threshold is None,
        "threshold": threshold, "candidate_spans": candidates,
        "model": {"checkpoint_sha256": checkpoint_sha256, "encoder": model_id,
                  "encoder_revision_used_for_inference": revision},
        "claim_status": "provisional_inference_not_certified_for_mining",
        "limitations": [
            "The original training manifest did not pin the encoder revision; this is the inference revision, not proof of training equivalence.",
            "The selected development threshold, if supplied, is not an independent precision guarantee.",
            "Predicted spans are not verified linker annotations or evidence of domain boundaries.",
            "The frozen checkpoint predates the prospective family-overlap gate and is not a new independent result.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--selection", required=True, type=Path)
    parser.add_argument("--artifact-root", type=Path, default=ROOT)
    parser.add_argument("--encoder-revision", required=True,
                        help="Exact cached encoder revision used for inference")
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--threshold", type=float,
                        help="Optional provisional span threshold; omitted means abstain")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    identifier, sequence = read_fasta(args.fasta)
    checkpoint, _, manifest = verified_checkpoint(args.selection, args.artifact_root)
    if args.threshold is not None and (not math.isfinite(args.threshold) or
                                       not 0 <= args.threshold <= 1):
        parser.error("--threshold must be finite and in [0,1]")
    scores = score_sequence(sequence, checkpoint, manifest["model_id"],
                            args.encoder_revision, args.device)
    report = make_report(identifier, sequence, scores,
                         checkpoint_sha256=sha256(checkpoint),
                         model_id=manifest["model_id"],
                         revision=args.encoder_revision, threshold=args.threshold)
    import torch
    import transformers
    config = {"window_residues": 510, "threshold": args.threshold,
              "encoder_revision": args.encoder_revision, "device": args.device}
    report["provenance"] = {
        "input_fasta_sha256": sha256(args.fasta),
        "selection_sha256": sha256(args.selection),
        "checkpoint_manifest_sha256": sha256(checkpoint.parent / "manifest.json"),
        "inference_config_sha256": hashlib.sha256(
            json.dumps(config, sort_keys=True).encode("utf-8")).hexdigest(),
        "source_sha256": sha256(Path(__file__)),
        "pipeline_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "software": {"python": sys.version.split()[0], "torch": torch.__version__,
                     "transformers": transformers.__version__},
        "evidence_quality": "uncalibrated_model_prediction_not_observation",
    }
    payload = json.dumps(report, indent=2) + "\n"
    if args.output:
        if args.output.exists():
            raise FileExistsError("Refusing to overwrite an inference output")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
