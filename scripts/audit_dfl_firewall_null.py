"""CPU-only null audit of the declared DP94 local-similarity screen.

Reads only candidate training sequences and the sealed DP94 sequence registry.
Shuffled controls preserve amino-acid composition (or short intact blocks),
but destroy long-range homology. This diagnoses specificity of the screening
rule; it never authorizes a weaker benchmark firewall or a training launch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path

if __package__:
    from scripts.preflight_dfl_family_training import RULE, local_hit
else:
    from preflight_dfl_family_training import RULE, local_hit


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shuffled(sequence: str, rng: random.Random, block_size: int) -> str:
    units = ([sequence[index:index + block_size]
              for index in range(0, len(sequence), block_size)]
             if block_size > 1 else list(sequence))
    rng.shuffle(units)
    return "".join(units)


def audit(source_root: Path, *, sample_size: int = 30, seed: int = 173) -> dict:
    if not 2 <= sample_size <= 40:
        raise ValueError("Sample size must be between 2 and 40")
    root = source_root.resolve()
    candidate_manifest_path = root / "data/processed/dfl_positive_candidates_2026_10_12.manifest.json"
    candidate_manifest = json.loads(candidate_manifest_path.read_text(encoding="utf-8"))
    candidate_path = root / "data/processed/dfl_positive_candidates_2026_10_12.jsonl"
    if sha(candidate_path) != candidate_manifest["output_sha256"]:
        raise ValueError("Candidate source hash mismatch")
    candidates = [(row["protein_id"], row["sequence"])
                  for row in map(json.loads, candidate_path.read_text(encoding="utf-8").splitlines())]
    if len(candidates) != candidate_manifest["records"]:
        raise ValueError("Candidate source count mismatch")
    dp_root = root / "data/processed/prospective_dfl_holdout_dp94_v1"
    dp_manifest_path = dp_root / "manifest.json"
    dp_manifest = json.loads(dp_manifest_path.read_text(encoding="utf-8"))
    if dp_manifest["status"] != "candidate_not_certified_independent":
        raise ValueError("Unexpected DP94 registry status")
    dp_path = dp_root / "sequences.fasta"
    raw = dp_path.read_text(encoding="utf-8").replace("\r\n", "\n")
    if hashlib.sha256(raw.encode()).hexdigest() != dp_manifest["sequence_fasta_sha256"]:
        raise ValueError("DP94 sequence hash mismatch")
    lines = raw.splitlines()
    if len(lines) != 2 * dp_manifest["holdout_proteins"]:
        raise ValueError("Unexpected DP94 record count")
    dp94 = [(lines[index][1:].split()[0], lines[index + 1])
            for index in range(0, len(lines), 2)]
    if [key for key, _ in dp94] != dp_manifest["holdout_ids"]:
        raise ValueError("DP94 ID roster mismatch")
    rng = random.Random(seed)
    candidates = rng.sample(candidates, sample_size)
    dp94 = rng.sample(dp94, sample_size)
    counts = {name: 0 for name in ("observed", "residue_shuffle", "block12_shuffle")}
    pairs = sample_size * sample_size
    for _, candidate_sequence in candidates:
        controls = {"observed": candidate_sequence,
                    "residue_shuffle": shuffled(candidate_sequence, rng, 1),
                    "block12_shuffle": shuffled(candidate_sequence, rng, 12)}
        for _, benchmark_sequence in dp94:
            for name, sequence in controls.items():
                counts[name] += bool(local_hit(sequence, benchmark_sequence))
    return {"status": "diagnostic_null_not_training_authorization",
            "source_sha256": sha(candidate_path), "source_manifest_sha256": sha(candidate_manifest_path),
            "dp94_sequence_manifest_sha256": sha(dp_manifest_path),
            "dp94_fasta_sha256_canonical_lf": dp_manifest["sequence_fasta_sha256"],
            "rule": RULE, "seed": seed, "sample_size_each_side": sample_size,
            "pairs_per_condition": pairs, "hits": counts,
            "rates": {name: count / pairs for name, count in counts.items()},
            "limitation": "Composition controls are not family labels; this diagnostic cannot relax the declared screen."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--sample-size", type=int, default=30)
    parser.add_argument("--seed", type=int, default=173)
    args = parser.parse_args()
    print(json.dumps(audit(args.source_root, sample_size=args.sample_size, seed=args.seed), indent=2))
