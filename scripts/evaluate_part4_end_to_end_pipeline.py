"""Evaluate Part 4 with non-circular tests and a demonstration on real coordinates.

Writes `reports/part4_end_to_end_ranker_and_sfgfp_mcherry_2026-10-19.json`.
Sections:
  scorecard_on_real_linkers      rule sensitivity on real linkers (geometry rules: CATH tier only)
  corruption_check               constructed corruptions (by construction; NOT evidence of quality)
  decoy_ranking                  true linker vs length/rigidity-matched real linkers (chance 1/(K+1))
  generator_comparison           six design strategies on 91 coordinate-complete natural linkers
  fusion_demo                    schematic assembly on real 2B3P / 2H5Q coordinates
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.dataset_v1 import (  # noqa: E402
    DEFAULT_DATASET_DIR,
    load_dataset_manifest,
    load_linker_records,
)
from linkerforge.pipeline import (  # noqa: E402
    LENGTH_POLICIES,
    AuthoritySequenceIndex,
    corruption_check,
    decoy_ranking_benchmark,
    fit_pairwise_compatibility_ranker,
    generator_comparison_on_cath,
    run_fusion_demo,
    scorecard_on_real_linkers,
)
from linkerforge.properties import load_library  # noqa: E402

MINED_SAMPLE = 600
FINAL_MINED_SAMPLE = 2000
SAMPLE_SEED = 20261019


def run(dataset_dir: Path, out_path: Path) -> dict:
    manifest = load_dataset_manifest(dataset_dir)
    manifest_sha = hashlib.sha256((dataset_dir / "manifest.json").read_bytes()).hexdigest()
    cath = list(load_linker_records(dataset_dir, tiers=["natural_cath_pdb_verified"]))
    disprot = list(load_linker_records(dataset_dir, tiers=["natural_disprot_curated_dfl"]))
    natural = cath + disprot
    mined_hi = list(load_linker_records(dataset_dir, tiers=["mined_high_confidence_dfl_call"]))
    mined_prov = list(load_linker_records(dataset_dir, tiers=["mined_provisional_dfl_candidate"]))
    templates = list(load_linker_records(dataset_dir, tiers=["literature_engineered_template"]))

    rng = random.Random(SAMPLE_SEED)
    shuffled = mined_hi[:]
    rng.shuffle(shuffled)
    mined_eval = shuffled[:MINED_SAMPLE]
    mined_decoys = shuffled[MINED_SAMPLE:]
    mined_final = shuffled[:FINAL_MINED_SAMPLE]

    print("1/5 rule sensitivity on real linkers", flush=True)
    sens = scorecard_on_real_linkers(
        {
            "natural_cath_pdb_verified": cath,
            "natural_disprot_curated_dfl": disprot,
            "literature_engineered_template": templates,
            "mined_high_confidence_dfl_call": mined_hi,
            "mined_provisional_dfl_candidate": mined_prov,
        }
    )

    print("2/5 constructed corruptions (by construction)", flush=True)
    corr = {
        "natural_148": corruption_check(natural),
        "mined_high_confidence_sample_600": corruption_check(mined_eval),
    }

    print("3/5 decoy ranking", flush=True)
    decoy = {
        "natural_148_queries_vs_mined_decoys": decoy_ranking_benchmark(natural, mined_hi),
        "mined_sample_600_queries_vs_unsampled_mined_decoys": decoy_ranking_benchmark(
            mined_eval, mined_decoys
        ),
    }

    print("4/5 generator comparison", flush=True)
    gen = {p: generator_comparison_on_cath(natural, length_policy=p) for p in LENGTH_POLICIES}

    print("5/5 fusion demo on real coordinates", flush=True)
    ranker = fit_pairwise_compatibility_ranker(natural + mined_final + templates)
    index = AuthoritySequenceIndex(natural + mined_hi + templates)
    demo = run_fusion_demo(
        ranker,
        index,
        load_library(),
        pdb_dir=ROOT / "data" / "raw" / "pdb",
        out_dir=ROOT / "reports" / "structures",
    )

    payload = {
        "schema_version": "2.0.0",
        "dataset_id": manifest["dataset_id"],
        "dataset_manifest_sha256": manifest_sha,
        "sampling": {
            "seed": SAMPLE_SEED,
            "mined_high_confidence_available": len(mined_hi),
            "mined_eval_sample": len(mined_eval),
            "mined_final_ranker_sample": len(mined_final),
        },
        "ranker": {
            "feature_names": list(ranker.feature_names),
            "weights": list(ranker.weights),
            "initialisation": "zero",
            "outputs_fusion_success_probability": False,
        },
        "scorecard_on_real_linkers": sens,
        "corruption_check": corr,
        "decoy_ranking": decoy,
        "generator_comparison": gen,
        "fusion_demo": demo,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR)
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "reports" / "part4_end_to_end_ranker_and_sfgfp_mcherry_2026-10-19.json",
    )
    args = ap.parse_args()
    res = run(args.dataset_dir, args.out)
    print("\nRULE SENSITIVITY")
    for tier, v in res["scorecard_on_real_linkers"].items():
        print(f"  {tier:36s} n={v['n_records']:6d} all-rules-pass={v['all_rules_pass_rate']}")
    for k, v in res["decoy_ranking"].items():
        print(f"\nDECOY {k}: n={v['n_queries']} skipped={v['n_skipped_too_few_decoys']} chance={v['chance']}")
        for a, m in v["arms"].items():
            print(f"  {a:28s} hit@1={m['hit_at_1']['mean']} mrr={m['mrr']['mean']} ci={m['mrr']['ci95']}")
    for p, v in res["generator_comparison"].items():
        print(f"\nGENERATORS ({p}) n={v['n_queries']}")
        for g, m in v["methods"].items():
            print(
                f"  {g:28s} closure={m['span_closure_pass']['mean']} rules={m['all_rules_pass']['mean']} "
                f"len<=2={m['length_within_2aa_of_natural']['mean']} cos={m['composition_cosine_to_natural']['mean']}"
            )
    d = res["fusion_demo"]
    print("\nDEMO placement:", d["placement"])


if __name__ == "__main__":
    main()
