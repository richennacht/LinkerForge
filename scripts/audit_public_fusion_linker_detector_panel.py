#!/usr/bin/env python3
"""Build and score a public-sequence linker-detector sensitivity panel.

This is not a biological-success benchmark. It asks whether a frozen linker
detector marks the known, deliberately inserted junction span and how often it
also marks residues inside the two public reference domains.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess


DOMAIN_PAIRS = {
    "1UBQ_2GB1": {
        "domain_a": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "domain_b": "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE",
        "sources": ["PDB 1UBQ chain A", "PDB 2GB1 chain A"],
    },
    "2B3P_2H5Q": {
        "domain_a": (
            "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKR"
            "HDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNG"
            "IKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        ),
        "domain_b": (
            "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKH"
            "PADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDG"
            "ALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK"
        ),
        "sources": ["PDB 2B3P sfGFP sequence used by LinkerForge", "PDB 2H5Q mCherry sequence used by LinkerForge"],
    },
    "1CRN_4PTI": {
        "domain_a": "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN",
        "domain_b": "RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCGGA",
        "sources": ["PDB 1CRN chain A (RCSB FASTA)", "PDB 4PTI chain A (RCSB FASTA)"],
    },
    "1PGA_1LMB": {
        "domain_a": "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE",
        "domain_b": "STKKKPLTQEQLEDARRLKAIYEKKKNELGLSQESVADKMGMGQSGVGALFNGINALNAYNAALLAKILKVSVEEFSPSIAREIYEMYEAVS",
        "sources": ["PDB 1PGA chain A (RCSB FASTA)", "PDB 1LMB protein chain C/D (RCSB FASTA)"],
    },
}

LINKERS = {
    "direct_0": {"sequence": "", "class": "direct_control", "source": "LinkerForge required baseline"},
    "ggggs_5": {"sequence": "GGGGS", "class": "flexible_template", "source": "Chen et al. 2013 Table 2"},
    "ggggs_10": {"sequence": "GGGGSGGGGS", "class": "flexible_template", "source": "Chen et al. 2013 Table 2"},
    "ggggs_15": {"sequence": "GGGGSGGGGSGGGGS", "class": "fixed_control", "source": "Chen et al. 2013 Tables 2-3"},
    "waldo_12": {"sequence": "GSAGSAAGSGEF", "class": "flexible_template", "source": "Chen et al. 2013 Table 3, attributing Waldo et al."},
    "bird_18": {"sequence": "KESGSVSSEQLAQFRSLD", "class": "flexible_template", "source": "Chen et al. 2013 Table 3, attributing Bird et al."},
    "eaaak_15": {"sequence": "EAAAKEAAAKEAAAK", "class": "rigid_control", "source": "Chen et al. 2013 Table 3"},
}


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def build_records() -> list[dict]:
    records = []
    for pair_id, pair in DOMAIN_PAIRS.items():
        for linker_id, linker in LINKERS.items():
            sequence = pair["domain_a"] + linker["sequence"] + pair["domain_b"]
            records.append(
                {
                    "construct_id": f"{pair_id}__{linker_id}",
                    "domain_pair": pair_id,
                    "domain_sources": pair["sources"],
                    "domain_a_length": len(pair["domain_a"]),
                    "linker_id": linker_id,
                    "linker_sequence": linker["sequence"],
                    "linker_length": len(linker["sequence"]),
                    "linker_start_1based": len(pair["domain_a"]) + 1 if linker["sequence"] else None,
                    "linker_end_1based": len(pair["domain_a"]) + len(linker["sequence"]) if linker["sequence"] else None,
                    "linker_class": linker["class"],
                    "linker_source": linker["source"],
                    "sequence": sequence,
                    "sequence_sha256": sha256_text(sequence),
                }
            )
    return records


def parse_dflpred(path: Path) -> dict[str, list[float]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) % 3:
        raise ValueError("DFLpred output must contain three non-empty lines per sequence")
    parsed = {}
    for index in range(0, len(lines), 3):
        identifier = lines[index].removeprefix(">")
        sequence = lines[index + 1]
        scores = [float(value) for value in lines[index + 2].split(",")]
        if len(sequence) != len(scores):
            raise ValueError(f"DFLpred length mismatch for {identifier}")
        parsed[identifier] = scores
    return parsed


def summarize(records: list[dict], scores: dict[str, list[float]], threshold: float) -> list[dict]:
    output = []
    for record in records:
        values = scores[record["construct_id"]]
        start = record["domain_a_length"]
        end = start + record["linker_length"]
        linker_values = values[start:end]
        domain_values = values[:start] + values[end:]
        output.append(
            {
                "construct_id": record["construct_id"],
                "domain_pair": record["domain_pair"],
                "linker_id": record["linker_id"],
                "linker_class": record["linker_class"],
                "linker_length": record["linker_length"],
                "linker_mean_score": round(sum(linker_values) / len(linker_values), 6) if linker_values else None,
                "linker_max_score": round(max(linker_values), 6) if linker_values else None,
                "linker_residues_called": sum(value >= threshold for value in linker_values),
                "linker_recall": round(sum(value >= threshold for value in linker_values) / len(linker_values), 6) if linker_values else None,
                "domain_residues_called": sum(value >= threshold for value in domain_values),
                "domain_call_fraction": round(sum(value >= threshold for value in domain_values) / len(domain_values), 6),
                "threshold": threshold,
            }
        )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dflpred-output", type=Path)
    parser.add_argument("--dflpred-revision", default=None)
    parser.add_argument("--threshold", type=float, default=0.18)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    commit_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    worktree_modified = bool(subprocess.run(
        ["git", "status", "--porcelain"], cwd=root, check=True,
        capture_output=True, text=True,
    ).stdout.strip())
    records = build_records()
    fasta = "".join(f">{record['construct_id']}\n{record['sequence']}\n" for record in records)
    (args.output_dir / "panel.fasta").write_text(fasta, encoding="utf-8")
    panel_config_sha = sha256_text(json.dumps(
        {"domain_pairs": DOMAIN_PAIRS, "linkers": LINKERS},
        sort_keys=True,
        separators=(",", ":"),
    ))
    manifest = {
        "schema_version": "1.0.0",
        "record_type": "public_fusion_linker_detector_sensitivity_panel",
        "evidence_quality": "public_sequence_constructed_sensitivity_panel",
        "source_ids": sorted({source for pair in DOMAIN_PAIRS.values() for source in pair["sources"]}),
        "provenance": {
            "commit_sha": commit_sha,
            "worktree_modified": worktree_modified,
            "config_sha256": panel_config_sha,
            "versions": {"python": platform.python_version(), "panel_builder": "1.0.0"},
        },
        "purpose": "Detector recognition and domain-overcall audit; not biological fusion validation",
        "records": records,
        "limitations": [
            "Inserted template spans are known by construction, not experimentally annotated DFL ground truth in these domain contexts.",
            "Detector scores measure DFL propensity, not compatibility, expression, function, or structural success.",
            "The four contexts are a sensitivity panel and cannot estimate population-level accuracy.",
        ],
    }
    (args.output_dir / "panel_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if args.dflpred_output:
        detector_config_sha = sha256_text(json.dumps({
            "detector": "DFLpred",
            "detector_revision": args.dflpred_revision,
            "threshold": args.threshold,
            "panel_config_sha256": panel_config_sha,
        }, sort_keys=True, separators=(",", ":")))
        result = {
            "schema_version": "1.0.0",
            "detector": "DFLpred",
            "detector_revision": args.dflpred_revision,
            "threshold": args.threshold,
            "input_fasta_sha256": sha256_text(fasta),
            "detector_output_sha256": hashlib.sha256(args.dflpred_output.read_bytes()).hexdigest(),
            "evidence_quality": "orthogonal_public_detector_sensitivity_not_compatibility_validation",
            "source_ids": [
                f"panel-fasta-sha256:{sha256_text(fasta)}",
                f"dflpred-revision:{args.dflpred_revision}",
                f"dflpred-output-sha256:{hashlib.sha256(args.dflpred_output.read_bytes()).hexdigest()}",
            ],
            "provenance": {
                "commit_sha": commit_sha,
                "worktree_modified": worktree_modified,
                "config_sha256": detector_config_sha,
                "versions": {"python": platform.python_version(), "audit_script": "1.0.0"},
            },
            "results": summarize(records, parse_dflpred(args.dflpred_output), args.threshold),
            "claim_boundary": "Sensitivity audit on deliberately inserted template spans; not an independent benchmark.",
        }
        (args.output_dir / "dflpred_summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
