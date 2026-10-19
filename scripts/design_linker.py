"""Command-line Part 2 Linker Design Engine (Geometry + Templates + Mined/Natural Authorities).

Answers the three LinkerForge design questions (Q1 feasibility, Q2 attachment/flanking
context, Q3 candidate linker sequences) under explicit rigidity, cleavability,
and topology (`end_to_end_tandem` vs `internal_insertion_two_tether`) constraints.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from linkerforge.dataset_v1 import DEFAULT_DATASET_DIR  # noqa: E402
from linkerforge.design import DesignRequest, design_linker_candidates  # noqa: E402
from linkerforge.geometry.structure import read_chain  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance-angstrom", type=float, required=True)
    parser.add_argument("--rigidity", choices=("flexible", "rigid"), default="flexible")
    parser.add_argument(
        "--cleavability", choices=("non_cleavable", "cleavable"), default="non_cleavable"
    )
    parser.add_argument("--protease", type=str, default=None)
    parser.add_argument(
        "--topology",
        choices=("end_to_end_tandem", "internal_insertion_two_tether"),
        default="end_to_end_tandem",
    )
    parser.add_argument("--requested-residues", type=int, default=None)
    parser.add_argument("--return-distance-angstrom", type=float, default=None)
    parser.add_argument("--left-flank-sequence", type=str, default="")
    parser.add_argument("--right-flank-sequence", type=str, default="")
    parser.add_argument("--left-domain-id", type=str, default=None)
    parser.add_argument("--right-domain-id", type=str, default=None)
    parser.add_argument("--structure-a", type=Path, default=None)
    parser.add_argument("--chain-a", type=str, default=None)
    parser.add_argument("--anchor-a", type=int, default=None)
    parser.add_argument("--structure-b", type=Path, default=None)
    parser.add_argument("--chain-b", type=str, default=None)
    parser.add_argument("--anchor-b", type=int, default=None)
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR)
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    domain_a = (
        read_chain(args.structure_a, args.chain_a)
        if args.structure_a and args.chain_a
        else None
    )
    struct_b = args.structure_b or args.structure_a
    domain_b = (
        read_chain(struct_b, args.chain_b)
        if struct_b and args.chain_b
        else None
    )

    req = DesignRequest(
        distance_angstrom=args.distance_angstrom,
        rigidity=args.rigidity,
        cleavability=args.cleavability,
        protease=args.protease,
        topology=args.topology,
        requested_residues=args.requested_residues,
        return_distance_angstrom=args.return_distance_angstrom,
        left_flank_sequence=args.left_flank_sequence,
        right_flank_sequence=args.right_flank_sequence,
        left_domain_id=args.left_domain_id,
        right_domain_id=args.right_domain_id,
    )
    report = design_linker_candidates(
        req,
        dataset_dir=args.dataset_dir,
        domain_a=domain_a,
        anchor_a=args.anchor_a,
        domain_b=domain_b,
        anchor_b=args.anchor_b,
        top_k=args.top_k,
    )
    text = json.dumps(report.to_dict(), indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)


if __name__ == "__main__":
    main()
