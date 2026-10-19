"""Part 2 Flanking-Context and Geometry Retrieval Baseline (`X2`).

Implements Step 4 of `CURRENT_DIRECTION_2026-10-18.md` and `CHECKLIST.md` item `X2`:
evaluates whether template retrieval and nearest-neighbour retrieval from
flanking domain context (and 3D port geometry when available) reconstruct
held-out natural linkers across strictly disjoint `homology_cluster_id` folds.
"""

from .baseline import (
     evaluate_x2_retrieval_benchmark,
    flanking_pair_similarity,
    sequence_nw_identity,
)

__all__ = [
    "evaluate_x2_retrieval_benchmark",
    "flanking_pair_similarity",
    "sequence_nw_identity",
]
