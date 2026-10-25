"""Part 4 pipeline package: scorecard, learned ranker, non-circular evaluation, real-coordinate demo."""

from .evaluation import (
    GENERATORS,
    LENGTH_POLICIES,
    corruption_check,
    decoy_ranking_benchmark,
    generator_comparison_on_cath,
    geometry_length,
    scorecard_on_real_linkers,
)
from .fusion_demo import run_fusion_demo
from .ranker import (
    BASELINE_NAMES,
    SYNTHETIC_NEGATIVE_TYPES,
    AuthoritySequenceIndex,
    HardRejectionScorecard,
    PairwiseCompatibilityRanker,
    evaluate_hard_rejection_scorecard,
    fit_pairwise_compatibility_ranker,
    generate_synthetic_negatives,
    record_query_context,
)

__all__ = [
    "AuthoritySequenceIndex",
    "BASELINE_NAMES",
    "GENERATORS",
    "HardRejectionScorecard",
    "LENGTH_POLICIES",
    "PairwiseCompatibilityRanker",
    "SYNTHETIC_NEGATIVE_TYPES",
    "corruption_check",
    "decoy_ranking_benchmark",
    "evaluate_hard_rejection_scorecard",
    "fit_pairwise_compatibility_ranker",
    "generate_synthetic_negatives",
    "generator_comparison_on_cath",
    "geometry_length",
    "record_query_context",
    "run_fusion_demo",
    "scorecard_on_real_linkers",
]
