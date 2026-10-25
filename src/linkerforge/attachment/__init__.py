"""Part 3 (Component 4): position-controlled attachment-site scorer and geometric site-pair ranker."""

from linkerforge.attachment.scorer import (
    AttachmentSiteCandidate,
    AttachmentSitePairRanking,
    TwoTowerAttachmentScorer,
    diagnose_original_task_leakage,
    evaluate_attachment_scorer_benchmark,
    fit_two_tower_attachment_scorer,
    matched_offset_positions,
    rank_attachment_site_pairs,
    score_domain_attachment_sites,
)

__all__ = [
    "AttachmentSiteCandidate",
    "AttachmentSitePairRanking",
    "TwoTowerAttachmentScorer",
    "diagnose_original_task_leakage",
    "evaluate_attachment_scorer_benchmark",
    "fit_two_tower_attachment_scorer",
    "matched_offset_positions",
    "rank_attachment_site_pairs",
    "score_domain_attachment_sites",
]
