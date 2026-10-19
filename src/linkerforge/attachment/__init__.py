"""Part 3 (Component 4): Two-Tower Attachment-Site Scorer and Geometric Site-Pair Ranker."""

from linkerforge.attachment.scorer import (
    AttachmentSiteCandidate,
    AttachmentSitePairRanking,
    TwoTowerAttachmentScorer,
    evaluate_attachment_scorer_benchmark,
    fit_two_tower_attachment_scorer,
    rank_attachment_site_pairs,
    score_domain_attachment_sites,
)

__all__ = [
    "AttachmentSiteCandidate",
    "AttachmentSitePairRanking",
    "TwoTowerAttachmentScorer",
    "evaluate_attachment_scorer_benchmark",
    "fit_two_tower_attachment_scorer",
    "rank_attachment_site_pairs",
    "score_domain_attachment_sites",
]
