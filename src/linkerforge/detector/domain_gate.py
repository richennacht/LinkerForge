"""Conservative domain-boundary candidate gate for an inter-domain DFL task.

Coordinates are 1-based inclusive. Domain matches and their sequence identity
must be verified by the caller; gaps are *candidates*, never DFL labels.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Sequence


@dataclass(frozen=True)
class DomainMatch:
    start: int
    end: int
    accession: str
    kind: str


@dataclass(frozen=True)
class CandidateSpan:
    start: int
    end: int
    left_accession: str
    right_accession: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def domain_gap_candidates(
    sequence: str,
    matches: Sequence[DomainMatch],
    *,
    sequence_verified: bool,
    min_length: int = 5,
    max_length: int = 60,
) -> tuple[list[CandidateSpan], str | None]:
    """Return gaps between adjacent distinct, nonoverlapping domain matches.

    Fails closed for missing sequence verification, duplicate/nested/overlapping
    domain calls, family/repeat/motif-only calls, or invalid coordinates.
    This intentionally sacrifices recall. The reported reason is an abstention,
    not evidence that the protein lacks a DFL.
    """

    if not sequence or any(aa not in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence):
        return [], "invalid_sequence"
    if not sequence_verified:
        return [], "sequence_not_verified"
    if not 1 <= min_length <= max_length:
        raise ValueError("Invalid candidate length limits")
    if len(matches) < 2:
        return [], "fewer_than_two_domains"
    ordered = sorted(matches, key=lambda match: (match.start, match.end, match.accession))
    for match in ordered:
        if match.kind.lower() != "domain" or not match.accession:
            return [], "non_domain_or_unidentified_match"
        if not 1 <= match.start <= match.end <= len(sequence):
            return [], "invalid_domain_coordinates"
    if any(left.end >= right.start for left, right in zip(ordered, ordered[1:])):
        return [], "ambiguous_overlapping_domains"
    candidates = []
    for left, right in zip(ordered, ordered[1:]):
        start, end = left.end + 1, right.start - 1
        if left.accession != right.accession and min_length <= end - start + 1 <= max_length:
            candidates.append(CandidateSpan(start, end, left.accession, right.accession))
    return candidates, None if candidates else "no_eligible_gap"


def rank_candidate_span(
    residue_scores: Sequence[float], candidate: CandidateSpan, *, minimum_score: float | None = None
) -> dict:
    """Score a candidate without promoting it to a confirmed linker.

    A missing threshold always abstains. A supplied threshold is an explicit
    experimental operating point, never an automatic precision certificate.
    """

    if not 1 <= candidate.start <= candidate.end <= len(residue_scores):
        raise ValueError("Candidate outside score vector")
    values = [float(value) for value in residue_scores[candidate.start - 1:candidate.end]]
    if not all(isfinite(value) and 0 <= value <= 1 for value in values):
        raise ValueError("Scores must be finite probabilities")
    mean_score = sum(values) / len(values)
    if minimum_score is not None and not 0 <= minimum_score <= 1:
        raise ValueError("Threshold must be a probability")
    return {
        "start": candidate.start,
        "end": candidate.end,
        "mean_score": mean_score,
        "decision": "abstain_unvalidated_threshold" if minimum_score is None else (
            "candidate_for_review" if mean_score >= minimum_score else "abstain_below_threshold"
        ),
    }
