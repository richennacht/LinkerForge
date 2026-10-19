"""Two-stage domain-gated conformal risk control and abstaining DFL span caller.

Implements checklist items B2 (conformal risk control with abstention under
family-disjoint exchangeability) and B3 (specificity verification on an
adjudicated no-linker / single-domain panel) on top of:
  - Stage 1: `linkerforge.detector.domain_gate.domain_gap_candidates`
    (restricts candidate generation to 5..60 aa gaps between adjacent,
    non-overlapping, distinct structural domains `left.accession != right.accession`,
    deterministically eliminating single-domain interiors, N/C-terminal tails,
    and homo-repeat loops).
  - Stage 2: Frozen Seed-31 `ThreeLayerDFLCNN` residue probabilities aggregated
    over candidate spans, calibrated via Learn-Then-Test (LTT) / Clopper-Pearson
    upper risk bounds with explicit abstention.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Sequence

from linkerforge.detector.domain_gate import CandidateSpan, DomainMatch, domain_gap_candidates
from linkerforge.detector.selective import clopper_pearson_lower


def clopper_pearson_upper_bound(failures: int, trials: int, alpha: float = 0.05) -> float:
    """One-sided (1 - alpha) exact Clopper-Pearson upper bound on failure rate."""
    if trials <= 0:
        return 1.0
    successes = trials - failures
    return float(1.0 - clopper_pearson_lower(successes, trials, alpha=alpha))


@dataclass(frozen=True)
class ConformalCalibrationResult:
    target_fdr: float
    alpha: float
    calibrated_threshold: float | None
    attained: bool
    empirical_fdr: float
    ucb_fdr: float
    empirical_recall: float
    abstention_rate: float
    true_positives: int
    false_positives: int
    total_positives: int
    total_negatives: int
    exchangeability_note: str


def calibrate_conformal_threshold(
    labels: Sequence[int],
    scores: Sequence[float],
    *,
    target_fdr: float = 0.10,
    alpha: float = 0.05,
    min_true_positives: int = 5,
) -> ConformalCalibrationResult:
    """Learn-Then-Test (LTT) finite-sample risk control using exact binomial UCB.

    Finds the lowest threshold `tau` on a family-disjoint calibration set such
    that the one-sided `1 - alpha` Clopper-Pearson upper confidence bound on the
    false discovery rate `FP / (TP + FP)` is at most `target_fdr`, subject to
    `TP >= min_true_positives`. If no threshold satisfies the bound, returns
    `calibrated_threshold=None` and `attained=False` (mandatory abstention).
    """
    if len(labels) != len(scores) or not labels:
        raise ValueError("labels and scores must be non-empty and equal length")
    if not (0.0 < target_fdr < 1.0) or not (0.0 < alpha < 1.0):
        raise ValueError("target_fdr and alpha must be in (0, 1)")
    if min_true_positives < 1:
        raise ValueError("min_true_positives must be >= 1")

    clean_labels: list[int] = []
    clean_scores: list[float] = []
    for y, s in zip(labels, scores):
        if y not in (0, 1):
            raise ValueError(f"Calibration labels must be binary 0/1, got {y!r}")
        val = float(s)
        if not isfinite(val) or not (0.0 <= val <= 1.0):
            raise ValueError(f"Calibration score must be finite in [0, 1], got {val!r}")
        clean_labels.append(int(y))
        clean_scores.append(val)

    total_pos = sum(clean_labels)
    total_neg = len(clean_labels) - total_pos
    if total_pos == 0 or total_neg == 0:
        raise ValueError("Calibration set must contain both positive (1) and negative (0) examples")

    unique_cutoffs = sorted(set(clean_scores))
    best: ConformalCalibrationResult | None = None
    note = (
        "Finite-sample Clopper-Pearson (1 - alpha) upper bound holds under "
        "exchangeability across family-disjoint calibration clusters (B1 V2)."
    )

    for cutoff in unique_cutoffs:
        tp = sum(1 for y, s in zip(clean_labels, clean_scores) if s >= cutoff and y == 1)
        fp = sum(1 for y, s in zip(clean_labels, clean_scores) if s >= cutoff and y == 0)
        called = tp + fp
        if tp < min_true_positives or called == 0:
            continue
        emp_fdr = fp / called
        ucb_fdr = clopper_pearson_upper_bound(fp, called, alpha=alpha)
        if ucb_fdr <= target_fdr:
            recall = tp / total_pos
            abstention = 1.0 - (called / len(clean_labels))
            candidate = ConformalCalibrationResult(
                target_fdr=target_fdr,
                alpha=alpha,
                calibrated_threshold=float(cutoff),
                attained=True,
                empirical_fdr=emp_fdr,
                ucb_fdr=ucb_fdr,
                empirical_recall=recall,
                abstention_rate=abstention,
                true_positives=tp,
                false_positives=fp,
                total_positives=total_pos,
                total_negatives=total_neg,
                exchangeability_note=note,
            )
            if best is None or candidate.empirical_recall > best.empirical_recall:
                best = candidate

    if best is not None:
        return best

    return ConformalCalibrationResult(
        target_fdr=target_fdr,
        alpha=alpha,
        calibrated_threshold=None,
        attained=False,
        empirical_fdr=1.0,
        ucb_fdr=1.0,
        empirical_recall=0.0,
        abstention_rate=1.0,
        true_positives=0,
        false_positives=0,
        total_positives=total_pos,
        total_negatives=total_neg,
        exchangeability_note=note + " Target risk not attained; caller must abstain.",
    )


def call_two_stage_interdomain_spans(
    sequence: str,
    domain_matches: Sequence[DomainMatch],
    residue_scores: Sequence[float],
    *,
    sequence_verified: bool = True,
    high_confidence_threshold: float = 0.7351981997489929,
    provisional_threshold: float = 0.50,
    min_length: int = 5,
    max_length: int = 60,
) -> dict:
    """Apply Stage 1 domain-boundary gate followed by Stage 2 CNN span scoring.

    Single-domain proteins (`< 2` domains) and homo-repeat domain pairs
    (`left_accession == right_accession`) deterministically abstain at Stage 1
    with zero called spans.
    """
    if len(residue_scores) != len(sequence):
        raise ValueError("residue_scores length must match sequence length")
    candidates, abstention_reason = domain_gap_candidates(
        sequence,
        domain_matches,
        sequence_verified=sequence_verified,
        min_length=min_length,
        max_length=max_length,
    )
    if abstention_reason is not None or not candidates:
        return {
            "stage1_passed": False,
            "stage1_abstention_reason": abstention_reason or "no_eligible_gap",
            "candidate_spans": [],
            "called_high_confidence_spans": 0,
            "called_provisional_spans": 0,
        }

    scored_spans = []
    high_conf_count = 0
    provisional_count = 0
    for span in candidates:
        span_scores = [float(x) for x in residue_scores[span.start - 1 : span.end]]
        if not all(isfinite(x) and 0.0 <= x <= 1.0 for x in span_scores):
            raise ValueError("Residue scores must be finite probabilities in [0, 1]")
        mean_score = sum(span_scores) / len(span_scores)
        max_score = max(span_scores)
        frac_above_high = sum(1 for x in span_scores if x >= high_confidence_threshold) / len(span_scores)
        if mean_score >= high_confidence_threshold or (max_score >= high_confidence_threshold and mean_score >= provisional_threshold):
            tier = "mined_high_confidence_dfl_call"
            high_conf_count += 1
        elif mean_score >= provisional_threshold:
            tier = "mined_provisional_dfl_candidate"
            provisional_count += 1
        else:
            tier = "mined_domain_gap_abstained"
        scored_spans.append(
            {
                "start": span.start,
                "end": span.end,
                "length": span.length,
                "sequence": sequence[span.start - 1 : span.end],
                "left_domain_accession": span.left_accession,
                "right_domain_accession": span.right_accession,
                "cnn_mean_score": round(mean_score, 6),
                "cnn_max_score": round(max_score, 6),
                "cnn_fraction_above_high_threshold": round(frac_above_high, 6),
                "calling_tier": tier,
            }
        )

    return {
        "stage1_passed": True,
        "stage1_abstention_reason": None,
        "candidate_spans": scored_spans,
        "called_high_confidence_spans": high_conf_count,
        "called_provisional_spans": provisional_count,
    }
