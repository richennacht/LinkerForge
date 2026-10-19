"""Fail-closed tests for domain-gated inter-domain candidate generation."""

import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from linkerforge.detector.domain_gate import CandidateSpan, DomainMatch, domain_gap_candidates, rank_candidate_span
from scripts.audit_b1_domain_gate import parse_fasta, parse_pfam


class DomainGateTests(unittest.TestCase):
    def setUp(self):
        self.sequence = "A" * 30 + "G" * 10 + "C" * 30
        self.domains = [DomainMatch(1, 30, "PF1", "domain"), DomainMatch(41, 70, "PF2", "domain")]

    def test_valid_interdomain_gap_is_only_candidate(self):
        spans, reason = domain_gap_candidates(self.sequence, self.domains, sequence_verified=True)
        self.assertIsNone(reason)
        self.assertEqual(spans, [CandidateSpan(31, 40, "PF1", "PF2")])
        self.assertEqual(rank_candidate_span([0.9] * 70, spans[0])["decision"], "abstain_unvalidated_threshold")

    def test_unverified_or_ambiguous_annotations_abstain(self):
        for matches, verified, expected in [
            (self.domains, False, "sequence_not_verified"),
            ([self.domains[0]], True, "fewer_than_two_domains"),
            ([self.domains[0], DomainMatch(25, 70, "PF2", "domain")], True, "ambiguous_overlapping_domains"),
            ([self.domains[0], DomainMatch(41, 70, "PF2", "family")], True, "non_domain_or_unidentified_match"),
            ([self.domains[0], DomainMatch(41, 70, "PF1", "domain")], True, "no_eligible_gap"),
        ]:
            with self.subTest(expected=expected):
                spans, reason = domain_gap_candidates(self.sequence, matches, sequence_verified=verified)
                self.assertEqual(spans, [])
                self.assertEqual(reason, expected)

    def test_source_parsers_refuse_discontinuous_match(self):
        self.assertEqual(parse_fasta(b">id\nAAAA\nBBBB\n"), "AAAABBBB")
        response = {"next": None, "results": [{
            "metadata": {"source_database": "pfam", "accession": "PF1", "type": "domain"},
            "proteins": [{"accession": "test", "protein_length": 70,
                          "entry_protein_locations": [{"fragments": [
                              {"start": 1, "end": 30, "dc-status": "CONTINUOUS"},
                              {"start": 41, "end": 70, "dc-status": "CONTINUOUS"}]}]}],
        }]}
        with self.assertRaises(ValueError):
            parse_pfam(json.dumps(response).encode(), "TEST", 70)

    def test_score_validation(self):
        span = CandidateSpan(2, 3, "A", "B")
        self.assertEqual(rank_candidate_span([0.2, 0.8, 0.9], span, minimum_score=0.75)["decision"],
                         "candidate_for_review")
        with self.assertRaises(ValueError):
            rank_candidate_span([0.2, float("nan"), 0.9], span)


if __name__ == "__main__":
    unittest.main()
