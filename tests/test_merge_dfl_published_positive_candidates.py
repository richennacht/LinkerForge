"""Synthetic tests for deduplication mask semantics and survey triage."""

import unittest

from scripts.merge_dfl_published_positive_candidates import union_positive_masks
from scripts.survey_dfl_literature_sources import classify, key_of


class PublishedCandidateTests(unittest.TestCase):
    def test_union_observed_positives_discards_unverified_zeros(self):
        self.assertEqual(union_positive_masks("ACDE", ["1xxx", "x1xx"]), "11xx")

    def test_source_one_zero_disagreement_is_masked(self):
        self.assertEqual(union_positive_masks("ACDE", ["10xx", "x10x"]), "1xxx")

    def test_invalid_mask_fails(self):
        with self.assertRaises(ValueError):
            union_positive_masks("ACDE", ["1xx"])

    def test_same_title_preprint_and_paper_deduplicate(self):
        a = {"title": "A DFL dataset.", "doi": "10.1/preprint"}
        b = {"title": "A DFL dataset", "doi": "10.1/paper"}
        self.assertEqual(key_of(a), key_of(b))

    def test_apod_neurology_false_positive_is_excluded(self):
        score, lane, _ = classify({"title": "Differences in neuroanatomical sites of apoD elevation", "abstractText": "Schizophrenia and bipolar disorder."})
        self.assertEqual(lane, "other_or_false_positive")


if __name__ == "__main__":
    unittest.main()
