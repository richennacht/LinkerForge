"""Synthetic-only checks of the published nonpositive-record classifier."""

import unittest

from scripts.audit_dfl_published_negative_capacity import is_nonpositive_published


class PublishedNegativeCapacityTests(unittest.TestCase):
    def test_known_zero_but_no_positive_is_candidate_not_gold(self):
        self.assertTrue(is_nonpositive_published({"source": "TransDFL_TR166", "task": "DFL",
                                                   "protein_id": "synthetic", "sequence": "ACDX",
                                                   "labels": "00xx"}))
        self.assertFalse(is_nonpositive_published({"source": "TransDFL_TR166", "task": "DFL",
                                                    "protein_id": "synthetic", "sequence": "ACDX",
                                                    "labels": "1xxx"}))
        self.assertFalse(is_nonpositive_published({"source": "TransDFL_TR166", "task": "DFL",
                                                    "protein_id": "synthetic", "sequence": "ACDX",
                                                    "labels": "xxxx"}))

    def test_wrong_source_or_malformed_masks_fail_closed(self):
        for source, labels in (("other", "0000"), ("TransDFL_TR166", "000?")):
            with self.assertRaises(ValueError):
                is_nonpositive_published({"source": source, "task": "DFL",
                                          "protein_id": "synthetic", "sequence": "ACDX",
                                          "labels": labels})


if __name__ == "__main__":
    unittest.main()
