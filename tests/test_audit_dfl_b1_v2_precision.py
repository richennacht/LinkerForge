import unittest

from scripts.audit_dfl_b1_v2_precision import precision_frontier


class PrecisionFrontierTests(unittest.TestCase):
    def test_tied_negative_cannot_be_skipped(self):
        labels = [1] * 50 + [0, 1]
        scores = [0.9] * 52
        self.assertFalse(precision_frontier(labels, scores)["attained_on_known_development_labels"])

    def test_minimum_support(self):
        labels = [1] * 49 + [0]
        scores = [0.9] * 49 + [0.1]
        self.assertFalse(precision_frontier(labels, scores)["attained_on_known_development_labels"])

    def test_zero_false_positives(self):
        labels = [1] * 50 + [0]
        scores = [0.9] * 50 + [0.1]
        result = precision_frontier(labels, scores)
        self.assertEqual((result["tp"], result["fp"], result["fn"]), (50, 0, 0))
        self.assertLess(result["one_sided_95pct_binomial_lower_bound_descriptive"], 0.99)


if __name__ == "__main__":
    unittest.main()
