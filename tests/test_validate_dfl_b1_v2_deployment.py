import unittest

from scripts.validate_dfl_b1_v2_deployment import all_windows


class AllWindowTests(unittest.TestCase):
    def test_unknown_only_windows_are_scored(self):
        rows = [{"protein_id": "p", "sequence": "ACDEFG", "labels": "xxx10x"}]
        self.assertEqual(all_windows(rows, 3), [("ACD", "xxx", "p"), ("EFG", "10x", "p")])

    def test_window_order_does_not_depend_on_labels(self):
        rows = [{"protein_id": "p", "sequence": "ACDEFG", "labels": "xxx10x"}]
        changed = [{**rows[0], "labels": "111111"}]
        self.assertEqual([w[0] for w in all_windows(rows, 3)], [w[0] for w in all_windows(changed, 3)])


if __name__ == "__main__":
    unittest.main()
