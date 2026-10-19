import unittest

from scripts.train_dfl_b1_v2_three_layer_cnn import windows


class WindowTests(unittest.TestCase):
    def test_only_observed_windows_enter_loss(self):
        rows = [{"protein_id": "p", "sequence": "ACDEFG", "labels": "xxx10x"}]
        self.assertEqual(windows(rows, 3), [("EFG", "10x", "p")])

    def test_unknown_is_never_converted_to_negative(self):
        rows = [{"protein_id": "p", "sequence": "ACD", "labels": "1xx"}]
        self.assertEqual(windows(rows, 3)[0][1], "1xx")


if __name__ == "__main__":
    unittest.main()
