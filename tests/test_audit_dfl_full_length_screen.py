import unittest

from scripts.audit_dfl_full_length_screen import full_length_hit


class FullLengthScreenTests(unittest.TestCase):
    def test_exact_and_contained_sequence_is_hard_hit(self):
        sequence = "ACDEFGHIKLMNPQRSTVWY" * 3
        self.assertTrue(full_length_hit(sequence, sequence))
        self.assertTrue(full_length_hit(sequence, "MM" + sequence + "KK"))

    def test_short_shared_segment_is_not_full_length_homology(self):
        shared = "ACDEFGHIKLMNPQRSTVWY" * 3
        left = "A" * 90 + shared + "C" * 90
        right = "W" * 90 + shared + "Y" * 90
        self.assertFalse(full_length_hit(left, right))

    def test_long_identical_region_with_short_terminal_change(self):
        shared = "ACDEFGHIKLMNPQRSTVWY" * 5
        left = shared + "AAAAA"
        right = shared + "CCCCC"
        self.assertTrue(full_length_hit(left, right))


if __name__ == "__main__":
    unittest.main()
