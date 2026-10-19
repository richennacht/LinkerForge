import random
import unittest
from collections import Counter

from scripts.audit_dfl_firewall_null import shuffled


class FirewallNullTests(unittest.TestCase):
    def test_residue_shuffle_is_deterministic_and_preserves_composition(self):
        sequence = "ACDEFGHIKLMNPQRSTVWY" * 3
        left = shuffled(sequence, random.Random(173), 1)
        right = shuffled(sequence, random.Random(173), 1)
        self.assertEqual(left, right)
        self.assertEqual(Counter(left), Counter(sequence))
        self.assertNotEqual(left, sequence)

    def test_block_shuffle_preserves_short_blocks_and_length(self):
        sequence = "AAAABBBBCCCCDDDDEEEEFFFF"
        result = shuffled(sequence, random.Random(173), 4)
        self.assertEqual(len(result), len(sequence))
        self.assertEqual(Counter(result[index:index + 4] for index in range(0, len(result), 4)),
                         Counter(sequence[index:index + 4] for index in range(0, len(sequence), 4)))


if __name__ == "__main__":
    unittest.main()
