"""Synthetic annotation semantics and family-graph checks."""

import unittest

from scripts.build_dfl_curated_evidence_pool import family_components, mask_entry


def region(term, start, end, key, *, unpublished=False):
    return {"term_id": term, "start": start, "end": end, "region_id": key,
            "validated": {"curator_id": "synthetic"}, "unpublished": unpublished,
            "reference_id": "synthetic-paper", "reference_source": "pmid",
            "released": "synthetic-release", "ec_id": "synthetic-evidence"}


class CuratedEvidenceTests(unittest.TestCase):
    def test_only_explicit_linker_and_order_regions_are_known(self):
        entry = {"disprot_id": "DPTEST", "acc": "SYNTHETIC", "sequence": "ACDEFGHIKL",
                 "regions": [region("IDPO:0000033", 2, 4, "p"),
                             region("IDPO:0000006", 7, 8, "o"),
                             region("GO:0005515", 5, 6, "other") ]}
        row = mask_entry(entry)
        self.assertEqual(row["labels"], "x111xx00xx")
        self.assertEqual(len(row["evidence"]["positive"]), 1)

    def test_unpublished_and_conflicting_residues_are_not_fabricated(self):
        entry = {"disprot_id": "DPTEST", "sequence": "ACDEFGHIKL",
                 "regions": [region("IDPO:0000033", 2, 4, "p"),
                             region("IDPO:0000006", 4, 7, "o"),
                             region("IDPO:0000033", 8, 9, "hidden", unpublished=True)]}
        row = mask_entry(entry)
        self.assertEqual(row["labels"], "x11x000xxx")
        self.assertEqual(row["conflicting_residues_masked_x"], 1)
        reviewed = mask_entry(entry, include_unpublished_flag=True)
        self.assertEqual(reviewed["labels"], "x11x00011x")
        self.assertTrue(reviewed["evidence"]["positive"][-1]["unpublished_flag"])

    def test_uncited_or_unreleased_regions_remain_unknown(self):
        uncited = region("IDPO:0000033", 2, 3, "uncited")
        uncited["reference_id"] = None
        unreleased = region("IDPO:0000006", 6, 7, "unreleased")
        unreleased["released"] = None
        entry = {"disprot_id": "DPTEST", "sequence": "ACDEFGHIKL",
                 "regions": [uncited, unreleased]}
        self.assertIsNone(mask_entry(entry))
        self.assertIsNone(mask_entry(entry, include_unpublished_flag=True))

    def test_family_edges_are_transitive(self):
        def item(pfam, gene3d):
            return {"families": {"pfam": pfam, "gene3D": gene3d, "uniref50": []}}
        membership = family_components({"a": item(["P1"], []), "b": item(["P1"], ["G1"]),
                                        "c": item([], ["G1"]), "d": item(["P2"], [])})
        self.assertEqual(membership["a"], membership["c"])
        self.assertNotEqual(membership["a"], membership["d"])


if __name__ == "__main__":
    unittest.main()
