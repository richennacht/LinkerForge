"""Guard the source-only family projection used by the B1 audit."""

import unittest

from scripts.audit_dfl_disprot_family_metadata import project_entry


class FamilyMetadataProjectionTests(unittest.TestCase):
    def test_projects_only_sequence_and_family_metadata(self):
        entry = {
            "disprot_id": "DP_TEST", "acc": "P_TEST", "sequence": "ACDEFG",
            "features": {
                "pfam": [{"id": "PF00001", "start": 1, "end": 3}],
                "gene3D": [{"id": "G3D_TEST", "start": 2, "end": 6}],
            },
            "uniref50": "UniRef50_P_TEST",
            "regions": [{"start": 1, "end": 2, "annotation": "sealed"}],
            "disprot_consensus": "sealed",
        }
        result = project_entry(entry)
        self.assertEqual(result["families"]["pfam"], ["PF00001"])
        self.assertEqual(result["families"]["gene3D"], ["G3D_TEST"])
        self.assertEqual(result["families"]["uniref50"], ["UniRef50_P_TEST"])
        self.assertNotIn("regions", result)
        self.assertNotIn("disprot_consensus", result)

    def test_rejects_invalid_family_coordinate(self):
        entry = {"disprot_id": "DP_TEST", "sequence": "ACDE",
                 "features": {"pfam": [{"id": "PF00001", "start": 1, "end": 5}]}}
        with self.assertRaises(ValueError):
            project_entry(entry)

    def test_rejects_malformed_feature(self):
        entry = {"disprot_id": "DP_TEST", "sequence": "ACDE",
                 "features": {"pfam": ["PF00001"]}}
        with self.assertRaises(ValueError):
            project_entry(entry)


if __name__ == "__main__":
    unittest.main()
