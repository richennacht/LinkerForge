import unittest

from scripts.prepare_uniprot_explicit_flexible_candidates import make_record


BASE = {"catalog_sha256": "a" * 64, "source_url": "synthetic",
        "retrieved_at_utc": "2026-10-18T00:00:00+00:00", "pipeline_commit": "test",
        "script_sha256": "b" * 64, "schema_sha256": "c" * 64,
        "software": {"python": "synthetic"}}


class ExplicitFlexibleSourceTests(unittest.TestCase):
    def test_only_region_is_positive_and_remainder_unknown(self):
        candidate = {"accession": "TEST", "peptide": "GGGS", "start": 3, "end": 6,
                     "note": "Flexible linker", "evidence": "ECO:0000269|PubMed:1"}
        record = make_record(candidate, "AAGGGSAA", BASE)
        self.assertEqual(record["labels"], "xx1111xx")
        self.assertEqual(record["annotation"]["evidence_tier"],
                         "pubmed_supported_region_not_DFL_adjudicated")
        self.assertEqual(record["training_eligibility"], "candidate_not_DFL_gold_or_negative_adjudicated")

    def test_release_sequence_change_is_quarantinable_error(self):
        candidate = {"accession": "TEST", "peptide": "GGGS", "start": 3, "end": 6,
                     "note": "Flexible linker", "evidence": "ECO:0000256|Rule"}
        with self.assertRaisesRegex(ValueError, "does not match"):
            make_record(candidate, "AAGGASAA", BASE)


if __name__ == "__main__":
    unittest.main()
