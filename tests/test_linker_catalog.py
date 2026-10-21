import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_linker_catalog", ROOT / "scripts/build_linker_catalog.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def catalog_sources_available():
    config = json.loads(module.CONFIG.read_text(encoding="utf-8"))
    paths = [config["sources"]["disprot"], config["sources"]["caid3"]]
    paths.extend(config["sources"]["dld"])
    return all((ROOT / path).is_file() for path in paths)


class CatalogTests(unittest.TestCase):
    def test_runs_are_one_based_inclusive(self):
        self.assertEqual(list(module.runs([0, 1, 1, 0, 1], 1)), [(2, 3), (5, 5)])

    @unittest.skipUnless(
        catalog_sources_available(),
        "requires non-versioned protected CAID3, DisProt and DLD source artifacts",
    )
    def test_reproducible_catalog_and_protected_partition(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            manifest = module.build(out=out)
            observations = [json.loads(s) for s in (out / "observations.v1.jsonl").read_text().splitlines()]
            unique = [json.loads(s) for s in (out / "unique_sequences.v1.jsonl").read_text().splitlines()]
            self.assertEqual(manifest["observation_count"], len(observations))
            self.assertEqual(manifest["unique_sequence_count"], len(unique))
            self.assertEqual(sum(len(row["observation_ids"]) for row in unique), len(observations))
            self.assertTrue(all(row["usage_status"] == "protected_benchmark_only"
                                for row in observations if row["source"] == "CAID3_Linkers"))
            self.assertTrue(all(row["usage_status"] == "quarantine_pending_overlap_audit"
                                for row in observations if row["source"] == "DisProt_2026_06"))
            self.assertTrue(all(row["end_1based"] - row["start_1based"] + 1 == len(row["sequence"])
                                for row in observations))
            first_hash = manifest["outputs_sha256"]
            self.assertEqual(first_hash, module.build(out=out)["outputs_sha256"])
            try:
                import jsonschema
            except ImportError:
                return
            for rows, schema_file in [(observations, "linker_catalog_observation.v1.schema.json"),
                                      (unique, "linker_catalog_unique.v1.schema.json")]:
                schema = json.loads((ROOT / "schemas" / schema_file).read_text())
                validator = jsonschema.Draft202012Validator(schema)
                for row in rows:
                    validator.validate(row)


if __name__ == "__main__":
    unittest.main()
