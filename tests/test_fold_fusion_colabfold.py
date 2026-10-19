import importlib.util
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "fold_fusion_colabfold.py"
spec = importlib.util.spec_from_file_location("fold_fusion_colabfold", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)


class FullFusionFoldPreflightTests(unittest.TestCase):
    def test_single_fasta_only(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "fusion.fasta"
            fasta.write_text(">fusion\nACDGGGGSEFG\n", encoding="utf-8")
            self.assertEqual(module.read_one_fasta(fasta), "ACDGGGGSEFG")
            fasta.write_text(">one\nACD\n>two\nEFG\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "exactly one"):
                module.read_one_fasta(fasta)

    def test_missing_predictor_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "fusion.fasta"
            output = root / "output"
            fasta.write_text(">fusion\nACDGGGGSEFG\n", encoding="utf-8")
            with self.assertRaises(FileNotFoundError):
                module.run(fasta, output, "not-a-real-colabfold-executable-xyz", 1)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
