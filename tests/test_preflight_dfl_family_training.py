"""Synthetic-only checks: no benchmark labels or real training records."""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from scripts.preflight_dfl_family_training import RULE, verify


def write(path, value):
    raw = (json.dumps(value, sort_keys=True).encode() if isinstance(value, dict)
           else value.encode())
    path.write_bytes(raw)
    return hashlib.sha256(raw).hexdigest()


def fixture(base):
    registries = []
    for name, key, sequence in (("caid3", "B1", "M" * 80), ("dp94", "B2", "W" * 80)):
        fasta = f">{key}\n{sequence}\n"
        fasta_hash = write(base / f"{name}.fasta", fasta)
        registry = {"schema_version": "1.0.0",
                    "status": "sequence_only_protected" if name == "caid3" else "sequence_only_labels_sealed",
                    "label_seal_state": "retrospective_spent_exclusion_only" if name == "caid3" else "unopened",
                    "source_release": "synthetic-v1",
                    "fasta_path": f"{name}.fasta", "fasta_sha256_canonical_lf": fasta_hash,
                    "ids": [key]}
        registries.append({"name": name, "path": f"{name}.json",
                           "sha256": write(base / f"{name}.json", registry)})
    family = {"schema_version": "1.0.0", "status": "family_clustered",
              "method": "synthetic", "version": "1", "source_sha256": "synthetic-input-hash",
              "membership": {"B1": "b1", "B2": "b2", "TR": "train", "DV": "dev", "AX": "aux"}}
    family_hash = write(base / "family.json", family)
    files = {}
    for role, key, residue in (("train", "TR", "A"),
                               ("development", "DV", "C"),
                               ("auxiliary", "AX", "G")):
        labels = ("1" + "x" * 79) if role != "auxiliary" else "x" * 80
        row = {"protein_id": key, "sequence": residue * 80, "labels": labels,
               "source": "synthetic", "task": "DFL_positive_unlabeled",
               "provenance": {"source_id": key, "input_sha256": "synthetic"}}
        files[role] = {"path": f"{role}.jsonl", "records": 1,
                       "sha256": write(base / f"{role}.jsonl", json.dumps(row) + "\n")}
    split = {"schema_version": "1.0.0", "status": "candidate_family_split_not_certified",
             "local_similarity_rule": RULE, "benchmark_registries": registries,
             "family_manifest": {"path": "family.json", "sha256": family_hash},
             "source_policies": {"synthetic": {"task": "DFL_positive_unlabeled",
                                               "zero_semantics": "none"}}, "files": files}
    write(base / "split.json", split)
    return base / "split.json"


def edit_split(path, change):
    split = json.loads(path.read_text())
    change(split)
    write(path, split)


def edit_row(path, role, change):
    split = json.loads(path.read_text())
    source = path.parent / split["files"][role]["path"]
    row = json.loads(source.read_text())
    change(row)
    split["files"][role]["sha256"] = write(source, json.dumps(row) + "\n")
    write(path, split)


def edit_family(path, change):
    split = json.loads(path.read_text())
    source = path.parent / split["family_manifest"]["path"]
    family = json.loads(source.read_text())
    change(family)
    split["family_manifest"]["sha256"] = write(source, family)
    write(path, split)


class DflFamilyPreflightTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = fixture(Path(self.temp.name))

    def test_honest_registry_states_pass_declared_checks_only(self):
        result = verify(self.path)
        self.assertEqual(result["status"], "pass_declared_sequence_and_family_manifest_only")
        self.assertEqual(result["counts"], dict.fromkeys(("train", "development", "auxiliary"), 1))

    def test_old_split_status_rejected(self):
        edit_split(self.path, lambda s: s.update(status="sequence_screened_not_full_homology_certification"))
        with self.assertRaisesRegex(ValueError, "old split is ineligible"):
            verify(self.path)

    def test_required_registry_cannot_be_omitted(self):
        edit_split(self.path, lambda s: s["benchmark_registries"].pop())
        with self.assertRaisesRegex(ValueError, "registry omitted"):
            verify(self.path)

    def test_changed_registry_bytes_rejected(self):
        (self.path.parent / "dp94.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "Hash mismatch"):
            verify(self.path)

    def test_changed_registry_label_seal_rejected(self):
        def change(split):
            entry = split["benchmark_registries"][0]
            source = self.path.parent / entry["path"]
            registry = json.loads(source.read_text())
            registry["label_seal_state"] = "opened"
            entry["sha256"] = write(source, registry)
        edit_split(self.path, change)
        with self.assertRaisesRegex(ValueError, "Registry seal/status mismatch"):
            verify(self.path)

    def test_caid3_cannot_be_falsely_marked_unopened(self):
        def change(split):
            entry = split["benchmark_registries"][0]
            source = self.path.parent / entry["path"]
            registry = json.loads(source.read_text())
            registry["status"] = "sequence_only_labels_sealed"
            registry["label_seal_state"] = "unopened"
            entry["sha256"] = write(source, registry)
        edit_split(self.path, change)
        with self.assertRaisesRegex(ValueError, "Registry seal/status mismatch"):
            verify(self.path)

    def test_sequence_collision_rejected_even_with_refreshed_hash(self):
        edit_row(self.path, "train", lambda r: r.update(sequence="M" * 80))
        with self.assertRaisesRegex(ValueError, "Benchmark sequence collision"):
            verify(self.path)

    def test_unverified_zero_rejected(self):
        edit_row(self.path, "train", lambda r: r.update(labels="0" * 80))
        with self.assertRaisesRegex(ValueError, "Unverified negative mask"):
            verify(self.path)

    def test_train_without_observed_positives_rejected(self):
        edit_row(self.path, "train", lambda r: r.update(labels="x" * 80))
        with self.assertRaisesRegex(ValueError, "observed DFL positives"):
            verify(self.path)

    def test_family_overlap_rejected(self):
        edit_family(self.path, lambda f: f["membership"].update(DV="train"))
        with self.assertRaisesRegex(ValueError, "Cross-split family overlap"):
            verify(self.path)

    def test_benchmark_family_overlap_rejected(self):
        edit_family(self.path, lambda f: f["membership"].update(TR="b1"))
        with self.assertRaisesRegex(ValueError, "Benchmark-family collision"):
            verify(self.path)

    def test_dld_substitution_rejected(self):
        edit_row(self.path, "train", lambda r: r.update(task="DLD"))
        with self.assertRaisesRegex(ValueError, "source policy"):
            verify(self.path)


if __name__ == "__main__":
    unittest.main()
