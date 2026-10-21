"""Synthetic geometric invariants and contract tests, not biological validation."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, Residue, Structure, PDBIO
from jsonschema import ValidationError

from scripts.assess_fusion_surface_burial import main as surface_main
from scripts.assess_fusion_surface_burial import measure, read_protein
from scripts.audit_fusion_structure_evidence import main as evidence_main
from scripts.assess_reduce_probe_overlaps import count_explicit_atoms, parse_probe_records, tool_command


def structure_at(points):
    structure = Structure.Structure("test")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    structure.add(model)
    model.add(chain)
    for index, xyz in enumerate(points, 1):
        residue = Residue.Residue((" ", index, " "), "ALA", " ")
        residue.add(Atom.Atom("CA", np.array(xyz, dtype=float), 0, 1, " ", " CA ", index, element="C"))
        chain.add(residue)
    return structure


def full_atom_structure(points):
    structure = Structure.Structure("test")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    structure.add(model)
    model.add(chain)
    serial = 1
    for index, xyz in enumerate(points, 1):
        residue = Residue.Residue((" ", index, " "), "ALA", " ")
        center = np.array(xyz, dtype=float)
        for name, offset, fullname, element in [
            ("N", (-1.2, 0.1, 0.0), " N  ", "N"),
            ("CA", (0.0, 0.0, 0.0), " CA ", "C"),
            ("C", (1.2, -0.1, 0.0), " C  ", "C"),
            ("O", (1.8, -0.8, 0.1), " O  ", "O"),
        ]:
            residue.add(Atom.Atom(name, center + offset, 0, 1, " ", fullname, serial, element=element))
            serial += 1
        chain.add(residue)
    return structure


def write_pdb(path, structure):
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(path))


def test_direct_fusion_has_no_linker_attributed_burial():
    rows = measure(structure_at([(0, 0, 0), (3, 0, 0)]), 1, 0, 240)
    assert all(row["linker_attributed_burial_angstrom2"] == 0 for row in rows)
    assert all(row["total_added_context_burial_angstrom2"] > 0 for row in rows)


def test_far_apart_regions_have_no_burial():
    rows = measure(structure_at([(0, 0, 0), (50, 0, 0), (100, 0, 0)]), 1, 1, 240)
    assert all(row["total_added_context_burial_angstrom2"] == 0 for row in rows)


def test_near_linker_occludes_domain_and_mapping_is_correct():
    rows = measure(structure_at([(0, 0, 0), (3, 0, 0), (100, 0, 0)]), 1, 1, 240)
    assert rows[0]["linker_attributed_burial_angstrom2"] > 0
    assert rows[1]["domain"] == "B"
    assert rows[1]["domain_residue"] == 1 and rows[1]["fusion_residue"] == 3
    assert rows[1]["total_added_context_burial_angstrom2"] == 0


def test_invalid_boundary_rejected():
    with pytest.raises(ValueError, match="boundaries"):
        measure(structure_at([(0, 0, 0), (3, 0, 0)]), 2, 0, 100)


def test_probe_records_deduplicate_reciprocal_pairs_and_map_regions():
    text = "\n".join([
        ":1->1:bo: A   2 ALA  CA  : A   4 GLY  HA2 :3:-0.450:-0.410:0:0:0:0:0:C:C:0:0:0:10:20",
        ":1->1:bo: A   4 GLY  HA2 : A   2 ALA  CA  :4:-0.510:-0.430:0:0:0:0:0:C:C:0:0:0:20:10",
        ":1->1:wo: A   4 GLY  HA3 : A   7 ALA  CB  :2:-0.620:-0.500:0:0:0:0:0:C:C:0:0:0:20:30",
    ])
    records = parse_probe_records(text, a_length=3, linker_length=2, total_length=8)
    assert len(records) == 2
    assert records[0]["minimum_surface_gap_angstrom"] == -0.62
    assert records[0]["region_pair"] == "domain_b--linker"
    assert records[1]["minimum_surface_gap_angstrom"] == -0.51
    assert records[1]["region_pair"] == "domain_a--linker"


def test_explicit_atom_count_ignores_duplicate_altloc_and_nonprotein():
    text = "\n".join([
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C",
        "ATOM      2  CA AALA A   1       0.000   0.000   0.000  1.00 10.00           C",
        "HETATM    3  O   HOH A   2       1.000   0.000   0.000  1.00 10.00           O",
    ])
    assert count_explicit_atoms(text) == 1


def test_tool_command_records_pinned_manifest_and_no_shell():
    command = tool_command(Path("M:/project/pixi.toml"), "reduce", "-version")
    assert command[:5] == [
        "wsl.exe", "--exec", "/home/linkerforge/.pixi/bin/pixi", "run", "--manifest-path",
    ]
    assert command[-2:] == ["reduce", "-version"]


def test_sequence_mismatch_and_missing_backbone_are_rejected(tmp_path):
    pdb = tmp_path / "partial.pdb"
    writer = PDBIO()
    writer.set_structure(structure_at([(0, 0, 0), (3, 0, 0)]))
    writer.save(str(pdb))
    with pytest.raises(ValueError, match="sequence differs"):
        read_protein(pdb, "GG")
    with pytest.raises(ValueError, match="Backbone atoms missing"):
        read_protein(pdb, "AA")


def test_report_generation_and_evidence_linkage_bind_boundaries(tmp_path, monkeypatch):
    points = [(0, 0, 0), (3.8, 0.5, 0.2), (7.6, -0.3, 0.7),
              (11.4, 0.8, -0.2), (15.2, -0.5, 0.6), (19.0, 0.4, -0.6), (22.8, 0, 0.4)]
    fusion_pdb = tmp_path / "fusion.pdb"
    domain_a_pdb = tmp_path / "domain-a.pdb"
    domain_b_pdb = tmp_path / "domain-b.pdb"
    write_pdb(fusion_pdb, full_atom_structure(points))
    write_pdb(domain_a_pdb, full_atom_structure(points[:3]))
    write_pdb(domain_b_pdb, full_atom_structure(points[4:]))
    sequence = "A" * 7
    sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
    pdb_hash = hashlib.sha256(fusion_pdb.read_bytes()).hexdigest()
    comparison = tmp_path / "comparison.json"
    manifest = tmp_path / "manifest.json"
    comparison.write_text(json.dumps({
        "schema_version": "linkerforge.public-multipair-structure-comparison.v1",
        "status": "complete",
        "records": [{
            "domain_pair": "TEST", "linker_id": "linker_1", "linker_length": 1,
            "pdb": str(fusion_pdb), "pdb_sha256": pdb_hash,
        }],
    }), encoding="utf-8")
    manifest.write_text(json.dumps({"records": [{
        "construct_id": "TEST__linker_1", "domain_sources": ["synthetic A", "synthetic B"],
        "domain_a_length": 3, "linker_length": 1,
        "linker_start_1based": 4, "linker_end_1based": 4,
        "sequence": sequence, "sequence_sha256": sequence_hash,
    }]}), encoding="utf-8")
    surface_dir = tmp_path / "surface"
    monkeypatch.setattr(sys, "argv", [
        "assess_fusion_surface_burial.py", "--comparison", str(comparison),
        "--sequence-manifest", str(manifest), "--output-dir", str(surface_dir),
        "--n-points", "120",
    ])
    surface_main()
    surface_json = surface_dir / "surface_burial.json"
    surface = json.loads(surface_json.read_text(encoding="utf-8"))
    assert surface["records"][0]["boundaries"] == {
        "domain_a_length": 3, "linker_length": 1, "domain_b_length": 3,
        "linker_start_1based": 4, "linker_end_1based": 4,
    }
    assert surface["source_ids"][0].startswith("matched-comparison-sha256:")

    score_json = tmp_path / "scores.json"
    score_json.write_text(json.dumps({
        "plddt": [80.0] * 7, "pae": [[5.0] * 7 for _ in range(7)],
        "ptm": 0.5, "max_pae": 5.0,
    }), encoding="utf-8")
    evidence_json = tmp_path / "evidence.json"
    base_args = [
        "audit_fusion_structure_evidence.py", "--fusion-pdb", str(fusion_pdb),
        "--score-json", str(score_json), "--domain-a-reference", str(domain_a_pdb),
        "--domain-b-reference", str(domain_b_pdb), "--domain-a-length", "3",
        "--linker-length", "1", "--surface-burial", str(surface_json),
    ]
    monkeypatch.setattr(sys, "argv", base_args + ["--output", str(evidence_json)])
    evidence_main()
    evidence = json.loads(evidence_json.read_text(encoding="utf-8"))
    assert evidence["linked_evidence"]["surface_burial"]["measurement_count"] == 6
    assert "fixed-coordinate domain surface-burial analysis" not in evidence["not_performed"]

    malformed = json.loads(surface_json.read_text(encoding="utf-8"))
    malformed["records"][0]["domain_surface_buried_by_linker_angstrom2"] = "not-a-number"
    malformed_path = tmp_path / "malformed-surface.json"
    malformed_path.write_text(json.dumps(malformed), encoding="utf-8")
    monkeypatch.setattr(sys, "argv", base_args[:-1] + [str(malformed_path), "--output", str(tmp_path / "malformed-evidence.json")])
    with pytest.raises(ValidationError):
        evidence_main()

    wrong_boundary = json.loads(surface_json.read_text(encoding="utf-8"))
    wrong_boundary["records"][0]["boundaries"]["domain_a_length"] = 2
    wrong_boundary_path = tmp_path / "wrong-boundary.json"
    wrong_boundary_path.write_text(json.dumps(wrong_boundary), encoding="utf-8")
    monkeypatch.setattr(sys, "argv", base_args[:-1] + [str(wrong_boundary_path), "--output", str(tmp_path / "wrong-boundary-evidence.json")])
    with pytest.raises(ValueError, match="boundaries"):
        evidence_main()

    wrong_mapping = json.loads(surface_json.read_text(encoding="utf-8"))
    wrong_mapping["records"][0]["measurements"][0]["domain_residue"] = 2
    wrong_mapping_path = tmp_path / "wrong-mapping.json"
    wrong_mapping_path.write_text(json.dumps(wrong_mapping), encoding="utf-8")
    monkeypatch.setattr(sys, "argv", base_args[:-1] + [str(wrong_mapping_path), "--output", str(tmp_path / "wrong-mapping-evidence.json")])
    with pytest.raises(ValueError, match="identity or domain mapping"):
        evidence_main()
