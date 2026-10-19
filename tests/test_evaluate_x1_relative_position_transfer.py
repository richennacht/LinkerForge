from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.evaluate_x1_relative_position_transfer import (
    REL_POS_3D_FEATURES,
    compute_feature_vector,
    evaluate_x1,
    extract_relative_position_geometry,
)


def _format_pdb_atom(serial: int, atom_name: str, res_name: str, chain_id: str, res_seq: int, x: float, y: float, z: float) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name:<4s} {res_name:>3s} {chain_id:1s}"
        f"{res_seq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
    )


class EvaluateX1RelativePositionTransferTest(unittest.TestCase):
    def test_extract_relative_position_geometry_is_se3_invariant(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rcsb_dir = Path(tmp)
            coords = np.array(
                [[float(i), float((i % 4) * 2.0), float(i * 1.3)] for i in range(1, 31)],
                dtype=np.float64,
            )

            def write_pdb(pdb_id: str, pts: np.ndarray) -> str:
                lines = [
                    _format_pdb_atom(idx + 1, "CA", "ALA", "A", idx + 1, float(pt[0]), float(pt[1]), float(pt[2]))
                    for idx, pt in enumerate(pts)
                ]
                raw = ("".join(lines) + "END\n").encode("utf-8")
                (rcsb_dir / f"{pdb_id.upper()}.pdb").write_bytes(raw)
                return hashlib.sha256(raw).hexdigest()

            sha_base = write_pdb("1AAA", coords)
            rot = np.array(
                [
                    [0.0, -1.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
                dtype=np.float64,
            )
            shift = np.array([25.0, -13.0, 7.5], dtype=np.float64)
            sha_rot = write_pdb("1BBB", coords @ rot.T + shift)

            def make_rec(pdb_id: str, sha: str) -> dict:
                return {
                    "example_id": f"{pdb_id}_A",
                    "source": {
                        "pdb_id": pdb_id,
                        "chain_id": "A",
                        "coordinate_sha256": sha,
                        "cath_author_numbered_range": [1, 30],
                    },
                    "domain_a": {"start": 1, "end": 10},
                    "domain_b": {"start": 16, "end": 30},
                    "connection_features": {
                        "port_a": {"accessibility_proxy": 0.5, "contact_density": 0.4},
                        "port_b": {"accessibility_proxy": 0.6, "contact_density": 0.3},
                    },
                }

            geom_base = extract_relative_position_geometry(make_rec("1AAA", sha_base), rcsb_dir)
            geom_rot = extract_relative_position_geometry(make_rec("1BBB", sha_rot), rcsb_dir)
            for key in geom_base:
                self.assertAlmostEqual(geom_base[key], geom_rot[key], places=4)

            vec_base = compute_feature_vector(geom_base, "GGGGSGGGGS", REL_POS_3D_FEATURES)
            vec_rot = compute_feature_vector(geom_rot, "GGGGSGGGGS", REL_POS_3D_FEATURES)
            np.testing.assert_allclose(vec_base, vec_rot, atol=1e-4)

    def test_evaluate_x1_runs_held_out_group_ablation_on_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rcsb_dir = root / "rcsb"
            rcsb_dir.mkdir()
            records_path = root / "records.jsonl"

            rows = []
            for idx in range(10):
                pdb_id = f"1T{idx:02d}"
                pdb_path = rcsb_dir / f"{pdb_id}.pdb"
                lines = []
                for res_seq in range(1, 31):
                    x = float(res_seq) + 0.2 * idx
                    y = float((res_seq % 5) * (idx + 1))
                    z = float(res_seq * 1.5)
                    lines.append(_format_pdb_atom(res_seq, "CA", "ALA", "A", res_seq, x, y, z))
                raw_bytes = ("".join(lines) + "END\n").encode("utf-8")
                pdb_path.write_bytes(raw_bytes)
                coord_sha = hashlib.sha256(raw_bytes).hexdigest()

                fam_a = f"1.10.{10 + (idx % 5)}.10"
                fam_b = f"2.20.{20 + (idx % 5)}.10"
                rows.append(
                    {
                        "example_id": f"{pdb_id}_A",
                        "homology_cluster_id": f"{fam_a}__{fam_b}",
                        "source": {
                            "pdb_id": pdb_id,
                            "chain_id": "A",
                            "coordinate_sha256": coord_sha,
                            "cath_author_numbered_range": [1, 30],
                        },
                        "domain_a": {"family": fam_a, "start": 1, "end": 10},
                        "domain_b": {"family": fam_b, "start": 16, "end": 30},
                        "linker": {
                            "sequence": "GGGGS" * (1 + (idx % 2)) + ("A" if idx % 3 == 0 else "S"),
                            "length": 5 * (1 + (idx % 2)) + 1,
                        },
                        "connection_features": {
                            "port_a": {"accessibility_proxy": 0.5, "contact_density": 0.4},
                            "port_b": {"accessibility_proxy": 0.6, "contact_density": 0.3},
                        },
                    }
                )
            records_path.write_text(
                "\n".join(json.dumps(r) for r in rows) + "\n",
                encoding="utf-8",
            )
            summary = evaluate_x1(records_path, rcsb_dir, seeds=(7, 17))
            self.assertEqual(
                summary["status"],
                "completed_held_out_domain_pair_ablation",
            )
            self.assertEqual(summary["coordinate_complete_records"], 10)
            self.assertIn("sequence_only", summary["aggregate_metrics"])
            self.assertIn("sequence_plus_relative_position_3d", summary["aggregate_metrics"])


if __name__ == "__main__":
    unittest.main()
