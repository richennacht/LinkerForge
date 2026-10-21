from __future__ import annotations

from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts.preflight_openmm_production import GIB, build_preflight


class OpenMMProductionPreflightTest(unittest.TestCase):
    def test_projects_observed_storage_and_blocks_without_headroom(self) -> None:
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            smoke_dir = root / "smoke"
            smoke_dir.mkdir()
            started = datetime(2026, 1, 1, tzinfo=timezone.utc)
            replicates = []
            for number in range(1, 4):
                replicate_dir = smoke_dir / f"replicate-{number:03d}"
                replicate_dir.mkdir()
                (replicate_dir / "trajectory.dcd").write_bytes(b"x" * 1000)
                replicates.append(
                    {
                        "replicate": number,
                        "production_started_at": started.isoformat(),
                        "completed_at": (started + timedelta(seconds=10)).isoformat(),
                        "structural_metric_frames": 10,
                    }
                )
            smoke_manifest = root / "smoke.json"
            smoke_manifest.write_text(
                json.dumps(
                    {
                        "status": "complete",
                        "atom_count": 100,
                        "config": {
                            "domain_a_length": 10,
                            "linker_length": 5,
                            "domain_b_length": 10,
                            "production_steps": 100,
                        },
                        "replicates": replicates,
                    }
                ),
                encoding="utf-8",
            )
            config = root / "production.json"
            config.write_text(
                json.dumps(
                    {
                        "domain_a_length": 10,
                        "linker_length": 5,
                        "domain_b_length": 10,
                        "production_steps": 1000,
                        "equilibration_steps": 100,
                        "report_interval_steps": 10,
                        "replicates": 3,
                        "timestep_fs": 2.0,
                    }
                ),
                encoding="utf-8",
            )
            disk_usage = type("DiskUsage", (), {"total": 20 * GIB, "used": 19 * GIB, "free": GIB})()
            with patch("scripts.preflight_openmm_production.shutil.disk_usage", return_value=disk_usage):
                result = build_preflight(config, smoke_manifest, smoke_dir, root)

            self.assertEqual("blocked_insufficient_storage_headroom", result["status"])
            self.assertFalse(result["run_authorized"])
            self.assertEqual(30_000, result["projection"]["trajectory_bytes"])
            self.assertEqual(330.0, result["projection"]["runtime_seconds"])

    def test_independent_trajectory_interval_reduces_only_stored_frames(self) -> None:
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            smoke_dir = root / "smoke"
            smoke_dir.mkdir()
            started = datetime(2026, 1, 1, tzinfo=timezone.utc)
            replicate_dir = smoke_dir / "replicate-001"
            replicate_dir.mkdir()
            (replicate_dir / "trajectory.dcd").write_bytes(b"x" * 1000)
            smoke_manifest = root / "smoke.json"
            smoke_manifest.write_text(json.dumps({
                "status": "complete", "atom_count": 100,
                "config": {"domain_a_length": 10, "linker_length": 5, "domain_b_length": 10, "production_steps": 100},
                "replicates": [{"replicate": 1, "production_started_at": started.isoformat(),
                                "completed_at": (started + timedelta(seconds=10)).isoformat(),
                                "structural_metric_frames": 10, "trajectory_frames": 10}],
            }), encoding="utf-8")
            config = root / "production.json"
            config.write_text(json.dumps({
                "domain_a_length": 10, "linker_length": 5, "domain_b_length": 10,
                "production_steps": 1000, "equilibration_steps": 100,
                "report_interval_steps": 10, "trajectory_report_interval_steps": 50,
                "replicates": 1, "timestep_fs": 2.0,
            }), encoding="utf-8")
            disk_usage = type("DiskUsage", (), {"total": 20 * GIB, "used": 0, "free": 20 * GIB})()
            with patch("scripts.preflight_openmm_production.shutil.disk_usage", return_value=disk_usage):
                result = build_preflight(config, smoke_manifest, smoke_dir, root)
            self.assertEqual(20, result["projection"]["trajectory_frames_total"])
            self.assertEqual(100, result["projection"]["structural_metric_frames_total"])
            self.assertEqual(2_000, result["projection"]["trajectory_bytes"])


if __name__ == "__main__":
    unittest.main()
