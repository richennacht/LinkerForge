from __future__ import annotations

from pathlib import Path

from scripts.preflight_caid3_protected_training import failure_report


def test_missing_protected_reference_is_machine_readable_and_fail_closed() -> None:
    result = failure_report(Path("split/manifest.json"), FileNotFoundError("protected FASTA missing"))
    assert result["status"] == "failed_closed"
    assert result["training_authorized"] is False
    assert "protected FASTA missing" in result["reason"]
    assert "never be bypassed" in result["claim_boundary"]
