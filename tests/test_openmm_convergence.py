import pytest

from scripts.analyze_openmm_convergence import summarize_rows


def rows(values):
    return [{"step": str((index + 1) * 10), "metric": str(value)} for index, value in enumerate(values)]


def test_summarize_rows_reports_block_and_half_drift():
    summary = summarize_rows(rows([1, 1, 2, 2, 3, 3, 4, 4, 5, 5]), blocks=5)
    metric = summary["metrics"]["metric"]
    assert summary["frame_count"] == 10
    assert metric["block_means"] == [1.0, 2.0, 3.0, 4.0, 5.0]
    assert metric["first_half_mean"] == 1.8
    assert metric["second_half_mean"] == 4.2
    assert metric["half_mean_absolute_drift"] == 2.4


def test_summarize_rows_rejects_nondivisible_or_duplicate_steps():
    with pytest.raises(ValueError, match="not divisible"):
        summarize_rows(rows([1, 2, 3]), blocks=2)
    duplicate = rows([1, 2, 3, 4])
    duplicate[1]["step"] = duplicate[0]["step"]
    with pytest.raises(ValueError, match="unique and increasing"):
        summarize_rows(duplicate, blocks=2)
