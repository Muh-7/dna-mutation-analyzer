"""Tests for AlphaGenome splicing score utilities."""

import pandas as pd
import pytest

from src.models.alphagenome_client import (
    compute_merged_splicing_score,
)


def test_compute_merged_splicing_score() -> None:
    scores = pd.DataFrame(
        {
            "output_type": [
                "SPLICE_SITES",
                "SPLICE_SITES",
                "SPLICE_SITE_USAGE",
                "SPLICE_JUNCTIONS",
            ],
            "raw_score": [
                0.6,
                0.2,
                0.4,
                2.5,
            ],
        }
    )

    result = compute_merged_splicing_score(
        scores
    )

    assert result["splice_sites"] == pytest.approx(
        0.6
    )

    assert result[
        "splice_site_usage"
    ] == pytest.approx(0.4)

    assert result[
        "splice_junctions"
    ] == pytest.approx(2.5)

    assert result[
        "merged_splicing_score"
    ] == pytest.approx(1.5)


def test_merged_score_uses_absolute_values() -> None:
    scores = pd.DataFrame(
        {
            "output_type": [
                "SPLICE_SITES",
                "SPLICE_SITE_USAGE",
                "SPLICE_JUNCTIONS",
            ],
            "raw_score": [
                -0.8,
                -0.3,
                -1.5,
            ],
        }
    )

    result = compute_merged_splicing_score(
        scores
    )

    assert result["splice_sites"] == pytest.approx(
        0.8
    )

    assert result[
        "splice_site_usage"
    ] == pytest.approx(0.3)

    assert result[
        "splice_junctions"
    ] == pytest.approx(1.5)

    assert result[
        "merged_splicing_score"
    ] == pytest.approx(1.4)


def test_missing_required_columns_is_rejected() -> None:
    scores = pd.DataFrame(
        {
            "output_type": [
                "SPLICE_SITES"
            ]
        }
    )

    with pytest.raises(
        ValueError,
        match="missing columns",
    ):
        compute_merged_splicing_score(
            scores
        )