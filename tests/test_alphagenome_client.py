"""Tests for AlphaGenome client utilities."""

import pandas as pd
import pytest

from src.models.alphagenome_client import (
    build_alphagenome_variant,
    build_alphagenome_variant_from_input,
    build_prediction_interval,
    filter_expression_scores,
    summarize_expression_scores,
)
from src.schemas.variant_input import VariantInput


def test_build_alphagenome_variant() -> None:
    variant = build_alphagenome_variant(
        chromosome="chr22",
        position_1_based=36201698,
        reference_bases="a",
        alternate_bases="c",
    )

    assert variant.chromosome == "chr22"
    assert variant.position == 36201698
    assert variant.reference_bases == "A"
    assert variant.alternate_bases == "C"


def test_position_must_be_positive() -> None:
    with pytest.raises(
        ValueError,
        match="greater than or equal to 1",
    ):
        build_alphagenome_variant(
            chromosome="chr22",
            position_1_based=0,
            reference_bases="A",
            alternate_bases="C",
        )


def test_build_variant_from_project_input() -> None:
    variant_input = VariantInput(
        analysis_mode="genomic",
        reference_sequence="AACCGAACTG",
        mutated_sequence="AACCGCACTG",
        chromosome="22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="+",
    )

    variant = build_alphagenome_variant_from_input(
        variant_input
    )

    assert variant.chromosome == "chr22"
    assert variant.position == 36201698
    assert variant.reference_bases == "A"
    assert variant.alternate_bases == "C"


def test_raw_mode_is_rejected() -> None:
    variant_input = VariantInput(
        reference_sequence="ACGT",
        mutated_sequence="AGGT",
    )

    with pytest.raises(
        ValueError,
        match="requires genomic mode",
    ):
        build_alphagenome_variant_from_input(
            variant_input
        )


def test_prediction_interval_contains_variant() -> None:
    variant = build_alphagenome_variant(
        chromosome="chr22",
        position_1_based=36201698,
        reference_bases="A",
        alternate_bases="C",
    )

    interval = build_prediction_interval(variant)

    assert interval.width > 0
    assert variant.reference_overlaps(interval)


def test_filter_expression_scores() -> None:
    scores = pd.DataFrame(
        {
            "output_type": [
                "RNA_SEQ",
                "RNA_SEQ",
                "ATAC",
            ],
            "gene_name": [
                "APOL4",
                "APOL4",
                None,
            ],
            "biosample_name": [
                "colon tissue",
                "brain tissue",
                "colon tissue",
            ],
            "gtex_tissue": [
                "Colon",
                "Brain",
                None,
            ],
            "title": [
                "Colon RNA-seq",
                "Brain RNA-seq",
                "Colon ATAC",
            ],
            "raw_score": [
                -0.8,
                0.2,
                2.0,
            ],
            "quantile_score": [
                -0.92,
                0.40,
                0.99,
            ],
        }
    )

    result = filter_expression_scores(
        scores=scores,
        gene_name="APOL4",
        tissue_query="colon",
    )

    assert len(result) == 1
    assert result.iloc[0]["raw_score"] == -0.8


def test_summarize_expression_scores() -> None:
    scores = pd.DataFrame(
        {
            "gene_name": ["APOL4", "APOL4"],
            "biosample_name": ["colon", "colon"],
            "ontology_curie": [
                "UBERON:0001157",
                "UBERON:0001157",
            ],
            "track_name": ["track-1", "track-2"],
            "raw_score": [-0.7, 0.2],
            "quantile_score": [-0.95, 0.30],
        }
    )

    summary = summarize_expression_scores(scores)

    assert summary["status"] == "completed"
    assert summary["matched_tracks"] == 2
    assert summary["direction"] == "predicted_decrease"
    assert summary["raw_score"] == -0.7


def test_empty_score_summary() -> None:
    summary = summarize_expression_scores(
        pd.DataFrame()
    )

    assert summary["status"] == (
        "no_matching_expression_tracks"
    )
    assert summary["matched_tracks"] == 0