"""Tests for the complete DNA preprocessing pipeline."""

import pytest

from src.preprocessing.pipeline import prepare_variant_input
from src.preprocessing.pipeline import (
    prepare_variant_input,
    prepare_variant_model,
)
from src.schemas.variant_input import VariantInput



def test_prepare_variant_input_returns_complete_result() -> None:
    result = prepare_variant_input(
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
        flank_size=2,
    )

    assert result["reference_sequence"] == "AACCGTACGT"
    assert result["mutated_sequence"] == "AACCGGACGT"
    assert result["sequence_length"] == 10
    assert result["flank_size_requested"] == 2

    mutation = result["mutation"]

    assert mutation["mutation_type"] == "SNV"
    assert mutation["position_0_based"] == 5
    assert mutation["position_1_based"] == 6
    assert mutation["reference_base"] == "T"
    assert mutation["alternate_base"] == "G"

    window = result["window"]

    assert window["reference_window"] == "CGTAC"
    assert window["mutated_window"] == "CGGAC"
    assert window["mutation_index_in_window"] == 2


def test_pipeline_cleans_input_sequences() -> None:
    result = prepare_variant_input(
        reference_sequence="aac cg tacgt",
        mutated_sequence="AACCG\nGACGT",
        flank_size=2,
    )

    assert result["reference_sequence"] == "AACCGTACGT"
    assert result["mutated_sequence"] == "AACCGGACGT"


def test_pipeline_rejects_identical_sequences() -> None:
    with pytest.raises(ValueError, match="identical"):
        prepare_variant_input(
            reference_sequence="ACGT",
            mutated_sequence="ACGT",
        )


def test_pipeline_rejects_multiple_mutations() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        prepare_variant_input(
            reference_sequence="ACGT",
            mutated_sequence="TGGT",
        )


def test_pipeline_rejects_different_lengths() -> None:
    with pytest.raises(ValueError, match="same length"):
        prepare_variant_input(
            reference_sequence="ACGT",
            mutated_sequence="ACGTA",
        )

def test_prepare_variant_model_includes_metadata() -> None:
    variant = VariantInput(
        analysis_mode="genomic",
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
        chromosome="22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="+",
        flank_size=2,
    )

    result = prepare_variant_model(variant)

    assert result["mutation"]["reference_base"] == "T"
    assert result["mutation"]["alternate_base"] == "G"
    assert result["window"]["reference_window"] == "CGTAC"
    assert result["window"]["mutated_window"] == "CGGAC"

    metadata = result["metadata"]

    assert metadata["analysis_mode"] == "genomic"
    assert metadata["genome_build"] == "hg38"
    assert metadata["chromosome"] == "chr22"
    assert metadata["genomic_position"] == 36201698
    assert metadata["gene_name"] == "APOL4"
    assert metadata["tissue"] == "colon"
    assert metadata["strand"] == "+"