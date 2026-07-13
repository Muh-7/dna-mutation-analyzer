"""Tests for the DNA variant input schema."""

import pytest
from pydantic import ValidationError

from src.schemas.variant_input import (
    AnalysisMode,
    GenomeBuild,
    Strand,
    VariantInput,
)


def test_raw_sequence_input_is_valid() -> None:
    data = VariantInput(
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
    )

    assert data.analysis_mode == AnalysisMode.RAW_SEQUENCE
    assert data.genome_build == GenomeBuild.HG38
    assert data.reference_sequence == "AACCGTACGT"
    assert data.mutated_sequence == "AACCGGACGT"
    assert data.flank_size == 100


def test_sequences_are_cleaned() -> None:
    data = VariantInput(
        reference_sequence="aac cg tacgt",
        mutated_sequence="AACCG\nGACGT",
    )

    assert data.reference_sequence == "AACCGTACGT"
    assert data.mutated_sequence == "AACCGGACGT"


def test_valid_genomic_input() -> None:
    data = VariantInput(
        analysis_mode="genomic",
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
        chromosome="22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="+",
    )

    assert data.analysis_mode == AnalysisMode.GENOMIC
    assert data.chromosome == "chr22"
    assert data.genomic_position == 36201698
    assert data.strand == Strand.FORWARD


def test_chromosome_is_normalized() -> None:
    data = VariantInput(
        reference_sequence="ACGT",
        mutated_sequence="AGGT",
        chromosome="x",
    )

    assert data.chromosome == "chrX"


def test_identical_sequences_are_rejected() -> None:
    with pytest.raises(ValidationError, match="identical"):
        VariantInput(
            reference_sequence="ACGT",
            mutated_sequence="ACGT",
        )


def test_multiple_mutations_are_rejected() -> None:
    with pytest.raises(ValidationError, match="exactly one"):
        VariantInput(
            reference_sequence="ACGT",
            mutated_sequence="TGGT",
        )


def test_invalid_chromosome_is_rejected() -> None:
    with pytest.raises(ValidationError, match="Chromosome"):
        VariantInput(
            reference_sequence="ACGT",
            mutated_sequence="AGGT",
            chromosome="chr90",
        )


def test_genomic_mode_requires_metadata() -> None:
    with pytest.raises(
        ValidationError,
        match="Genomic analysis mode requires",
    ):
        VariantInput(
            analysis_mode="genomic",
            reference_sequence="ACGT",
            mutated_sequence="AGGT",
        )


def test_negative_flank_size_is_rejected() -> None:
    with pytest.raises(ValidationError):
        VariantInput(
            reference_sequence="ACGT",
            mutated_sequence="AGGT",
            flank_size=-1,
        )


def test_unknown_fields_are_rejected() -> None:
    with pytest.raises(ValidationError, match="Extra inputs"):
        VariantInput(
            reference_sequence="ACGT",
            mutated_sequence="AGGT",
            unknown_field="value",  # type: ignore[call-arg]
        )