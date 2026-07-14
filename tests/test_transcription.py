"""Tests for DNA-to-RNA transcription."""

import pytest

from src.rna.transcription import (
    RNAComparison,
    analyze_rna_change,
    analyze_variant_rna,
    dna_to_rna,
    reverse_complement_dna,
)
from src.schemas.variant_input import VariantInput


def test_reverse_complement_dna() -> None:
    result = reverse_complement_dna("ATGCC")

    assert result == "GGCAT"


def test_reverse_complement_cleans_sequence() -> None:
    result = reverse_complement_dna("at gcc")

    assert result == "GGCAT"


def test_forward_strand_transcription() -> None:
    result = dna_to_rna(
        sequence="ATGCC",
        strand="+",
    )

    assert result == "AUGCC"


def test_reverse_strand_transcription() -> None:
    result = dna_to_rna(
        sequence="ATGCC",
        strand="-",
    )

    assert result == "GGCAU"


def test_invalid_strand_is_rejected() -> None:
    with pytest.raises(ValueError, match="Strand"):
        dna_to_rna(
            sequence="ATGCC",
            strand="forward",
        )


def test_forward_strand_rna_change() -> None:
    result = analyze_rna_change(
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
        strand="+",
    )

    assert isinstance(result, RNAComparison)
    assert result.reference_rna == "AACCGUACGU"
    assert result.mutated_rna == "AACCGGACGU"
    assert result.position_0_based == 5
    assert result.position_1_based == 6
    assert result.reference_rna_base == "U"
    assert result.alternate_rna_base == "G"
    assert result.strand == "+"


def test_reverse_strand_rna_change() -> None:
    result = analyze_rna_change(
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
        strand="-",
    )

    assert result.reference_rna == "ACGUACGGUU"
    assert result.mutated_rna == "ACGUCCGGUU"
    assert result.position_0_based == 4
    assert result.position_1_based == 5
    assert result.reference_rna_base == "A"
    assert result.alternate_rna_base == "C"
    assert result.strand == "-"


def test_variant_model_rna_analysis() -> None:
    variant = VariantInput(
        analysis_mode="genomic",
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
        chromosome="22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="+",
    )

    result = analyze_variant_rna(variant)

    assert result.reference_rna == "AACCGUACGU"
    assert result.mutated_rna == "AACCGGACGU"
    assert result.position_1_based == 6


def test_raw_variant_without_strand_is_rejected() -> None:
    variant = VariantInput(
        reference_sequence="ACGT",
        mutated_sequence="AGGT",
    )

    with pytest.raises(
        ValueError,
        match="requires strand information",
    ):
        analyze_variant_rna(variant)


def test_n_is_rejected_by_default() -> None:
    with pytest.raises(ValueError, match="invalid characters"):
        dna_to_rna(
            sequence="ATGCN",
            strand="+",
        )


def test_n_can_be_allowed() -> None:
    result = dna_to_rna(
        sequence="ATGCN",
        strand="+",
        allow_n=True,
    )

    assert result == "AUGCN"