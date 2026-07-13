"""Tests for single-nucleotide variant detection."""

import pytest

from src.preprocessing.mutation_detector import (
    MutationRecord,
    detect_snv,
)


def test_detect_snv_returns_correct_mutation() -> None:
    result = detect_snv(
        reference_sequence="ACGTAC",
        mutated_sequence="ACGGAC",
    )

    assert isinstance(result, MutationRecord)
    assert result.mutation_type == "SNV"
    assert result.position_0_based == 3
    assert result.position_1_based == 4
    assert result.reference_base == "T"
    assert result.alternate_base == "G"
    assert result.sequence_length == 6


def test_detect_snv_cleans_sequences_before_comparison() -> None:
    result = detect_snv(
        reference_sequence="acg tac",
        mutated_sequence="ACG\nGAC",
    )

    assert result.position_1_based == 4
    assert result.reference_base == "T"
    assert result.alternate_base == "G"


def test_detect_snv_detects_mutation_at_first_position() -> None:
    result = detect_snv(
        reference_sequence="ACGT",
        mutated_sequence="GCGT",
    )

    assert result.position_0_based == 0
    assert result.position_1_based == 1
    assert result.reference_base == "A"
    assert result.alternate_base == "G"


def test_detect_snv_detects_mutation_at_last_position() -> None:
    result = detect_snv(
        reference_sequence="ACGT",
        mutated_sequence="ACGA",
    )

    assert result.position_0_based == 3
    assert result.position_1_based == 4
    assert result.reference_base == "T"
    assert result.alternate_base == "A"


def test_detect_snv_rejects_identical_sequences() -> None:
    with pytest.raises(ValueError, match="identical"):
        detect_snv("ACGT", "ACGT")


def test_detect_snv_rejects_different_lengths() -> None:
    with pytest.raises(ValueError, match="same length"):
        detect_snv("ACGT", "ACGTA")


def test_detect_snv_rejects_multiple_differences() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        detect_snv("ACGT", "TGGT")


def test_detect_snv_rejects_invalid_characters() -> None:
    with pytest.raises(ValueError, match="invalid characters"):
        detect_snv("ACGT", "ACGX")


def test_detect_snv_rejects_n_by_default() -> None:
    with pytest.raises(ValueError, match="invalid characters"):
        detect_snv("ACNT", "ACGT")


def test_detect_snv_rejects_n_at_mutation_site() -> None:
    with pytest.raises(ValueError, match="mutation site cannot contain N"):
        detect_snv(
            reference_sequence="ACNT",
            mutated_sequence="ACGT",
            allow_n=True,
        )