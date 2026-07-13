"""Tests for mutation-window extraction."""

import pytest

from src.preprocessing.window_extractor import (
    MutationWindow,
    extract_mutation_window,
)


def test_extract_centered_mutation_window() -> None:
    result = extract_mutation_window(
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
        flank_size=2,
    )

    assert isinstance(result, MutationWindow)
    assert result.reference_window == "CGTAC"
    assert result.mutated_window == "CGGAC"
    assert result.start_0_based == 3
    assert result.end_0_based_exclusive == 8
    assert result.mutation_position_0_based == 5
    assert result.mutation_position_1_based == 6
    assert result.mutation_index_in_window == 2
    assert result.left_flank_length == 2
    assert result.right_flank_length == 2


def test_extract_window_at_first_position() -> None:
    result = extract_mutation_window(
        reference_sequence="ACGTAC",
        mutated_sequence="GCGTAC",
        flank_size=3,
    )

    assert result.reference_window == "ACGT"
    assert result.mutated_window == "GCGT"
    assert result.start_0_based == 0
    assert result.end_0_based_exclusive == 4
    assert result.mutation_index_in_window == 0
    assert result.left_flank_length == 0
    assert result.right_flank_length == 3


def test_extract_window_at_last_position() -> None:
    result = extract_mutation_window(
        reference_sequence="ACGT",
        mutated_sequence="ACGA",
        flank_size=3,
    )

    assert result.reference_window == "ACGT"
    assert result.mutated_window == "ACGA"
    assert result.mutation_index_in_window == 3
    assert result.left_flank_length == 3
    assert result.right_flank_length == 0


def test_zero_flank_returns_only_mutation_base() -> None:
    result = extract_mutation_window(
        reference_sequence="AACCGT",
        mutated_sequence="AACCAT",
        flank_size=0,
    )

    assert result.reference_window == "G"
    assert result.mutated_window == "A"
    assert result.mutation_index_in_window == 0
    assert result.left_flank_length == 0
    assert result.right_flank_length == 0


def test_large_flank_returns_complete_sequence() -> None:
    result = extract_mutation_window(
        reference_sequence="ACGT",
        mutated_sequence="AGGT",
        flank_size=100,
    )

    assert result.reference_window == "ACGT"
    assert result.mutated_window == "AGGT"
    assert result.start_0_based == 0
    assert result.end_0_based_exclusive == 4


def test_sequences_are_cleaned_before_extraction() -> None:
    result = extract_mutation_window(
        reference_sequence="aa ccg tacgt",
        mutated_sequence="AACCG\nGACGT",
        flank_size=2,
    )

    assert result.reference_window == "CGTAC"
    assert result.mutated_window == "CGGAC"


def test_negative_flank_is_rejected() -> None:
    with pytest.raises(ValueError, match="cannot be negative"):
        extract_mutation_window(
            reference_sequence="ACGT",
            mutated_sequence="AGGT",
            flank_size=-1,
        )


def test_non_integer_flank_is_rejected() -> None:
    with pytest.raises(TypeError, match="must be an integer"):
        extract_mutation_window(
            reference_sequence="ACGT",
            mutated_sequence="AGGT",
            flank_size=2.5,  # type: ignore[arg-type]
        )


def test_identical_sequences_are_rejected() -> None:
    with pytest.raises(ValueError, match="identical"):
        extract_mutation_window(
            reference_sequence="ACGT",
            mutated_sequence="ACGT",
            flank_size=2,
        )


def test_multiple_mutations_are_rejected() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        extract_mutation_window(
            reference_sequence="ACGT",
            mutated_sequence="TGGT",
            flank_size=2,
        )