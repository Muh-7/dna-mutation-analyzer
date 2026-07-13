"""Tests for DNA sequence preprocessing."""

import pytest

from src.preprocessing.sequence_validator import (
    clean_sequence,
    validate_sequence,
)


def test_clean_sequence_converts_to_uppercase() -> None:
    result = clean_sequence("acgt")

    assert result == "ACGT"


def test_clean_sequence_removes_whitespace() -> None:
    sequence = "AC G\tT\nAAC"

    result = clean_sequence(sequence)

    assert result == "ACGTAAC"


def test_clean_sequence_rejects_non_string_input() -> None:
    with pytest.raises(TypeError, match="must be a string"):
        clean_sequence(12345)  # type: ignore[arg-type]


def test_validate_sequence_accepts_valid_dna() -> None:
    result = validate_sequence("ACGTACGT")

    assert result == "ACGTACGT"


def test_validate_sequence_accepts_n_by_default() -> None:
    result = validate_sequence("ACGTN")

    assert result == "ACGTN"


def test_validate_sequence_can_reject_n() -> None:
    with pytest.raises(ValueError, match="invalid characters"):
        validate_sequence("ACGTN", allow_n=False)


def test_validate_sequence_rejects_empty_sequence() -> None:
    with pytest.raises(ValueError, match="cannot be empty"):
        validate_sequence("   \n\t")


def test_validate_sequence_rejects_invalid_characters() -> None:
    with pytest.raises(ValueError, match="X"):
        validate_sequence("ACGTX")


def test_validate_sequence_cleans_before_validation() -> None:
    result = validate_sequence(" acg\nnt ")

    assert result == "ACGNT"