"""Tests for local reference-genome sequence handling."""

from pathlib import Path

import pytest

from src.reference.reference_genome import (
    apply_snv,
    build_variant_sequences,
    fetch_reference_base,
    fetch_reference_window,
    normalize_chromosome,
    open_reference_genome,
)


@pytest.fixture()
def tiny_reference_fasta(
    tmp_path: Path,
) -> Path:
    """Create a tiny FASTA file for isolated tests."""
    fasta_path = tmp_path / "tiny_reference.fa"

    fasta_path.write_text(
        ">chr1\n"
        "AACCGTACGT\n"
        ">chrM\n"
        "ATGC\n",
        encoding="utf-8",
    )

    return fasta_path


def test_normalize_chromosome() -> None:
    assert normalize_chromosome("1") == "chr1"
    assert normalize_chromosome("chr22") == "chr22"
    assert normalize_chromosome("x") == "chrX"
    assert normalize_chromosome("MT") == "chrM"


def test_fetch_reference_base(
    tiny_reference_fasta: Path,
) -> None:
    reference = open_reference_genome(
        tiny_reference_fasta
    )

    try:
        base = fetch_reference_base(
            reference=reference,
            chromosome="1",
            position_1_based=6,
        )
    finally:
        reference.close()

    assert base == "T"


def test_fetch_centered_reference_window(
    tiny_reference_fasta: Path,
) -> None:
    reference = open_reference_genome(
        tiny_reference_fasta
    )

    try:
        sequence, start, end, mutation_index = (
            fetch_reference_window(
                reference=reference,
                chromosome="chr1",
                position_1_based=6,
                flank_size=2,
            )
        )
    finally:
        reference.close()

    assert sequence == "CGTAC"
    assert start == 4
    assert end == 8
    assert mutation_index == 2


def test_window_near_chromosome_start(
    tiny_reference_fasta: Path,
) -> None:
    reference = open_reference_genome(
        tiny_reference_fasta
    )

    try:
        sequence, start, end, mutation_index = (
            fetch_reference_window(
                reference=reference,
                chromosome="chr1",
                position_1_based=1,
                flank_size=3,
            )
        )
    finally:
        reference.close()

    assert sequence == "AACC"
    assert start == 1
    assert end == 4
    assert mutation_index == 0


def test_apply_snv() -> None:
    result = apply_snv(
        reference_sequence="CGTAC",
        mutation_index_0_based=2,
        alternate_allele="G",
    )

    assert result == "CGGAC"


def test_build_variant_sequences(
    tiny_reference_fasta: Path,
) -> None:
    result = build_variant_sequences(
        fasta_path=tiny_reference_fasta,
        chromosome="1",
        position_1_based=6,
        reference_allele="T",
        alternate_allele="G",
        flank_size=2,
    )

    assert result.chromosome == "chr1"
    assert result.reference_sequence == "CGTAC"
    assert result.mutated_sequence == "CGGAC"
    assert result.reference_allele == "T"
    assert result.alternate_allele == "G"
    assert result.window_start_1_based == 4
    assert result.window_end_1_based == 8
    assert result.mutation_index_0_based == 2


def test_reference_allele_mismatch_is_rejected(
    tiny_reference_fasta: Path,
) -> None:
    with pytest.raises(
        ValueError,
        match="Reference allele mismatch",
    ):
        build_variant_sequences(
            fasta_path=tiny_reference_fasta,
            chromosome="1",
            position_1_based=6,
            reference_allele="A",
            alternate_allele="G",
            flank_size=2,
        )


def test_same_alternate_allele_is_rejected(
    tiny_reference_fasta: Path,
) -> None:
    with pytest.raises(
        ValueError,
        match="must differ",
    ):
        build_variant_sequences(
            fasta_path=tiny_reference_fasta,
            chromosome="1",
            position_1_based=6,
            reference_allele="T",
            alternate_allele="T",
            flank_size=2,
        )


def test_missing_chromosome_is_rejected(
    tiny_reference_fasta: Path,
) -> None:
    with pytest.raises(
        ValueError,
        match="was not found",
    ):
        build_variant_sequences(
            fasta_path=tiny_reference_fasta,
            chromosome="22",
            position_1_based=1,
            reference_allele="A",
            alternate_allele="G",
        )