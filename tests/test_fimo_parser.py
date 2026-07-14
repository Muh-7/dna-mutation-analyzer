"""Tests for parsing FIMO result files."""

from pathlib import Path

from src.tf_binding.fimo_parser import (
    filter_hits_by_sequence,
    filter_hits_overlapping_position,
    parse_fimo_tsv,
)


def test_parse_fimo_tsv(
    tmp_path: Path,
) -> None:
    result_path = tmp_path / "fimo.tsv"

    result_path.write_text(
        "motif_id\tmotif_alt_id\tsequence_name\t"
        "start\tstop\tstrand\tscore\tp-value\t"
        "q-value\tmatched_sequence\n"
        "MA0001.1\tTF1\treference\t48\t57\t+\t"
        "13.8\t1.7e-05\t0.005\tCTCACCCGAC\n"
        "# FIMO generated comment\n",
        encoding="utf-8",
    )

    hits = parse_fimo_tsv(result_path)

    assert len(hits) == 1

    hit = hits[0]

    assert hit.motif_id == "MA0001.1"
    assert hit.motif_alt_id == "TF1"
    assert hit.sequence_name == "reference"
    assert hit.start == 48
    assert hit.stop == 57
    assert hit.score == 13.8
    assert hit.overlaps_position(51)


def test_filter_fimo_hits(
    tmp_path: Path,
) -> None:
    result_path = tmp_path / "fimo.tsv"

    result_path.write_text(
        "motif_id\tmotif_alt_id\tsequence_name\t"
        "start\tstop\tstrand\tscore\tp-value\t"
        "q-value\tmatched_sequence\n"
        "MA1\tTF1\treference\t48\t57\t+\t"
        "13.8\t1e-05\t0.005\tCTCACCCGAC\n"
        "MA2\tTF2\tmutated\t10\t20\t+\t"
        "10.0\t2e-05\t0.010\tAAAAAAAAAAA\n",
        encoding="utf-8",
    )

    hits = parse_fimo_tsv(result_path)

    reference_hits = filter_hits_by_sequence(
        hits,
        "reference",
    )

    overlapping = filter_hits_overlapping_position(
        hits,
        position_1_based=51,
    )

    assert len(reference_hits) == 1
    assert len(overlapping) == 1
    assert overlapping[0].motif_alt_id == "TF1"
