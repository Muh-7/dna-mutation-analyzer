"""Tests for high-level TF-binding analysis."""

from pathlib import Path

import pytest

from src.tf_binding.analysis import (
    analyze_tf_binding_sequences,
)
from src.tf_binding.fimo_parser import FimoHit


def make_hit(
    motif_id: str,
    motif_name: str,
    sequence_name: str,
    score: float,
    matched_sequence: str,
) -> FimoHit:
    """Create one test FIMO hit."""
    return FimoHit(
        motif_id=motif_id,
        motif_alt_id=motif_name,
        sequence_name=sequence_name,
        start=3,
        stop=8,
        strand="+",
        score=score,
        p_value=1e-5,
        q_value=0.01,
        matched_sequence=matched_sequence,
    )


def test_analyze_tf_binding_sequences(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    motif_path = (
        tmp_path
        / "motifs.meme"
    )

    motif_path.write_text(
        "MEME version 4\n",
        encoding="utf-8",
    )

    fake_hits = [
        make_hit(
            motif_id="MA0595.1",
            motif_name="SREBF1",
            sequence_name="reference",
            score=13.8,
            matched_sequence="CCGAAC",
        ),
        make_hit(
            motif_id="MA1596.1",
            motif_name="ZNF460",
            sequence_name="mutated",
            score=12.9,
            matched_sequence="CCGCAC",
        ),
    ]

    monkeypatch.setattr(
        "src.tf_binding.analysis.run_fimo",
        lambda **kwargs: fake_hits,
    )

    result = analyze_tf_binding_sequences(
        reference_sequence="AACCGAACTG",
        mutated_sequence="AACCGCACTG",
        motif_database_path=motif_path,
    )

    assert result["status"] == "completed"
    assert result["mutation"]["position_1_based"] == 6
    assert result["total_fimo_hits"] == 2
    assert result["reference_hits"] == 1
    assert result["mutated_hits"] == 1
    assert result["counts"]["gained"] == 1
    assert result["counts"]["lost"] == 1
    assert result["changed_sites"] == 2
