"""Tests for reference-versus-mutated motif comparison."""

from src.tf_binding.fimo_parser import FimoHit
from src.tf_binding.motif_comparator import (
    TFBindingChangeType,
    compare_tf_binding_hits,
    summarize_tf_binding_changes,
)


def make_hit(
    motif_id: str,
    motif_name: str,
    sequence_name: str,
    start: int,
    stop: int,
    score: float,
    matched_sequence: str,
) -> FimoHit:
    return FimoHit(
        motif_id=motif_id,
        motif_alt_id=motif_name,
        sequence_name=sequence_name,
        start=start,
        stop=stop,
        strand="+",
        score=score,
        p_value=1e-5,
        q_value=0.01,
        matched_sequence=matched_sequence,
    )


def test_detects_gained_and_lost_sites() -> None:
    reference_hits = [
        make_hit(
            motif_id="MA0595.1",
            motif_name="SREBF1",
            sequence_name="reference",
            start=48,
            stop=57,
            score=13.8293,
            matched_sequence="CTCACCCGAC",
        )
    ]

    mutated_hits = [
        make_hit(
            motif_id="MA1596.1",
            motif_name="ZNF460",
            sequence_name="mutated",
            start=40,
            stop=55,
            score=12.9213,
            matched_sequence="GCCTCGGACTCCCCCG",
        )
    ]

    changes = compare_tf_binding_hits(
        reference_hits=reference_hits,
        mutated_hits=mutated_hits,
        mutation_position_1_based=51,
    )

    assert len(changes) == 2

    change_types = {
        change.change_type
        for change in changes
    }

    assert TFBindingChangeType.GAINED in change_types
    assert TFBindingChangeType.LOST in change_types


def test_detects_score_strengthening() -> None:
    reference_hit = make_hit(
        motif_id="MA0001.1",
        motif_name="TF1",
        sequence_name="reference",
        start=48,
        stop=57,
        score=8.0,
        matched_sequence="CTCACCCGAC",
    )

    mutated_hit = make_hit(
        motif_id="MA0001.1",
        motif_name="TF1",
        sequence_name="mutated",
        start=48,
        stop=57,
        score=10.0,
        matched_sequence="CTCCCCCGAC",
    )

    changes = compare_tf_binding_hits(
        reference_hits=[reference_hit],
        mutated_hits=[mutated_hit],
        mutation_position_1_based=51,
    )

    assert len(changes) == 1
    assert (
        changes[0].change_type
        == TFBindingChangeType.STRENGTHENED
    )
    assert changes[0].score_delta == 2.0


def test_summary_counts_changes() -> None:
    reference_hits = [
        make_hit(
            "MA0595.1",
            "SREBF1",
            "reference",
            48,
            57,
            13.8,
            "CTCACCCGAC",
        )
    ]

    mutated_hits = [
        make_hit(
            "MA1596.1",
            "ZNF460",
            "mutated",
            40,
            55,
            12.9,
            "GCCTCGGACTCCCCCG",
        )
    ]

    changes = compare_tf_binding_hits(
        reference_hits=reference_hits,
        mutated_hits=mutated_hits,
        mutation_position_1_based=51,
    )

    summary = summarize_tf_binding_changes(changes)

    assert summary["status"] == "completed"
    assert summary["counts"]["gained"] == 1
    assert summary["counts"]["lost"] == 1
    assert summary["changed_sites"] == 2
