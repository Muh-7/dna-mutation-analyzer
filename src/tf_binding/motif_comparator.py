"""Compare FIMO motif matches between reference and mutated DNA."""

from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any

from src.tf_binding.fimo_parser import (
    FimoHit,
    filter_hits_overlapping_position,
)


class TFBindingChangeType(str, Enum):
    """Possible changes in predicted motif matching."""

    GAINED = "gained"
    LOST = "lost"
    STRENGTHENED = "strengthened"
    WEAKENED = "weakened"
    UNCHANGED = "unchanged"


@dataclass(frozen=True)
class TFBindingChange:
    """Describe one predicted motif-match change."""

    change_type: TFBindingChangeType

    motif_id: str
    motif_alt_id: str

    start: int
    stop: int
    strand: str

    reference_score: float | None
    mutated_score: float | None
    score_delta: float | None

    reference_p_value: float | None
    mutated_p_value: float | None

    reference_q_value: float | None
    mutated_q_value: float | None

    reference_matched_sequence: str | None
    mutated_matched_sequence: str | None


def _hit_key(
    hit: FimoHit,
) -> tuple[str, int, int, str]:
    """
    Build a key for matching the same motif site.

    Since the normal and mutated sequences have equal lengths,
    their coordinates are directly comparable.
    """
    return (
        hit.motif_id,
        hit.start,
        hit.stop,
        hit.strand,
    )


def _select_best_hits(
    hits: list[FimoHit],
) -> dict[tuple[str, int, int, str], FimoHit]:
    """
    Keep the most significant hit for every motif-position key.
    """
    selected: dict[
        tuple[str, int, int, str],
        FimoHit,
    ] = {}

    for hit in hits:
        key = _hit_key(hit)
        current = selected.get(key)

        if current is None or hit.p_value < current.p_value:
            selected[key] = hit

    return selected


def compare_tf_binding_hits(
    reference_hits: list[FimoHit],
    mutated_hits: list[FimoHit],
    mutation_position_1_based: int,
    minimum_score_delta: float = 1e-6,
) -> list[TFBindingChange]:
    """
    Compare reference and mutated motif hits around one SNV.

    Only motif sites overlapping the mutation are compared.
    """
    if minimum_score_delta < 0:
        raise ValueError(
            "Minimum score delta cannot be negative."
        )

    reference_overlapping = filter_hits_overlapping_position(
        hits=reference_hits,
        position_1_based=mutation_position_1_based,
    )

    mutated_overlapping = filter_hits_overlapping_position(
        hits=mutated_hits,
        position_1_based=mutation_position_1_based,
    )

    reference_by_key = _select_best_hits(
        reference_overlapping
    )

    mutated_by_key = _select_best_hits(
        mutated_overlapping
    )

    all_keys = sorted(
        set(reference_by_key) | set(mutated_by_key)
    )

    changes: list[TFBindingChange] = []

    for key in all_keys:
        reference_hit = reference_by_key.get(key)
        mutated_hit = mutated_by_key.get(key)

        motif_id, start, stop, strand = key

        if reference_hit is None:
            change_type = TFBindingChangeType.GAINED
            score_delta = None

        elif mutated_hit is None:
            change_type = TFBindingChangeType.LOST
            score_delta = None

        else:
            score_delta = (
                mutated_hit.score - reference_hit.score
            )

            if score_delta > minimum_score_delta:
                change_type = (
                    TFBindingChangeType.STRENGTHENED
                )

            elif score_delta < -minimum_score_delta:
                change_type = (
                    TFBindingChangeType.WEAKENED
                )

            else:
                change_type = TFBindingChangeType.UNCHANGED

        representative_hit = reference_hit or mutated_hit

        if representative_hit is None:
            raise RuntimeError(
                "A motif comparison has no representative hit."
            )

        changes.append(
            TFBindingChange(
                change_type=change_type,
                motif_id=motif_id,
                motif_alt_id=(
                    representative_hit.motif_alt_id
                ),
                start=start,
                stop=stop,
                strand=strand,
                reference_score=(
                    reference_hit.score
                    if reference_hit is not None
                    else None
                ),
                mutated_score=(
                    mutated_hit.score
                    if mutated_hit is not None
                    else None
                ),
                score_delta=score_delta,
                reference_p_value=(
                    reference_hit.p_value
                    if reference_hit is not None
                    else None
                ),
                mutated_p_value=(
                    mutated_hit.p_value
                    if mutated_hit is not None
                    else None
                ),
                reference_q_value=(
                    reference_hit.q_value
                    if reference_hit is not None
                    else None
                ),
                mutated_q_value=(
                    mutated_hit.q_value
                    if mutated_hit is not None
                    else None
                ),
                reference_matched_sequence=(
                    reference_hit.matched_sequence
                    if reference_hit is not None
                    else None
                ),
                mutated_matched_sequence=(
                    mutated_hit.matched_sequence
                    if mutated_hit is not None
                    else None
                ),
            )
        )

    priority = {
        TFBindingChangeType.GAINED: 0,
        TFBindingChangeType.LOST: 1,
        TFBindingChangeType.STRENGTHENED: 2,
        TFBindingChangeType.WEAKENED: 3,
        TFBindingChangeType.UNCHANGED: 4,
    }

    return sorted(
        changes,
        key=lambda change: (
            priority[change.change_type],
            change.motif_alt_id,
            change.start,
        ),
    )


def summarize_tf_binding_changes(
    changes: list[TFBindingChange],
) -> dict[str, Any]:
    """Create a serializable TF-binding comparison summary."""
    counts = {
        change_type.value: 0
        for change_type in TFBindingChangeType
    }

    for change in changes:
        counts[change.change_type.value] += 1

    serialized_changes = []

    for change in changes:
        record = asdict(change)
        record["change_type"] = change.change_type.value
        serialized_changes.append(record)

    biologically_relevant_changes = [
        record
        for record in serialized_changes
        if record["change_type"] != "unchanged"
    ]

    return {
        "status": "completed",
        "analysis_type": "motif_match_comparison",
        "counts": counts,
        "total_compared_sites": len(changes),
        "changed_sites": len(
            biologically_relevant_changes
        ),
        "changes": biologically_relevant_changes,
        "scientific_note": (
            "FIMO measures statistical DNA-to-motif sequence "
            "matching. These results do not directly demonstrate "
            "transcription-factor binding inside living cells."
        ),
    }
