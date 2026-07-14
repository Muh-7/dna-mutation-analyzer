"""High-level transcription-factor motif analysis."""

from pathlib import Path
from typing import Any

from src.preprocessing.mutation_detector import (
    detect_snv,
)
from src.tf_binding.fimo_parser import (
    filter_hits_by_sequence,
    filter_hits_overlapping_position,
)
from src.tf_binding.fimo_runner import run_fimo
from src.tf_binding.motif_comparator import (
    compare_tf_binding_hits,
    summarize_tf_binding_changes,
)


def analyze_tf_binding_sequences(
    reference_sequence: str,
    mutated_sequence: str,
    motif_database_path: str | Path,
    fimo_executable: str = "fimo",
    p_value_threshold: float = 1e-4,
) -> dict[str, Any]:
    """
    Run FIMO and compare motifs overlapping one SNV.
    """
    mutation = detect_snv(
        reference_sequence=reference_sequence,
        mutated_sequence=mutated_sequence,
        allow_n=False,
    )

    all_hits = run_fimo(
        reference_sequence=reference_sequence,
        mutated_sequence=mutated_sequence,
        motif_database_path=motif_database_path,
        fimo_executable=fimo_executable,
        p_value_threshold=p_value_threshold,
    )

    reference_hits = filter_hits_by_sequence(
        hits=all_hits,
        sequence_name="reference",
    )

    mutated_hits = filter_hits_by_sequence(
        hits=all_hits,
        sequence_name="mutated",
    )

    reference_overlapping = (
        filter_hits_overlapping_position(
            hits=reference_hits,
            position_1_based=(
                mutation.position_1_based
            ),
        )
    )

    mutated_overlapping = (
        filter_hits_overlapping_position(
            hits=mutated_hits,
            position_1_based=(
                mutation.position_1_based
            ),
        )
    )

    changes = compare_tf_binding_hits(
        reference_hits=reference_hits,
        mutated_hits=mutated_hits,
        mutation_position_1_based=(
            mutation.position_1_based
        ),
    )

    comparison = summarize_tf_binding_changes(
        changes
    )

    return {
        "status": "completed",
        "tool": "FIMO",
        "analysis_type": (
            "reference_vs_mutated_motif_comparison"
        ),
        "motif_database": str(
            Path(motif_database_path)
        ),
        "threshold_type": "p_value",
        "p_value_threshold": float(
            p_value_threshold
        ),
        "sequence_length": mutation.sequence_length,
        "mutation": {
            "position_0_based": (
                mutation.position_0_based
            ),
            "position_1_based": (
                mutation.position_1_based
            ),
            "reference_base": (
                mutation.reference_base
            ),
            "alternate_base": (
                mutation.alternate_base
            ),
        },
        "total_fimo_hits": len(all_hits),
        "reference_hits": len(reference_hits),
        "mutated_hits": len(mutated_hits),
        "reference_hits_overlapping_mutation": len(
            reference_overlapping
        ),
        "mutated_hits_overlapping_mutation": len(
            mutated_overlapping
        ),
        "counts": comparison["counts"],
        "total_compared_sites": comparison[
            "total_compared_sites"
        ],
        "changed_sites": comparison[
            "changed_sites"
        ],
        "changes": comparison["changes"],
        "scientific_note": comparison[
            "scientific_note"
        ],
    }
