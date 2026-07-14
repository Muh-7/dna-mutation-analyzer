"""Compare APOL4 reference and mutated FIMO motif hits."""

from pprint import pprint

from src.tf_binding.fimo_parser import (
    filter_hits_by_sequence,
    parse_fimo_tsv,
)
from src.tf_binding.motif_comparator import (
    compare_tf_binding_hits,
    summarize_tf_binding_changes,
)


FIMO_RESULTS_PATH = "outputs/fimo_apol4/fimo.tsv"
MUTATION_POSITION_IN_SEQUENCE = 51


def main() -> None:
    """Parse and compare the APOL4 FIMO results."""
    all_hits = parse_fimo_tsv(
        FIMO_RESULTS_PATH
    )

    reference_hits = filter_hits_by_sequence(
        all_hits,
        sequence_name="reference",
    )

    mutated_hits = filter_hits_by_sequence(
        all_hits,
        sequence_name="mutated",
    )

    changes = compare_tf_binding_hits(
        reference_hits=reference_hits,
        mutated_hits=mutated_hits,
        mutation_position_1_based=(
            MUTATION_POSITION_IN_SEQUENCE
        ),
    )

    summary = summarize_tf_binding_changes(
        changes
    )

    print("\nTF-binding motif comparison:\n")
    pprint(summary)


if __name__ == "__main__":
    main()
