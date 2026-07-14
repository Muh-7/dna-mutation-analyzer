"""Verify one real variant against the local GRCh38 reference."""

from dataclasses import asdict
from pprint import pprint

from src.reference.reference_genome import (
    build_variant_sequences,
)


FASTA_PATH = (
    "data/reference/"
    "GRCh38.primary_assembly.genome.fa"
)


def main() -> None:
    """Verify the APOL4 example variant."""
    result = build_variant_sequences(
        fasta_path=FASTA_PATH,
        chromosome="chr22",
        position_1_based=36201698,
        reference_allele="A",
        alternate_allele="C",
        flank_size=5,
    )

    pprint(asdict(result))

    print("\nReference sequence:")
    print(result.reference_sequence)

    print("\nMutated sequence:")
    print(result.mutated_sequence)

    print("\nMutation visualization:")
    print(result.reference_sequence)
    print(" " * result.mutation_index_0_based + "^")
    print(
        f"{result.reference_allele}"
        f" → {result.alternate_allele}"
    )


if __name__ == "__main__":
    main()