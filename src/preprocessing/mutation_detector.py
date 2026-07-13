"""Detection of single-nucleotide variants between DNA sequences."""

from dataclasses import dataclass

from src.preprocessing.sequence_validator import validate_sequence


@dataclass(frozen=True)
class MutationRecord:
    """
    Store information about a single-nucleotide variant.

    Attributes
    ----------
    mutation_type:
        The detected mutation type. Currently always "SNV".

    position_0_based:
        Mutation position using Python indexing.

    position_1_based:
        Mutation position using biological display indexing.

    reference_base:
        The nucleotide in the reference sequence.

    alternate_base:
        The nucleotide in the mutated sequence.

    sequence_length:
        Length of both DNA sequences.
    """

    mutation_type: str
    position_0_based: int
    position_1_based: int
    reference_base: str
    alternate_base: str
    sequence_length: int


def detect_snv(
    reference_sequence: str,
    mutated_sequence: str,
    allow_n: bool = False,
) -> MutationRecord:
    """
    Detect one single-nucleotide variant between two DNA sequences.

    Both sequences are cleaned and validated before comparison.

    Parameters
    ----------
    reference_sequence:
        The normal/reference DNA sequence.

    mutated_sequence:
        The mutated DNA sequence.

    allow_n:
        Whether the unknown nucleotide N is allowed in the sequences.

    Returns
    -------
    MutationRecord
        Information about the detected SNV.

    Raises
    ------
    ValueError
        If the sequences have different lengths, are identical,
        contain multiple differences, or contain N at the mutation site.

    TypeError
        If either sequence is not a string.
    """
    reference = validate_sequence(reference_sequence, allow_n=allow_n)
    mutated = validate_sequence(mutated_sequence, allow_n=allow_n)

    if len(reference) != len(mutated):
        raise ValueError(
            "Reference and mutated sequences must have the same length "
            "for SNV detection."
        )

    differences = [
        index
        for index, (reference_base, mutated_base) in enumerate(
            zip(reference, mutated)
        )
        if reference_base != mutated_base
    ]

    if not differences:
        raise ValueError(
            "No mutation was detected because the sequences are identical."
        )

    if len(differences) > 1:
        raise ValueError(
            "SNV detection requires exactly one nucleotide difference, "
            f"but {len(differences)} differences were found."
        )

    mutation_position = differences[0]
    reference_base = reference[mutation_position]
    alternate_base = mutated[mutation_position]

    if "N" in {reference_base, alternate_base}:
        raise ValueError(
            "The mutation site cannot contain N because the original "
            "nucleotide is unknown."
        )

    return MutationRecord(
        mutation_type="SNV",
        position_0_based=mutation_position,
        position_1_based=mutation_position + 1,
        reference_base=reference_base,
        alternate_base=alternate_base,
        sequence_length=len(reference),
    )