"""Extract DNA sequence windows around a detected mutation."""

from dataclasses import dataclass

from src.preprocessing.mutation_detector import detect_snv
from src.preprocessing.sequence_validator import validate_sequence


@dataclass(frozen=True)
class MutationWindow:
    """
    Store reference and mutated sequence windows around an SNV.

    Attributes
    ----------
    reference_window:
        Extracted window from the reference sequence.

    mutated_window:
        Extracted window from the mutated sequence.

    start_0_based:
        Start position of the window in the original sequence.

    end_0_based_exclusive:
        End position of the window using Python's exclusive indexing.

    mutation_position_0_based:
        Mutation position in the complete sequence.

    mutation_position_1_based:
        Mutation position using biological indexing.

    mutation_index_in_window:
        Mutation position inside the extracted window.

    left_flank_length:
        Number of extracted nucleotides before the mutation.

    right_flank_length:
        Number of extracted nucleotides after the mutation.
    """

    reference_window: str
    mutated_window: str
    start_0_based: int
    end_0_based_exclusive: int
    mutation_position_0_based: int
    mutation_position_1_based: int
    mutation_index_in_window: int
    left_flank_length: int
    right_flank_length: int


def extract_mutation_window(
    reference_sequence: str,
    mutated_sequence: str,
    flank_size: int = 100,
    allow_n: bool = False,
) -> MutationWindow:
    """
    Extract reference and mutated DNA windows around one SNV.

    Parameters
    ----------
    reference_sequence:
        The normal/reference DNA sequence.

    mutated_sequence:
        The mutated DNA sequence.

    flank_size:
        Requested number of nucleotides before and after the mutation.

    allow_n:
        Whether unknown nucleotides represented by N are allowed.

    Returns
    -------
    MutationWindow
        Extracted windows and their position information.

    Raises
    ------
    TypeError
        If flank_size is not an integer.

    ValueError
        If flank_size is negative or SNV detection fails.
    """
    if isinstance(flank_size, bool) or not isinstance(flank_size, int):
        raise TypeError("Flank size must be an integer.")

    if flank_size < 0:
        raise ValueError("Flank size cannot be negative.")

    mutation = detect_snv(
        reference_sequence=reference_sequence,
        mutated_sequence=mutated_sequence,
        allow_n=allow_n,
    )

    reference = validate_sequence(reference_sequence, allow_n=allow_n)
    mutated = validate_sequence(mutated_sequence, allow_n=allow_n)

    mutation_position = mutation.position_0_based

    start = max(0, mutation_position - flank_size)

    end = min(
        len(reference),
        mutation_position + flank_size + 1,
    )

    reference_window = reference[start:end]
    mutated_window = mutated[start:end]

    mutation_index_in_window = mutation_position - start
    left_flank_length = mutation_index_in_window
    right_flank_length = end - mutation_position - 1

    return MutationWindow(
        reference_window=reference_window,
        mutated_window=mutated_window,
        start_0_based=start,
        end_0_based_exclusive=end,
        mutation_position_0_based=mutation_position,
        mutation_position_1_based=mutation.position_1_based,
        mutation_index_in_window=mutation_index_in_window,
        left_flank_length=left_flank_length,
        right_flank_length=right_flank_length,
    )