"""Complete preprocessing pipeline for one DNA SNV."""

from dataclasses import asdict
from typing import Any

from src.preprocessing.mutation_detector import detect_snv
from src.preprocessing.sequence_validator import validate_sequence
from src.preprocessing.window_extractor import extract_mutation_window


def prepare_variant_input(
    reference_sequence: str,
    mutated_sequence: str,
    flank_size: int = 100,
    allow_n: bool = False,
) -> dict[str, Any]:
    """
    Validate two DNA sequences, detect one SNV, and extract its window.

    Parameters
    ----------
    reference_sequence:
        The normal/reference DNA sequence.

    mutated_sequence:
        The DNA sequence containing one nucleotide substitution.

    flank_size:
        Number of nucleotides requested before and after the mutation.

    allow_n:
        Whether N is allowed outside the mutation site.

    Returns
    -------
    dict
        A structured preprocessing result.
    """
    reference = validate_sequence(
        reference_sequence,
        allow_n=allow_n,
    )

    mutated = validate_sequence(
        mutated_sequence,
        allow_n=allow_n,
    )

    mutation = detect_snv(
        reference_sequence=reference,
        mutated_sequence=mutated,
        allow_n=allow_n,
    )

    window = extract_mutation_window(
        reference_sequence=reference,
        mutated_sequence=mutated,
        flank_size=flank_size,
        allow_n=allow_n,
    )

    return {
        "reference_sequence": reference,
        "mutated_sequence": mutated,
        "sequence_length": len(reference),
        "flank_size_requested": flank_size,
        "mutation": asdict(mutation),
        "window": asdict(window),
    }