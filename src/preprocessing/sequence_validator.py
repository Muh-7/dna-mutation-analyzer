"""Utilities for cleaning and validating DNA sequences."""

# القواعد المقبولة في تسلسل DNA.
# N تعني أن القاعدة غير معروفة أو غير محددة.
VALID_DNA_BASES = frozenset({"A", "C", "G", "T", "N"})


def clean_sequence(sequence: str) -> str:
    """
    Clean a DNA sequence.

    The function:
    1. Removes spaces, tabs, and line breaks.
    2. Converts all characters to uppercase.

    Parameters
    ----------
    sequence:
        The original DNA sequence.

    Returns
    -------
    str
        The cleaned DNA sequence.

    Raises
    ------
    TypeError
        If sequence is not a string.
    """
    if not isinstance(sequence, str):
        raise TypeError("DNA sequence must be a string.")

    return "".join(sequence.split()).upper()


def validate_sequence(sequence: str, allow_n: bool = True) -> str:
    """
    Clean and validate a DNA sequence.

    Parameters
    ----------
    sequence:
        The DNA sequence to validate.

    allow_n:
        Whether the unknown nucleotide N is allowed.

    Returns
    -------
    str
        The cleaned and validated DNA sequence.

    Raises
    ------
    ValueError
        If the sequence is empty or contains invalid characters.
    TypeError
        If sequence is not a string.
    """
    cleaned_sequence = clean_sequence(sequence)

    if not cleaned_sequence:
        raise ValueError("DNA sequence cannot be empty.")

    allowed_bases = VALID_DNA_BASES if allow_n else frozenset(
        {"A", "C", "G", "T"}
    )

    invalid_bases = sorted(set(cleaned_sequence) - allowed_bases)

    if invalid_bases:
        invalid_text = ", ".join(invalid_bases)

        raise ValueError(
            f"DNA sequence contains invalid characters: {invalid_text}. "
            f"Allowed characters are: {', '.join(sorted(allowed_bases))}."
        )

    return cleaned_sequence