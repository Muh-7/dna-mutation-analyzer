"""DNA-to-RNA transcription and RNA mutation comparison."""

from dataclasses import dataclass

from src.preprocessing.mutation_detector import detect_snv
from src.preprocessing.sequence_validator import validate_sequence
from src.schemas.variant_input import VariantInput


DNA_COMPLEMENT_TABLE = str.maketrans(
    {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
    }
)


@dataclass(frozen=True)
class RNAComparison:
    """
    Store the RNA-level effect of one DNA SNV.

    Attributes
    ----------
    reference_rna:
        RNA transcribed from the reference DNA sequence.

    mutated_rna:
        RNA transcribed from the mutated DNA sequence.

    strand:
        Gene strand used during transcription.

    position_0_based:
        Changed RNA position using Python indexing.

    position_1_based:
        Changed RNA position using biological indexing.

    reference_rna_base:
        RNA nucleotide in the reference transcript.

    alternate_rna_base:
        RNA nucleotide in the mutated transcript.

    sequence_length:
        Length of both RNA sequences.
    """

    reference_rna: str
    mutated_rna: str
    strand: str
    position_0_based: int
    position_1_based: int
    reference_rna_base: str
    alternate_rna_base: str
    sequence_length: int


def reverse_complement_dna(
    sequence: str,
    allow_n: bool = False,
) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Example
    -------
    ATGC becomes GCAT.
    """
    cleaned_sequence = validate_sequence(
        sequence,
        allow_n=allow_n,
    )

    complement = cleaned_sequence.translate(DNA_COMPLEMENT_TABLE)

    return complement[::-1]


def dna_to_rna(
    sequence: str,
    strand: str = "+",
    allow_n: bool = False,
) -> str:
    """
    Transcribe genomic DNA into RNA according to gene strand.

    The input DNA is assumed to be written in the reference genome's
    forward 5'-to-3' orientation.

    Parameters
    ----------
    sequence:
        Genomic DNA sequence.

    strand:
        "+" for a forward-strand gene.
        "-" for a reverse-strand gene.

    allow_n:
        Whether unknown bases represented by N are allowed.

    Returns
    -------
    str
        Transcribed RNA sequence.

    Raises
    ------
    ValueError
        If strand is not "+" or "-".
    """
    if strand not in {"+", "-"}:
        raise ValueError("Strand must be either '+' or '-'.")

    cleaned_sequence = validate_sequence(
        sequence,
        allow_n=allow_n,
    )

    if strand == "-":
        cleaned_sequence = reverse_complement_dna(
            cleaned_sequence,
            allow_n=allow_n,
        )

    return cleaned_sequence.replace("T", "U")


def compare_rna_sequences(
    reference_rna: str,
    mutated_rna: str,
    strand: str,
) -> RNAComparison:
    """
    Compare two RNA sequences containing exactly one difference.

    Parameters
    ----------
    reference_rna:
        RNA generated from reference DNA.

    mutated_rna:
        RNA generated from mutated DNA.

    strand:
        Strand used to produce both RNA sequences.

    Returns
    -------
    RNAComparison
        Details of the RNA nucleotide change.
    """
    if len(reference_rna) != len(mutated_rna):
        raise ValueError(
            "Reference and mutated RNA sequences must have the same length."
        )

    differences = [
        index
        for index, (reference_base, mutated_base) in enumerate(
            zip(reference_rna, mutated_rna)
        )
        if reference_base != mutated_base
    ]

    if not differences:
        raise ValueError("No RNA difference was detected.")

    if len(differences) > 1:
        raise ValueError(
            "RNA comparison requires exactly one nucleotide difference, "
            f"but {len(differences)} differences were found."
        )

    position = differences[0]

    return RNAComparison(
        reference_rna=reference_rna,
        mutated_rna=mutated_rna,
        strand=strand,
        position_0_based=position,
        position_1_based=position + 1,
        reference_rna_base=reference_rna[position],
        alternate_rna_base=mutated_rna[position],
        sequence_length=len(reference_rna),
    )


def analyze_rna_change(
    reference_sequence: str,
    mutated_sequence: str,
    strand: str,
    allow_n: bool = False,
) -> RNAComparison:
    """
    Validate DNA, detect one SNV, transcribe both sequences, and compare RNA.
    """
    detect_snv(
        reference_sequence=reference_sequence,
        mutated_sequence=mutated_sequence,
        allow_n=allow_n,
    )

    reference_rna = dna_to_rna(
        sequence=reference_sequence,
        strand=strand,
        allow_n=allow_n,
    )

    mutated_rna = dna_to_rna(
        sequence=mutated_sequence,
        strand=strand,
        allow_n=allow_n,
    )

    return compare_rna_sequences(
        reference_rna=reference_rna,
        mutated_rna=mutated_rna,
        strand=strand,
    )


def analyze_variant_rna(
    variant: VariantInput,
) -> RNAComparison:
    """
    Analyze RNA change using a validated VariantInput object.

    RNA analysis requires known gene-strand information.
    """
    if variant.strand is None:
        raise ValueError(
            "RNA analysis requires strand information. "
            "Use genomic mode or provide the gene strand."
        )

    return analyze_rna_change(
        reference_sequence=variant.reference_sequence,
        mutated_sequence=variant.mutated_sequence,
        strand=variant.strand.value,
        allow_n=False,
    )