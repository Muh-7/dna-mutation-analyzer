"""Fetch and validate genomic sequences from a reference FASTA."""

from dataclasses import dataclass
from pathlib import Path

from pyfaidx import Fasta

from src.preprocessing.sequence_validator import validate_sequence


VALID_SNV_BASES = frozenset({"A", "C", "G", "T"})


@dataclass(frozen=True)
class GenomicVariantSequences:
    """
    Store reference and mutated sequences around one genomic SNV.

    All genomic positions are based on the forward reference genome.
    """

    genome_build: str
    fasta_path: str

    chromosome: str
    position_1_based: int

    reference_allele: str
    alternate_allele: str

    reference_sequence: str
    mutated_sequence: str

    window_start_1_based: int
    window_end_1_based: int
    mutation_index_0_based: int

    left_flank_length: int
    right_flank_length: int
    flank_size_requested: int


def normalize_chromosome(chromosome: str) -> str:
    """
    Normalize a chromosome name to GENCODE style.

    Examples
    --------
    22 -> chr22
    X  -> chrX
    MT -> chrM
    """
    if not isinstance(chromosome, str):
        raise TypeError("Chromosome must be a string.")

    normalized = chromosome.strip()

    if not normalized:
        raise ValueError("Chromosome cannot be empty.")

    if not normalized.lower().startswith("chr"):
        normalized = f"chr{normalized}"

    suffix = normalized[3:].upper()

    if suffix == "MT":
        suffix = "M"

    valid_suffixes = {
        *(str(number) for number in range(1, 23)),
        "X",
        "Y",
        "M",
    }

    if suffix not in valid_suffixes:
        raise ValueError(
            "Chromosome must be one of chr1-chr22, chrX, chrY, or chrM."
        )

    return f"chr{suffix}"


def open_reference_genome(
    fasta_path: str | Path,
) -> Fasta:
    """
    Open and index a reference FASTA file.

    pyfaidx creates a .fai index automatically when one is not present.
    """
    fasta = Path(fasta_path)

    if not fasta.exists():
        raise FileNotFoundError(
            f"Reference FASTA file was not found: {fasta}"
        )

    if fasta.suffix == ".gz":
        raise ValueError(
            "Use the uncompressed .fa file. Standard gzip FASTA files "
            "cannot be indexed reliably for random access."
        )

    return Fasta(
        str(fasta),
        as_raw=True,
        sequence_always_upper=True,
        strict_bounds=True,
    )


def _validate_position(
    reference: Fasta,
    chromosome: str,
    position_1_based: int,
) -> None:
    """Validate chromosome availability and genomic position."""
    if isinstance(position_1_based, bool) or not isinstance(
        position_1_based,
        int,
    ):
        raise TypeError("Genomic position must be an integer.")

    if position_1_based < 1:
        raise ValueError(
            "Genomic position must be greater than or equal to 1."
        )

    if chromosome not in reference.keys():
        raise ValueError(
            f"Chromosome {chromosome} was not found in the reference FASTA."
        )

    chromosome_length = len(reference[chromosome])

    if position_1_based > chromosome_length:
        raise ValueError(
            f"Position {position_1_based} exceeds the length of "
            f"{chromosome}: {chromosome_length}."
        )


def fetch_reference_base(
    reference: Fasta,
    chromosome: str,
    position_1_based: int,
) -> str:
    """
    Fetch one nucleotide from a 1-based genomic position.
    """
    normalized_chromosome = normalize_chromosome(chromosome)

    _validate_position(
        reference=reference,
        chromosome=normalized_chromosome,
        position_1_based=position_1_based,
    )

    position_0_based = position_1_based - 1

    base = str(
        reference[normalized_chromosome][
            position_0_based:position_0_based + 1
        ]
    ).upper()

    if len(base) != 1:
        raise RuntimeError(
            "Reference genome returned an unexpected nucleotide length."
        )

    return base


def fetch_reference_window(
    reference: Fasta,
    chromosome: str,
    position_1_based: int,
    flank_size: int = 100,
) -> tuple[str, int, int, int]:
    """
    Fetch a reference sequence window around a genomic position.

    Returns
    -------
    tuple
        Sequence, 1-based window start, 1-based window end,
        and 0-based mutation index inside the window.
    """
    if isinstance(flank_size, bool) or not isinstance(flank_size, int):
        raise TypeError("Flank size must be an integer.")

    if flank_size < 0:
        raise ValueError("Flank size cannot be negative.")

    normalized_chromosome = normalize_chromosome(chromosome)

    _validate_position(
        reference=reference,
        chromosome=normalized_chromosome,
        position_1_based=position_1_based,
    )

    chromosome_length = len(reference[normalized_chromosome])
    mutation_position_0_based = position_1_based - 1

    window_start_0_based = max(
        0,
        mutation_position_0_based - flank_size,
    )

    window_end_0_based_exclusive = min(
        chromosome_length,
        mutation_position_0_based + flank_size + 1,
    )

    sequence = str(
        reference[normalized_chromosome][
            window_start_0_based:window_end_0_based_exclusive
        ]
    ).upper()

    sequence = validate_sequence(
        sequence,
        allow_n=True,
    )

    window_start_1_based = window_start_0_based + 1
    window_end_1_based = window_end_0_based_exclusive

    mutation_index_0_based = (
        mutation_position_0_based - window_start_0_based
    )

    return (
        sequence,
        window_start_1_based,
        window_end_1_based,
        mutation_index_0_based,
    )


def validate_reference_allele(
    actual_reference_base: str,
    expected_reference_allele: str,
    chromosome: str,
    position_1_based: int,
) -> str:
    """
    Confirm that a supplied REF allele matches the reference genome.
    """
    actual = actual_reference_base.strip().upper()
    expected = expected_reference_allele.strip().upper()

    if len(expected) != 1 or expected not in VALID_SNV_BASES:
        raise ValueError(
            "The expected reference allele must be one of A, C, G, or T."
        )

    if actual == "N":
        raise ValueError(
            f"The reference genome contains N at "
            f"{chromosome}:{position_1_based}; the allele is unknown."
        )

    if actual not in VALID_SNV_BASES:
        raise ValueError(
            f"Unexpected reference base at "
            f"{chromosome}:{position_1_based}: {actual}."
        )

    if actual != expected:
        raise ValueError(
            "Reference allele mismatch: "
            f"{chromosome}:{position_1_based} contains {actual} "
            f"in the reference genome, but {expected} was provided."
        )

    return actual


def apply_snv(
    reference_sequence: str,
    mutation_index_0_based: int,
    alternate_allele: str,
) -> str:
    """
    Apply one SNV to an already extracted reference sequence.
    """
    sequence = validate_sequence(
        reference_sequence,
        allow_n=True,
    )

    if isinstance(mutation_index_0_based, bool) or not isinstance(
        mutation_index_0_based,
        int,
    ):
        raise TypeError("Mutation index must be an integer.")

    if not 0 <= mutation_index_0_based < len(sequence):
        raise ValueError(
            "Mutation index is outside the reference sequence."
        )

    alternate = alternate_allele.strip().upper()

    if len(alternate) != 1 or alternate not in VALID_SNV_BASES:
        raise ValueError(
            "The alternate allele must be one of A, C, G, or T."
        )

    current_base = sequence[mutation_index_0_based]

    if current_base == "N":
        raise ValueError(
            "An SNV cannot be applied at an unknown reference base N."
        )

    if current_base == alternate:
        raise ValueError(
            "The alternate allele must differ from the reference allele."
        )

    return (
        sequence[:mutation_index_0_based]
        + alternate
        + sequence[mutation_index_0_based + 1:]
    )


def build_variant_sequences(
    fasta_path: str | Path,
    chromosome: str,
    position_1_based: int,
    reference_allele: str,
    alternate_allele: str,
    flank_size: int = 100,
    genome_build: str = "GRCh38.p14",
) -> GenomicVariantSequences:
    """
    Fetch the true reference sequence and create its mutated version.
    """
    normalized_chromosome = normalize_chromosome(chromosome)

    reference = open_reference_genome(fasta_path)

    try:
        (
            reference_sequence,
            window_start_1_based,
            window_end_1_based,
            mutation_index_0_based,
        ) = fetch_reference_window(
            reference=reference,
            chromosome=normalized_chromosome,
            position_1_based=position_1_based,
            flank_size=flank_size,
        )

        actual_reference_base = reference_sequence[
            mutation_index_0_based
        ]

        validated_reference_allele = validate_reference_allele(
            actual_reference_base=actual_reference_base,
            expected_reference_allele=reference_allele,
            chromosome=normalized_chromosome,
            position_1_based=position_1_based,
        )

        mutated_sequence = apply_snv(
            reference_sequence=reference_sequence,
            mutation_index_0_based=mutation_index_0_based,
            alternate_allele=alternate_allele,
        )

    finally:
        reference.close()

    normalized_alternate = alternate_allele.strip().upper()

    return GenomicVariantSequences(
        genome_build=genome_build,
        fasta_path=str(Path(fasta_path)),
        chromosome=normalized_chromosome,
        position_1_based=position_1_based,
        reference_allele=validated_reference_allele,
        alternate_allele=normalized_alternate,
        reference_sequence=reference_sequence,
        mutated_sequence=mutated_sequence,
        window_start_1_based=window_start_1_based,
        window_end_1_based=window_end_1_based,
        mutation_index_0_based=mutation_index_0_based,
        left_flank_length=mutation_index_0_based,
        right_flank_length=(
            len(reference_sequence)
            - mutation_index_0_based
            - 1
        ),
        flank_size_requested=flank_size,
    )