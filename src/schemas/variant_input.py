"""Input schemas for DNA mutation analysis."""

from enum import Enum
from typing import Self

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    field_validator,
    model_validator,
)

from src.preprocessing.mutation_detector import detect_snv
from src.preprocessing.sequence_validator import validate_sequence


class AnalysisMode(str, Enum):
    """Supported project analysis modes."""

    RAW_SEQUENCE = "raw_sequence"
    GENOMIC = "genomic"


class GenomeBuild(str, Enum):
    """Supported human genome references."""

    HG38 = "hg38"


class Strand(str, Enum):
    """DNA strand direction."""

    FORWARD = "+"
    REVERSE = "-"


class VariantInput(BaseModel):
    """
    Store and validate all inputs required for mutation analysis.

    Raw-sequence mode requires only two DNA sequences.

    Genomic mode additionally requires chromosome, genomic position,
    gene name, tissue, and strand.
    """

    model_config = ConfigDict(
        extra="forbid",
        frozen=True,
        str_strip_whitespace=True,
    )

    reference_sequence: str
    mutated_sequence: str

    analysis_mode: AnalysisMode = AnalysisMode.RAW_SEQUENCE
    genome_build: GenomeBuild = GenomeBuild.HG38

    chromosome: str | None = None
    genomic_position: int | None = Field(default=None, ge=1)
    gene_name: str | None = None
    tissue: str | None = None
    strand: Strand | None = None

    flank_size: int = Field(default=100, ge=0)

    @field_validator(
        "reference_sequence",
        "mutated_sequence",
        mode="before",
    )
    @classmethod
    def validate_dna_sequence(cls, value: object) -> str:
        """Clean and validate each DNA sequence."""
        return validate_sequence(value, allow_n=False)  # type: ignore[arg-type]

    @field_validator("chromosome")
    @classmethod
    def normalize_chromosome(cls, value: str | None) -> str | None:
        """Normalize and validate a human chromosome name."""
        if value is None:
            return None

        chromosome = value.strip()

        if not chromosome.lower().startswith("chr"):
            chromosome = f"chr{chromosome}"

        suffix = chromosome[3:].upper()

        valid_suffixes = {
            *(str(number) for number in range(1, 23)),
            "X",
            "Y",
            "M",
            "MT",
        }

        if suffix not in valid_suffixes:
            raise ValueError(
                "Chromosome must be one of chr1-chr22, chrX, "
                "chrY, chrM, or chrMT."
            )

        if suffix == "MT":
            suffix = "M"

        return f"chr{suffix}"

    @model_validator(mode="after")
    def validate_complete_input(self) -> Self:
        """Validate the sequence pair and genomic-mode requirements."""
        detect_snv(
            reference_sequence=self.reference_sequence,
            mutated_sequence=self.mutated_sequence,
            allow_n=False,
        )

        if self.analysis_mode == AnalysisMode.GENOMIC:
            required_genomic_fields = {
                "chromosome": self.chromosome,
                "genomic_position": self.genomic_position,
                "gene_name": self.gene_name,
                "tissue": self.tissue,
                "strand": self.strand,
            }

            missing_fields = [
                field_name
                for field_name, value in required_genomic_fields.items()
                if value is None or value == ""
            ]

            if missing_fields:
                missing_text = ", ".join(missing_fields)

                raise ValueError(
                    "Genomic analysis mode requires these fields: "
                    f"{missing_text}."
                )

        return self
