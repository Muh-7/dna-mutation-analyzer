"""Complete preprocessing pipeline for one DNA SNV."""

from dataclasses import asdict
from typing import Any

from src.preprocessing.mutation_detector import detect_snv
from src.preprocessing.sequence_validator import validate_sequence
from src.preprocessing.window_extractor import extract_mutation_window

from src.schemas.variant_input import VariantInput

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
    
def prepare_variant_model(
    variant: VariantInput,
) -> dict[str, Any]:
    """
    Prepare a validated VariantInput object for downstream analysis.

    Parameters
    ----------
    variant:
        A validated DNA variant input model.

    Returns
    -------
    dict
        Preprocessed sequences, mutation information,
        extracted windows, and genomic metadata.
    """
    preprocessing_result = prepare_variant_input(
        reference_sequence=variant.reference_sequence,
        mutated_sequence=variant.mutated_sequence,
        flank_size=variant.flank_size,
        allow_n=False,
    )

    preprocessing_result["metadata"] = {
        "analysis_mode": variant.analysis_mode.value,
        "genome_build": variant.genome_build.value,
        "chromosome": variant.chromosome,
        "genomic_position": variant.genomic_position,
        "gene_name": variant.gene_name,
        "tissue": variant.tissue,
        "strand": (
            variant.strand.value
            if variant.strand is not None
            else None
        ),
    }

    return preprocessing_result