"""Main analysis pipeline for the DNA Mutation Analyzer project."""

from dataclasses import asdict
from typing import Any

from src.preprocessing.pipeline import prepare_variant_model
from src.rna.transcription import analyze_variant_rna
from src.schemas.variant_input import VariantInput


def run_analysis_pipeline(
    variant: VariantInput,
) -> dict[str, Any]:
    """
    Run the currently supported DNA mutation analyses.

    The pipeline currently performs:

    1. DNA sequence validation.
    2. SNV detection.
    3. Mutation-window extraction.
    4. Genomic metadata organization.
    5. Direct DNA-to-RNA transcription when strand is available.
    6. Reference-versus-mutated RNA comparison.

    Parameters
    ----------
    variant:
        A validated VariantInput object.

    Returns
    -------
    dict
        A structured result containing preprocessing information,
        RNA analysis, and scientific limitations.
    """
    preprocessing_result = prepare_variant_model(variant)

    result: dict[str, Any] = {
        "project_version": "0.2.0",
        "analysis_status": "completed",
        "preprocessing": preprocessing_result,
        "rna_analysis": {
            "status": "not_run",
            "reason": (
                "RNA analysis requires the gene strand to be provided."
            ),
        },
        "scientific_notes": [
            (
                "The RNA output currently represents direct transcription "
                "from the provided genomic DNA sequence."
            ),
            (
                "The generated RNA is not necessarily mature mRNA because "
                "exon selection, intron removal, and transcript annotation "
                "have not yet been applied."
            ),
            (
                "The results are computational outputs intended for "
                "research and educational use."
            ),
        ],
    }

    if variant.strand is not None:
        rna_result = analyze_variant_rna(variant)

        result["rna_analysis"] = {
            "status": "completed",
            "analysis_type": "direct_genomic_transcription",
            "is_mature_mrna": False,
            **asdict(rna_result),
        }

    return result