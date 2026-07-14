"""Main analysis pipeline for the DNA Mutation Analyzer project."""

from dataclasses import asdict
from pathlib import Path
from typing import Any

from src.annotation.genomic_annotator import annotate_genomic_position
from src.preprocessing.pipeline import prepare_variant_model
from src.rna.transcription import analyze_variant_rna
from src.schemas.variant_input import AnalysisMode, VariantInput


def run_analysis_pipeline(
    variant: VariantInput,
    annotation_database_path: str | Path | None = None,
) -> dict[str, Any]:
    """
    Run all currently supported DNA mutation analyses.

    The pipeline currently performs:

    1. DNA validation.
    2. SNV detection.
    3. Mutation-window extraction.
    4. Genomic metadata organization.
    5. Direct DNA-to-RNA transcription.
    6. RNA mutation comparison.
    7. Genomic region annotation when a database is provided.

    Parameters
    ----------
    variant:
        Validated DNA variant input.

    annotation_database_path:
        Path to the local GENCODE SQLite database.

    Returns
    -------
    dict
        Complete structured mutation-analysis result.
    """
    preprocessing_result = prepare_variant_model(variant)

    result: dict[str, Any] = {
        "project_version": "0.3.0",
        "analysis_status": "completed",
        "preprocessing": preprocessing_result,
        "rna_analysis": {
            "status": "not_run",
            "reason": (
                "RNA analysis requires the gene strand to be provided."
            ),
        },
        "genomic_annotation": {
            "status": "not_run",
            "reason": (
                "Genomic annotation requires genomic mode and "
                "an annotation database."
            ),
        },
        "scientific_notes": [
            (
                "The RNA output represents direct transcription from "
                "the provided genomic DNA sequence."
            ),
            (
                "The RNA output is not necessarily mature mRNA because "
                "transcript selection, exon joining, and intron removal "
                "have not yet been applied."
            ),
            (
                "Genomic-region labels are derived from the selected "
                "GENCODE annotation database."
            ),
            (
                "The outputs are computational results intended for "
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

    if (
        variant.analysis_mode == AnalysisMode.GENOMIC
        and annotation_database_path is not None
    ):
        annotation_result = annotate_genomic_position(
            database_path=annotation_database_path,
            chromosome=variant.chromosome,
            position_1_based=variant.genomic_position,
            gene_name=variant.gene_name,
        )

        result["genomic_annotation"] = {
            "status": "completed",
            **asdict(annotation_result),
        }

    return result