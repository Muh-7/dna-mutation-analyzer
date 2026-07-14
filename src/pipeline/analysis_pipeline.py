"""Main analysis pipeline for the DNA Mutation Analyzer project."""

from src.annotation.genomic_annotator import annotate_genomic_position
from src.preprocessing.pipeline import prepare_variant_model
from src.rna.transcription import analyze_rna_change
from src.schemas.variant_input import AnalysisMode, VariantInput


from dataclasses import asdict
from pathlib import Path
from typing import Any

from src.models.alphagenome_client import (
    analyze_variant_expression,
    create_alphagenome_client,
)



from src.tf_binding.analysis import (
    analyze_tf_binding_sequences,
)





def run_analysis_pipeline(
    variant: VariantInput,
    annotation_database_path: str | Path | None = None,
    run_expression_analysis: bool = False,
    alphagenome_model: Any | None = None,
    run_tf_binding_analysis: bool = False,
    motif_database_path: str | Path | None = None,
    fimo_executable: str = "fimo",
    tf_binding_p_value_threshold: float = 1e-4,
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
        "project_version": "0.5.0",
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
        
        "expression_analysis": {
    "status": "not_run",
    "reason": (
        "AlphaGenome expression analysis was not requested."
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
            (
                "AlphaGenome expression scores compare predicted ALT and "
                "REF RNA-seq signals and are not direct laboratory measurements."
            ),
            (
                "FIMO results describe statistical sequence-to-motif "
                "matches and do not directly prove transcription-factor "
                "binding in living cells."
            ),            
        ],
        
        "tf_binding_analysis": {
                "status": "not_run",
                "reason": (
                "FIMO TF-binding analysis was not requested."
        ),
        },
        
        
        "warnings": [],
        
    }
    
    
    annotation_strand: str | None = None

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

        annotation_data = asdict(annotation_result)

        result["genomic_annotation"] = {
            "status": "completed",
            **annotation_data,
        }

        annotation_strand = annotation_result.strand    
    
        
    user_strand = (
        variant.strand.value
        if variant.strand is not None
        else None
    )

    effective_strand = annotation_strand or user_strand

    if (
        annotation_strand is not None
        and user_strand is not None
        and annotation_strand != user_strand
    ):
        result["warnings"].append(
            {
                "type": "strand_mismatch",
                "message": (
                    f"The input strand was {user_strand}, but GENCODE "
                    f"annotates {variant.gene_name} on strand "
                    f"{annotation_strand}. The GENCODE strand was used "
                    "for RNA transcription."
                ),
            }
        )

    if effective_strand is not None:
        rna_result = analyze_rna_change(
            reference_sequence=variant.reference_sequence,
            mutated_sequence=variant.mutated_sequence,
            strand=effective_strand,
            allow_n=False,
        )

        result["rna_analysis"] = {
            "status": "completed",
            "analysis_type": "direct_genomic_transcription",
            "is_mature_mrna": False,
            "strand_source": (
                "GENCODE"
                if annotation_strand is not None
                else "user_input"
            ),
            **asdict(rna_result),
        }
    
    
    

    
    if run_expression_analysis:
        if variant.analysis_mode != AnalysisMode.GENOMIC:
            result["expression_analysis"] = {
                "status": "not_run",
                "reason": (
                    "AlphaGenome expression analysis requires "
                    "genomic mode."
                ),
            }

        else:
            model = (
                alphagenome_model
                if alphagenome_model is not None
                else create_alphagenome_client()
            )

            result["expression_analysis"] = (
                analyze_variant_expression(
                    model=model,
                    variant_input=variant,
                )
            )
    
    
    if run_tf_binding_analysis:
        if motif_database_path is None:
            raise ValueError(
                "motif_database_path is required when "
                "run_tf_binding_analysis=True."
            )

        preprocessing_data = result["preprocessing"]
        window_data = preprocessing_data["window"]

        result["tf_binding_analysis"] = (
            analyze_tf_binding_sequences(
                reference_sequence=(
                    window_data["reference_window"]
                ),
                mutated_sequence=(
                    window_data["mutated_window"]
                ),
                motif_database_path=(
                    motif_database_path
                ),
                fimo_executable=fimo_executable,
                p_value_threshold=(
                    tf_binding_p_value_threshold
                ),
            )
        )    
        
        
    
    return result