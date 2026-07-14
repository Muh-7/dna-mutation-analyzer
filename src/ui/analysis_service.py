"""Services used by the Gradio user interface."""

from pathlib import Path
from typing import Any

from src.pipeline.analysis_pipeline import (
    run_analysis_pipeline,
)
from src.schemas.variant_input import VariantInput


ANNOTATION_DATABASE_PATH = Path(
    "data/annotations/"
    "gencode.v50.basic.annotation.db"
)

MOTIF_DATABASE_PATH = Path(
    "data/motifs/"
    "JASPAR2026_CORE_vertebrates.meme"
)


def _optional_text(
    value: str | None,
) -> str | None:
    """Normalize an optional text input."""
    if value is None:
        return None

    cleaned = value.strip()

    return cleaned or None


def _optional_positive_integer(
    value: int | float | str | None,
    field_name: str,
) -> int | None:
    """Normalize an optional positive integer."""
    if value is None or value == "":
        return None

    try:
        numeric_value = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"{field_name} must be an integer."
        ) from error

    if not numeric_value.is_integer():
        raise ValueError(
            f"{field_name} must be an integer."
        )

    integer_value = int(numeric_value)

    if integer_value < 1:
        raise ValueError(
            f"{field_name} must be greater than or equal to 1."
        )

    return integer_value


def build_variant_input_from_ui(
    reference_sequence: str,
    mutated_sequence: str,
    analysis_mode: str,
    chromosome: str | None,
    genomic_position: int | float | str | None,
    gene_name: str | None,
    tissue: str | None,
    strand: str | None,
    flank_size: int | float,
) -> VariantInput:
    """Create and validate a VariantInput from UI values."""
    normalized_mode = analysis_mode.strip()

    if normalized_mode not in {
        "raw_sequence",
        "genomic",
    }:
        raise ValueError(
            "Analysis mode must be raw_sequence or genomic."
        )

    normalized_strand = _optional_text(strand)

    if normalized_strand not in {
        None,
        "+",
        "-",
    }:
        raise ValueError(
            "Strand must be +, -, or empty."
        )

    normalized_flank_size = (
        _optional_positive_integer(
            flank_size,
            "Flank size",
        )
    )

    if normalized_flank_size is None:
        normalized_flank_size = 100

    if normalized_mode == "raw_sequence":
        return VariantInput(
            analysis_mode="raw_sequence",
            reference_sequence=reference_sequence,
            mutated_sequence=mutated_sequence,
            strand=normalized_strand,
            flank_size=normalized_flank_size,
        )

    return VariantInput(
        analysis_mode="genomic",
        reference_sequence=reference_sequence,
        mutated_sequence=mutated_sequence,
        chromosome=_optional_text(chromosome),
        genomic_position=_optional_positive_integer(
            genomic_position,
            "Genomic position",
        ),
        gene_name=_optional_text(gene_name),
        tissue=_optional_text(tissue),
        strand=normalized_strand,
        flank_size=normalized_flank_size,
    )


def _result_section(
    result: dict[str, Any],
    section_name: str,
) -> dict[str, Any]:
    """Return a pipeline section in a UI-safe form."""
    section = result.get(section_name)

    if isinstance(section, dict):
        return section

    return {
        "status": "not_available",
        "reason": (
            f"The {section_name} section was not produced."
        ),
    }


def build_status_markdown(
    result: dict[str, Any],
) -> str:
    """Create a readable overview of the analysis."""
    mutation = result.get(
        "preprocessing",
        {},
    ).get(
        "mutation",
        {},
    )

    reference_base = mutation.get(
        "reference_base",
        "?",
    )

    alternate_base = mutation.get(
        "alternate_base",
        "?",
    )

    position = mutation.get(
        "position_1_based",
        "?",
    )

    lines = [
        "## Analysis completed",
        "",
        (
            f"**Detected mutation:** "
            f"`{reference_base}>{alternate_base}` "
            f"at sequence position `{position}`"
        ),
    ]

    annotation = _result_section(
        result,
        "genomic_annotation",
    )

    if annotation.get("status") == "completed":
        lines.extend(
            [
                "",
                (
                    f"**Gene:** "
                    f"`{annotation.get('gene_name', 'Unknown')}`"
                ),
                (
                    f"**Genomic region:** "
                    f"`{annotation.get('primary_region', 'Unknown')}`"
                ),
                (
                    f"**Gene strand:** "
                    f"`{annotation.get('strand', 'Unknown')}`"
                ),
            ]
        )

    expression = _result_section(
        result,
        "expression_analysis",
    )

    if expression.get("status") == "completed":
        summary = expression.get(
            "summary",
            {},
        )

        lines.extend(
            [
                "",
                (
                    "**Expression prediction:** "
                    f"`{summary.get('direction', 'Unknown')}`"
                ),
                (
                    "**Strongest expression score:** "
                    f"`{summary.get('raw_score', 'Unknown')}`"
                ),
            ]
        )

    splicing = _result_section(
        result,
        "splicing_analysis",
    )

    if splicing.get("status") == "completed":
        splicing_summary = (
            splicing.get("target_gene_summary")
            or splicing.get("global_summary")
            or {}
        )

        merged_score = splicing_summary.get(
            "merged_splicing_score"
        )

        lines.extend(
            [
                "",
                (
                    "**Merged splicing score:** "
                    f"`{merged_score}`"
                ),
            ]
        )

    tf_binding = _result_section(
        result,
        "tf_binding_analysis",
    )

    if tf_binding.get("status") == "completed":
        counts = tf_binding.get(
            "counts",
            {},
        )

        lines.extend(
            [
                "",
                "**TF motif changes:**",
                (
                    f"- Gained: "
                    f"`{counts.get('gained', 0)}`"
                ),
                (
                    f"- Lost: "
                    f"`{counts.get('lost', 0)}`"
                ),
                (
                    f"- Strengthened: "
                    f"`{counts.get('strengthened', 0)}`"
                ),
                (
                    f"- Weakened: "
                    f"`{counts.get('weakened', 0)}`"
                ),
            ]
        )

    warnings = result.get(
        "warnings",
        [],
    )

    if warnings:
        lines.extend(
            [
                "",
                "### Warnings",
            ]
        )

        for warning in warnings:
            message = warning.get(
                "message",
                str(warning),
            )

            lines.append(
                f"- {message}"
            )

    lines.extend(
        [
            "",
            "> These results are computational predictions "
            "and are not a medical diagnosis.",
        ]
    )

    return "\n".join(lines)


def run_analysis_from_ui(
    reference_sequence: str,
    mutated_sequence: str,
    analysis_mode: str,
    chromosome: str | None,
    genomic_position: int | float | str | None,
    gene_name: str | None,
    tissue: str | None,
    strand: str | None,
    flank_size: int | float,
    run_expression: bool,
    run_splicing: bool,
    run_tf_binding: bool,
) -> tuple[
    str,
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
]:
    """Run the complete project pipeline for the UI."""
    variant = build_variant_input_from_ui(
        reference_sequence=reference_sequence,
        mutated_sequence=mutated_sequence,
        analysis_mode=analysis_mode,
        chromosome=chromosome,
        genomic_position=genomic_position,
        gene_name=gene_name,
        tissue=tissue,
        strand=strand,
        flank_size=flank_size,
    )

    if (
        variant.analysis_mode.value
        == "raw_sequence"
        and (
            run_expression
            or run_splicing
        )
    ):
        raise ValueError(
            "AlphaGenome expression and splicing analyses "
            "require Genomic Context mode."
        )

    annotation_database_path: Path | None = None

    if (
        variant.analysis_mode.value == "genomic"
        and ANNOTATION_DATABASE_PATH.exists()
    ):
        annotation_database_path = (
            ANNOTATION_DATABASE_PATH
        )

    motif_database_path: Path | None = None

    if run_tf_binding:
        if not MOTIF_DATABASE_PATH.exists():
            raise FileNotFoundError(
                "The JASPAR motif database was not found at: "
                f"{MOTIF_DATABASE_PATH}"
            )

        motif_database_path = MOTIF_DATABASE_PATH

    result = run_analysis_pipeline(
        variant=variant,
        annotation_database_path=(
            annotation_database_path
        ),
        run_expression_analysis=run_expression,
        run_splicing_analysis=run_splicing,
        run_tf_binding_analysis=run_tf_binding,
        motif_database_path=motif_database_path,
    )

    mutation_section = result.get(
        "preprocessing",
        {},
    ).get(
        "mutation",
        {},
    )

    return (
        build_status_markdown(result),
        mutation_section,
        _result_section(
            result,
            "genomic_annotation",
        ),
        _result_section(
            result,
            "rna_analysis",
        ),
        _result_section(
            result,
            "expression_analysis",
        ),
        _result_section(
            result,
            "splicing_analysis",
        ),
        _result_section(
            result,
            "tf_binding_analysis",
        ),
        result,
    )