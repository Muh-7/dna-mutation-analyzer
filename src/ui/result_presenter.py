"""Prepare scientific pipeline results for the Gradio interface."""

from html import escape
from typing import Any

import pandas as pd


def _section(
    result: dict[str, Any],
    section_name: str,
) -> dict[str, Any]:
    """Return a pipeline section safely."""
    section = result.get(section_name)

    if isinstance(section, dict):
        return section

    return {}


def _is_completed(section: dict[str, Any]) -> bool:
    """Check whether an analysis section completed."""
    return section.get("status") == "completed"


def _format_number(
    value: Any,
    digits: int = 4,
) -> str:
    """Format a numeric result for human-readable output."""
    if value is None:
        return "Not available"

    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return str(value)

    return f"{numeric_value:.{digits}f}"


def _truncate_sequence(
    sequence: str | None,
    maximum_length: int = 180,
) -> str:
    """Shorten very long sequences for Markdown display."""
    if not sequence:
        return "Not available"

    if len(sequence) <= maximum_length:
        return sequence

    half = maximum_length // 2

    return (
        sequence[:half]
        + "..."
        + sequence[-half:]
    )


def _metric_card(
    title: str,
    value: str,
    subtitle: str,
    accent_class: str = "",
) -> str:
    """Create one HTML summary card."""
    classes = "metric-card"

    if accent_class:
        classes += f" {accent_class}"

    return f"""
    <div class="{classes}">
        <div class="metric-title">{escape(title)}</div>
        <div class="metric-value">{escape(value)}</div>
        <div class="metric-subtitle">{escape(subtitle)}</div>
    </div>
    """


def build_overview_html(
    result: dict[str, Any],
) -> str:
    """Build the visual summary cards shown above the result tabs."""
    preprocessing = _section(
        result,
        "preprocessing",
    )

    mutation = preprocessing.get(
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

    mutation_position = mutation.get(
        "position_1_based",
        "?",
    )

    mutation_card = _metric_card(
        title="DNA mutation",
        value=f"{reference_base} → {alternate_base}",
        subtitle=(
            f"Sequence position {mutation_position}"
        ),
        accent_class="accent-primary",
    )

    annotation = _section(
        result,
        "genomic_annotation",
    )

    if _is_completed(annotation):
        gene_name = annotation.get(
            "gene_name",
            "Unknown",
        )

        primary_region = annotation.get(
            "primary_region",
            "Unknown region",
        )

        genomic_card = _metric_card(
            title="Genomic context",
            value=str(gene_name),
            subtitle=str(primary_region).replace(
                "_",
                " ",
            ),
            accent_class="accent-genomic",
        )
    else:
        genomic_card = _metric_card(
            title="Genomic context",
            value="Sequence only",
            subtitle="No genomic annotation",
        )

    expression = _section(
        result,
        "expression_analysis",
    )

    if _is_completed(expression):
        expression_summary = expression.get(
            "summary",
            {},
        )

        direction = expression_summary.get(
            "direction",
            "Unknown",
        )

        readable_direction = str(
            direction
        ).replace(
            "_",
            " ",
        ).title()

        raw_score = _format_number(
            expression_summary.get(
                "raw_score"
            )
        )

        expression_card = _metric_card(
            title="Gene expression",
            value=readable_direction,
            subtitle=f"Raw score: {raw_score}",
            accent_class=(
                "accent-decrease"
                if direction == "predicted_decrease"
                else "accent-increase"
            ),
        )
    else:
        expression_card = _metric_card(
            title="Gene expression",
            value="Not run",
            subtitle="Optional AlphaGenome analysis",
        )

    splicing = _section(
        result,
        "splicing_analysis",
    )

    if _is_completed(splicing):
        splicing_summary = (
            splicing.get(
                "target_gene_summary"
            )
            or splicing.get(
                "global_summary"
            )
            or {}
        )

        merged_score = _format_number(
            splicing_summary.get(
                "merged_splicing_score"
            )
        )

        splicing_card = _metric_card(
            title="Splicing effect",
            value=merged_score,
            subtitle="Merged AlphaGenome score",
            accent_class="accent-splicing",
        )
    else:
        splicing_card = _metric_card(
            title="Splicing effect",
            value="Not run",
            subtitle="Optional AlphaGenome analysis",
        )

    tf_binding = _section(
        result,
        "tf_binding_analysis",
    )

    if _is_completed(tf_binding):
        counts = tf_binding.get(
            "counts",
            {},
        )

        gained = counts.get(
            "gained",
            0,
        )

        lost = counts.get(
            "lost",
            0,
        )

        tf_card = _metric_card(
            title="TF motif changes",
            value=f"{gained} gained",
            subtitle=f"{lost} lost motif matches",
            accent_class="accent-tf",
        )
    else:
        tf_card = _metric_card(
            title="TF motif changes",
            value="Not run",
            subtitle="Optional JASPAR/FIMO analysis",
        )

    return f"""
    <section class="result-overview">
        <div class="success-banner">
            <span class="success-icon">✓</span>
            <div>
                <strong>Analysis completed successfully</strong>
                <p>
                    The most important predicted effects are
                    summarized below.
                </p>
            </div>
        </div>

        <div class="metrics-grid">
            {mutation_card}
            {genomic_card}
            {expression_card}
            {splicing_card}
            {tf_card}
        </div>
    </section>
    """


def _highlight_sequence(
    sequence: str,
    mutation_index: int | None,
    label: str,
) -> list[tuple[str, str | None]]:
    """Highlight one mutated nucleotide inside a sequence."""
    if not sequence:
        return []

    if (
        mutation_index is None
        or mutation_index < 0
        or mutation_index >= len(sequence)
    ):
        return [(sequence, None)]

    parts: list[tuple[str, str | None]] = []

    prefix = sequence[:mutation_index]
    nucleotide = sequence[mutation_index]
    suffix = sequence[mutation_index + 1:]

    if prefix:
        parts.append(
            (prefix, None)
        )

    parts.append(
        (nucleotide, label)
    )

    if suffix:
        parts.append(
            (suffix, None)
        )

    return parts


def build_sequence_highlights(
    result: dict[str, Any],
) -> tuple[
    list[tuple[str, str | None]],
    list[tuple[str, str | None]],
]:
    """Build reference and mutated sequence visualizations."""
    preprocessing = _section(
        result,
        "preprocessing",
    )

    window = preprocessing.get(
        "window",
        {},
    )

    reference_sequence = (
        window.get("reference_window")
        or preprocessing.get(
            "reference_sequence"
        )
        or ""
    )

    mutated_sequence = (
        window.get("mutated_window")
        or preprocessing.get(
            "mutated_sequence"
        )
        or ""
    )

    mutation_index = window.get(
        "mutation_index_in_window"
    )

    if mutation_index is None:
        mutation = preprocessing.get(
            "mutation",
            {},
        )

        mutation_index = mutation.get(
            "position_0_based"
        )

    return (
        _highlight_sequence(
            sequence=reference_sequence,
            mutation_index=mutation_index,
            label="Reference base",
        ),
        _highlight_sequence(
            sequence=mutated_sequence,
            mutation_index=mutation_index,
            label="Mutated base",
        ),
    )


def build_key_findings_markdown(
    result: dict[str, Any],
) -> str:
    """Create the short scientific interpretation section."""
    findings: list[str] = []

    final_report = result.get(
        "final_report",
        {},
    )

    report_findings = final_report.get(
        "key_findings",
        [],
    )

    if isinstance(report_findings, list):
        findings.extend(
            str(finding)
            for finding in report_findings
            if finding
        )

    if not findings:
        preprocessing = _section(
            result,
            "preprocessing",
        )

        mutation = preprocessing.get(
            "mutation",
            {},
        )

        findings.append(
            "A single-nucleotide variant "
            f"`{mutation.get('reference_base', '?')}"
            f">{mutation.get('alternate_base', '?')}` "
            "was detected at sequence position "
            f"`{mutation.get('position_1_based', '?')}`."
        )

        annotation = _section(
            result,
            "genomic_annotation",
        )

        if _is_completed(annotation):
            findings.append(
                "The variant is associated with "
                f"`{annotation.get('gene_name', 'Unknown')}` "
                "and is annotated as "
                f"`{annotation.get('primary_region', 'Unknown')}`."
            )

        expression = _section(
            result,
            "expression_analysis",
        )

        if _is_completed(expression):
            summary = expression.get(
                "summary",
                {},
            )

            findings.append(
                "AlphaGenome predicted "
                f"`{str(summary.get('direction', 'unknown')).replace('_', ' ')}` "
                "for the target gene."
            )

        splicing = _section(
            result,
            "splicing_analysis",
        )

        if _is_completed(splicing):
            summary = (
                splicing.get(
                    "target_gene_summary"
                )
                or splicing.get(
                    "global_summary"
                )
                or {}
            )

            findings.append(
                "The merged predicted splicing score is "
                f"`{_format_number(summary.get('merged_splicing_score'))}`."
            )

        tf_binding = _section(
            result,
            "tf_binding_analysis",
        )

        if _is_completed(tf_binding):
            counts = tf_binding.get(
                "counts",
                {},
            )

            findings.append(
                "FIMO detected "
                f"`{counts.get('gained', 0)}` gained and "
                f"`{counts.get('lost', 0)}` lost motif matches "
                "overlapping the mutation."
            )

    markdown_lines = [
        "## Key findings",
        "",
    ]

    markdown_lines.extend(
        f"- {finding}"
        for finding in findings
    )

    warnings = result.get(
        "warnings",
        [],
    )

    if warnings:
        markdown_lines.extend(
            [
                "",
                "### Warnings",
            ]
        )

        for warning in warnings:
            if isinstance(warning, dict):
                message = warning.get(
                    "message",
                    str(warning),
                )
            else:
                message = str(warning)

            markdown_lines.append(
                f"- {message}"
            )

    markdown_lines.extend(
        [
            "",
            "> Results are computational predictions for research "
            "and educational use, not a medical diagnosis.",
        ]
    )

    return "\n".join(
        markdown_lines
    )


def build_genomic_markdown(
    result: dict[str, Any],
) -> str:
    """Create a compact genomic-context table."""
    annotation = _section(
        result,
        "genomic_annotation",
    )

    if not _is_completed(annotation):
        reason = annotation.get(
            "reason",
            "Genomic annotation was not requested.",
        )

        return (
            "### Genomic context\n\n"
            f"{reason}"
        )

    transcript_ids = annotation.get(
        "transcript_ids",
        (),
    )

    transcript_count = len(
        transcript_ids
    )

    region_labels = annotation.get(
        "region_labels",
        (),
    )

    readable_regions = ", ".join(
        str(region).replace(
            "_",
            " ",
        )
        for region in region_labels
    )

    return f"""
### Genomic context

| Field | Result |
|---|---|
| Gene | `{annotation.get("gene_name", "Unknown")}` |
| Gene ID | `{annotation.get("gene_id", "Unknown")}` |
| Gene type | `{annotation.get("gene_type", "Unknown")}` |
| Chromosome | `{annotation.get("chromosome", "Unknown")}` |
| Position | `{annotation.get("position_1_based", "Unknown")}` |
| Strand | `{annotation.get("strand", "Unknown")}` |
| Primary region | `{str(annotation.get("primary_region", "Unknown")).replace("_", " ")}` |
| Region labels | `{readable_regions or "Unknown"}` |
| Distance to TSS | `{annotation.get("distance_to_tss", "Unknown")}` |
| Matched transcripts | `{transcript_count}` |
"""


def build_rna_markdown(
    result: dict[str, Any],
) -> str:
    """Create a concise RNA result."""
    rna = _section(
        result,
        "rna_analysis",
    )

    if not _is_completed(rna):
        reason = rna.get(
            "reason",
            "RNA analysis was not run because the strand is unknown.",
        )

        return (
            "### RNA transcription\n\n"
            f"{reason}"
        )

    reference_rna = _truncate_sequence(
        rna.get(
            "reference_rna"
        )
    )

    mutated_rna = _truncate_sequence(
        rna.get(
            "mutated_rna"
        )
    )

    return f"""
### RNA transcription

| Field | Result |
|---|---|
| RNA change | `{rna.get("reference_rna_base", "?")} → {rna.get("alternate_rna_base", "?")}` |
| RNA position | `{rna.get("position_1_based", "Unknown")}` |
| Gene strand | `{rna.get("strand", "Unknown")}` |
| Strand source | `{rna.get("strand_source", "Unknown")}` |
| Output type | Direct genomic transcription |

**Reference RNA**

```text
{reference_rna}
