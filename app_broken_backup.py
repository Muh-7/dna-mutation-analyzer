"""Modern Gradio interface for the DNA Mutation Analyzer."""

from html import escape
from pathlib import Path
from typing import Any

import gradio as gr
import pandas as pd
from Bio import SeqIO

from src.ui.analysis_service import run_analysis_from_ui


# =========================================================
# Theme and CSS
# =========================================================

APP_THEME = gr.themes.Soft(
    primary_hue="emerald",
    secondary_hue="teal",
    neutral_hue="slate",
)


APP_CSS = """
.gradio-container {
    min-height: 100vh;
    background:
        radial-gradient(
            circle at 10% 0%,
            rgba(16, 185, 129, 0.10),
            transparent 30rem
        ),
        radial-gradient(
            circle at 95% 10%,
            rgba(6, 182, 212, 0.08),
            transparent 30rem
        );
}

#main-container {
    max-width: 1380px;
    margin: 0 auto;
    padding: 20px 14px 60px;
}

.hero {
    margin-bottom: 22px;
    padding: 34px;
    border-radius: 25px;
    color: white;
    background:
        linear-gradient(
            135deg,
            #064e3b 0%,
            #0f766e 55%,
            #0e7490 100%
        );
    box-shadow:
        0 22px 50px rgba(6, 78, 59, 0.20);
}

.hero h1 {
    margin: 0;
    font-size: clamp(2rem, 5vw, 3.4rem);
    line-height: 1.05;
}

.hero p {
    max-width: 800px;
    margin: 16px 0 0;
    font-size: 1.02rem;
    line-height: 1.7;
    opacity: 0.92;
}

.hero-tags {
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
    margin-top: 22px;
}

.hero-tag {
    padding: 7px 12px;
    border: 1px solid rgba(255, 255, 255, 0.25);
    border-radius: 999px;
    background: rgba(255, 255, 255, 0.11);
    font-size: 0.82rem;
    font-weight: 650;
}

.input-card {
    padding: 8px;
    border-radius: 20px !important;
    border:
        1px solid
        var(--border-color-primary) !important;
    box-shadow:
        0 12px 32px rgba(15, 23, 42, 0.05);
}

#run-button button {
    min-height: 52px;
    font-size: 1rem;
    font-weight: 750;
    box-shadow:
        0 12px 26px rgba(5, 150, 105, 0.22);
}

.result-banner {
    display: flex;
    align-items: center;
    gap: 14px;
    margin: 22px 0 16px;
    padding: 17px 19px;
    border: 1px solid rgba(16, 185, 129, 0.30);
    border-radius: 18px;
    background: rgba(16, 185, 129, 0.08);
}

.result-banner.error {
    border-color: rgba(239, 68, 68, 0.35);
    background: rgba(239, 68, 68, 0.08);
}

.result-icon {
    display: grid;
    width: 40px;
    height: 40px;
    flex: 0 0 40px;
    place-items: center;
    border-radius: 50%;
    color: white;
    background: #059669;
    font-size: 1.2rem;
    font-weight: 800;
}

.result-banner.error .result-icon {
    background: #dc2626;
}

.result-banner p {
    margin: 4px 0 0;
    opacity: 0.76;
}

.metrics-grid {
    display: grid;
    grid-template-columns:
        repeat(auto-fit, minmax(190px, 1fr));
    gap: 14px;
    margin-bottom: 18px;
}

.metric-card {
    position: relative;
    overflow: hidden;
    min-height: 130px;
    padding: 18px;
    border:
        1px solid
        var(--border-color-primary);
    border-radius: 18px;
    background: var(--block-background-fill);
    box-shadow:
        0 8px 23px rgba(15, 23, 42, 0.05);
}

.metric-card::before {
    content: "";
    position: absolute;
    inset: 0 0 auto;
    height: 4px;
    background: #94a3b8;
}

.metric-card.mutation::before {
    background: #059669;
}

.metric-card.genomic::before {
    background: #0891b2;
}

.metric-card.expression::before {
    background: #2563eb;
}

.metric-card.splicing::before {
    background: #7c3aed;
}

.metric-card.tf::before {
    background: #db2777;
}

.metric-title {
    margin-bottom: 12px;
    color: var(--body-text-color-subdued);
    font-size: 0.78rem;
    font-weight: 750;
    letter-spacing: 0.07em;
    text-transform: uppercase;
}

.metric-value {
    margin-bottom: 8px;
    font-size: 1.35rem;
    font-weight: 800;
    line-height: 1.2;
}

.metric-subtitle {
    color: var(--body-text-color-subdued);
    font-size: 0.86rem;
    line-height: 1.4;
}

.sequence-grid {
    display: grid;
    grid-template-columns:
        repeat(auto-fit, minmax(300px, 1fr));
    gap: 16px;
}

.sequence-box {
    padding: 18px;
    border:
        1px solid
        var(--border-color-primary);
    border-radius: 17px;
    background: var(--block-background-fill);
}

.sequence-title {
    margin-bottom: 12px;
    font-weight: 750;
}

.sequence-code {
    overflow-wrap: anywhere;
    font-family:
        "JetBrains Mono",
        "Fira Code",
        monospace;
    font-size: 1rem;
    line-height: 1.9;
    letter-spacing: 0.06em;
}

.reference-base {
    display: inline-block;
    padding: 1px 6px;
    border-radius: 6px;
    color: white;
    background: #059669;
    font-weight: 850;
}

.mutated-base {
    display: inline-block;
    padding: 1px 6px;
    border-radius: 6px;
    color: white;
    background: #e11d48;
    font-weight: 850;
}

footer {
    display: none !important;
}

@media (max-width: 720px) {
    #main-container {
        padding: 10px 8px 40px;
    }

    .hero {
        padding: 25px 20px;
        border-radius: 20px;
    }

    .hero h1 {
        font-size: 2.2rem;
    }
}
"""


# =========================================================
# Initial UI values
# =========================================================

INITIAL_SUMMARY = """
<section>
    <div class="result-banner">
        <span class="result-icon">→</span>

        <div>
            <strong>Ready for analysis</strong>

            <p>
                Enter the two DNA sequences, select the analyses,
                and press Run analysis.
            </p>
        </div>
    </div>
</section>
"""


INITIAL_FINDINGS = """
## Key findings

The scientific interpretation will appear here after running
the analysis.
"""


INITIAL_SEQUENCES = """
<div class="sequence-grid">
    <div class="sequence-box">
        <div class="sequence-title">Reference DNA</div>
        <div>No sequence analyzed yet.</div>
    </div>

    <div class="sequence-box">
        <div class="sequence-title">Mutated DNA</div>
        <div>No sequence analyzed yet.</div>
    </div>
</div>
"""


EXPRESSION_COLUMNS = [
    "Biosample",
    "GTEx tissue",
    "Track",
    "Raw score",
    "Quantile score",
]


SPLICING_COLUMNS = [
    "Output",
    "Biosample",
    "GTEx tissue",
    "Track",
    "Raw score",
    "Quantile score",
]


TF_COLUMNS = [
    "Change",
    "Factor",
    "Motif",
    "Start",
    "Stop",
    "Strand",
    "Reference score",
    "Mutated score",
    "Score delta",
    "Reference p-value",
    "Mutated p-value",
]


# =========================================================
# General helper functions
# =========================================================

def get_section(
    result: dict[str, Any],
    name: str,
) -> dict[str, Any]:
    """Safely return one pipeline section."""
    section = result.get(name)

    if isinstance(section, dict):
        return section

    return {}


def is_completed(
    section: dict[str, Any],
) -> bool:
    """Check whether an analysis section completed."""
    return section.get("status") == "completed"


def format_number(
    value: Any,
    digits: int = 4,
) -> str:
    """Format numeric values for the interface."""
    if value is None:
        return "Not available"

    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def create_metric_card(
    title: str,
    value: str,
    subtitle: str,
    css_class: str,
) -> str:
    """Create one result summary card."""
    return f"""
    <div class="metric-card {css_class}">
        <div class="metric-title">
            {escape(title)}
        </div>

        <div class="metric-value">
            {escape(value)}
        </div>

        <div class="metric-subtitle">
            {escape(subtitle)}
        </div>
    </div>
    """


# =========================================================
# Result summary cards
# =========================================================

def build_summary_html(
    result: dict[str, Any],
) -> str:
    """Create the main visual result summary."""
    preprocessing = get_section(
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

    mutation_card = create_metric_card(
        title="DNA mutation",
        value=f"{reference_base} → {alternate_base}",
        subtitle=f"Sequence position {mutation_position}",
        css_class="mutation",
    )

    annotation = get_section(
        result,
        "genomic_annotation",
    )

    if is_completed(annotation):
        genomic_card = create_metric_card(
            title="Genomic context",
            value=str(
                annotation.get(
                    "gene_name",
                    "Unknown",
                )
            ),
            subtitle=str(
                annotation.get(
                    "primary_region",
                    "Unknown region",
                )
            ).replace(
                "_",
                " ",
            ),
            css_class="genomic",
        )
    else:
        genomic_card = create_metric_card(
            title="Genomic context",
            value="Sequence only",
            subtitle="Genomic annotation was not run",
            css_class="genomic",
        )

    expression = get_section(
        result,
        "expression_analysis",
    )

    if is_completed(expression):
        expression_summary = expression.get(
            "summary",
            {},
        )

        direction = str(
            expression_summary.get(
                "direction",
                "Unknown",
            )
        ).replace(
            "_",
            " ",
        ).title()

        expression_card = create_metric_card(
            title="Gene expression",
            value=direction,
            subtitle=(
                "Raw score: "
                + format_number(
                    expression_summary.get(
                        "raw_score"
                    )
                )
            ),
            css_class="expression",
        )
    else:
        expression_card = create_metric_card(
            title="Gene expression",
            value="Not run",
            subtitle="Optional AlphaGenome analysis",
            css_class="expression",
        )

    splicing = get_section(
        result,
        "splicing_analysis",
    )

    if is_completed(splicing):
        splicing_summary = (
            splicing.get(
                "target_gene_summary"
            )
            or splicing.get(
                "global_summary"
            )
            or {}
        )

        splicing_card = create_metric_card(
            title="Splicing score",
            value=format_number(
                splicing_summary.get(
                    "merged_splicing_score"
                )
            ),
            subtitle="Merged AlphaGenome score",
            css_class="splicing",
        )
    else:
        splicing_card = create_metric_card(
            title="Splicing score",
            value="Not run",
            subtitle="Optional AlphaGenome analysis",
            css_class="splicing",
        )

    tf_binding = get_section(
        result,
        "tf_binding_analysis",
    )

    if is_completed(tf_binding):
        counts = tf_binding.get(
            "counts",
            {},
        )

        tf_card = create_metric_card(
            title="TF motif changes",
            value=(
                f"{counts.get('gained', 0)} gained"
            ),
            subtitle=(
                f"{counts.get('lost', 0)} lost motif matches"
            ),
            css_class="tf",
        )
    else:
        tf_card = create_metric_card(
            title="TF motif changes",
            value="Not run",
            subtitle="Optional JASPAR/FIMO analysis",
            css_class="tf",
        )

    return f"""
    <section>
        <div class="result-banner">
            <span class="result-icon">✓</span>

            <div>
                <strong>
                    Analysis completed successfully
                </strong>

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


# =========================================================
# DNA sequence visualization
# =========================================================

def highlight_nucleotide(
    sequence: str,
    index: int | None,
    css_class: str,
) -> str:
    """Highlight one nucleotide in a DNA sequence."""
    if not sequence:
        return "Sequence not available."

    safe_sequence = escape(sequence)

    if (
        index is None
        or index < 0
        or index >= len(sequence)
    ):
        return safe_sequence

    prefix = escape(
        sequence[:index]
    )

    nucleotide = escape(
        sequence[index]
    )

    suffix = escape(
        sequence[index + 1:]
    )

    return (
        prefix
        + f'<span class="{css_class}">'
        + nucleotide
        + "</span>"
        + suffix
    )


def build_sequence_html(
    result: dict[str, Any],
) -> str:
    """Visualize the reference and mutated sequences."""
    preprocessing = get_section(
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
        mutation_index = preprocessing.get(
            "mutation",
            {},
        ).get(
            "position_0_based"
        )

    highlighted_reference = highlight_nucleotide(
        sequence=reference_sequence,
        index=mutation_index,
        css_class="reference-base",
    )

    highlighted_mutated = highlight_nucleotide(
        sequence=mutated_sequence,
        index=mutation_index,
        css_class="mutated-base",
    )

    return f"""
    <div class="sequence-grid">
        <div class="sequence-box">
            <div class="sequence-title">
                Normal / reference DNA
            </div>

            <div class="sequence-code">
                {highlighted_reference}
            </div>
        </div>

        <div class="sequence-box">
            <div class="sequence-title">
                Mutated DNA
            </div>

            <div class="sequence-code">
                {highlighted_mutated}
            </div>
        </div>
    </div>
    """


# =========================================================
# Key findings
# =========================================================

def build_findings_markdown(
    result: dict[str, Any],
) -> str:
    """Create a short scientific interpretation."""
    findings: list[str] = []

    preprocessing = get_section(
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

    annotation = get_section(
        result,
        "genomic_annotation",
    )

    if is_completed(annotation):
        findings.append(
            "The variant is associated with gene "
            f"`{annotation.get('gene_name', 'Unknown')}` "
            "and is located in a predicted "
            f"`{annotation.get('primary_region', 'Unknown')}` "
            "region."
        )

    expression = get_section(
        result,
        "expression_analysis",
    )

    if is_completed(expression):
        expression_summary = expression.get(
            "summary",
            {},
        )

        direction = str(
            expression_summary.get(
                "direction",
                "Unknown",
            )
        ).replace(
            "_",
            " ",
        )

        findings.append(
            "AlphaGenome predicted a "
            f"`{direction}` in expression of "
            f"`{expression.get('gene_name', 'the target gene')}`."
        )

    splicing = get_section(
        result,
        "splicing_analysis",
    )

    if is_completed(splicing):
        splicing_summary = (
            splicing.get(
                "target_gene_summary"
            )
            or splicing.get(
                "global_summary"
            )
            or {}
        )

        findings.append(
            "The predicted merged splicing score is "
            f"`{format_number(splicing_summary.get('merged_splicing_score'))}`."
        )

    tf_binding = get_section(
        result,
        "tf_binding_analysis",
    )

    if is_completed(tf_binding):
        counts = tf_binding.get(
            "counts",
            {},
        )

        findings.append(
            "FIMO detected "
            f"`{counts.get('gained', 0)}` gained, "
            f"`{counts.get('lost', 0)}` lost, "
            f"`{counts.get('strengthened', 0)}` strengthened, "
            "and "
            f"`{counts.get('weakened', 0)}` weakened "
            "motif matches overlapping the mutation."
        )

    lines = [
        "## Key findings",
        "",
    ]

    lines.extend(
        f"- {finding}"
        for finding in findings
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
            if isinstance(warning, dict):
                message = warning.get(
                    "message",
                    str(warning),
                )
            else:
                message = str(warning)

            lines.append(
                f"- {message}"
            )

    lines.extend(
        [
            "",
            "> These outputs are computational predictions "
            "for research and educational use. They are not "
            "a medical diagnosis.",
        ]
    )

    return "\n".join(lines)


# =========================================================
# Genomic and RNA formatting
# =========================================================

def build_genomic_markdown(
    result: dict[str, Any],
) -> str:
    """Display genomic annotation as a compact table."""
    annotation = get_section(
        result,
        "genomic_annotation",
    )

    if not is_completed(annotation):
        reason = annotation.get(
            "reason",
            "Genomic annotation was not requested.",
        )

        return (
            "### Genomic context\n\n"
            f"{reason}"
        )

    regions = annotation.get(
        "region_labels",
        (),
    )

    region_text = ", ".join(
        str(region).replace(
            "_",
            " ",
        )
        for region in regions
    )

    transcripts = annotation.get(
        "transcript_ids",
        (),
    )

    return f"""
### Genomic context

| Field | Result |
|---|---|
| Gene | `{annotation.get("gene_name", "Unknown")}` |
| Gene ID | `{annotation.get("gene_id", "Unknown")}` |
| Gene type | `{annotation.get("gene_type", "Unknown")}` |
| Chromosome | `{annotation.get("chromosome", "Unknown")}` |
| Genomic position | `{annotation.get("position_1_based", "Unknown")}` |
| Strand | `{annotation.get("strand", "Unknown")}` |
| Primary region | `{str(annotation.get("primary_region", "Unknown")).replace("_", " ")}` |
| Region labels | `{region_text or "Unknown"}` |
| Distance to TSS | `{annotation.get("distance_to_tss", "Unknown")}` |
| Matched transcripts | `{len(transcripts)}` |
"""


def shorten_sequence(
    sequence: str | None,
    limit: int = 180,
) -> str:
    """Shorten long RNA sequences for display."""
    if not sequence:
        return "Not available"

    if len(sequence) <= limit:
        return sequence

    half = limit // 2

    return (
        sequence[:half]
        + "..."
        + sequence[-half:]
    )


def build_rna_markdown(
    result: dict[str, Any],
) -> str:
    """Display RNA transcription results."""
    rna = get_section(
        result,
        "rna_analysis",
    )

    if not is_completed(rna):
        reason = rna.get(
            "reason",
            (
                "RNA analysis was not run because "
                "the strand is unknown."
            ),
        )

        return (
            "### RNA transcription\n\n"
            f"{reason}"
        )

    reference_rna = shorten_sequence(
        rna.get("reference_rna")
    )

    mutated_rna = shorten_sequence(
        rna.get("mutated_rna")
    )

    return (
        "### RNA transcription\n\n"
        "| Field | Result |\n"
        "|---|---|\n"
        f"| RNA change | "
        f"`{rna.get('reference_rna_base', '?')} → "
        f"{rna.get('alternate_rna_base', '?')}` |\n"
        f"| RNA position | "
        f"`{rna.get('position_1_based', 'Unknown')}` |\n"
        f"| Gene strand | "
        f"`{rna.get('strand', 'Unknown')}` |\n"
        f"| Strand source | "
        f"`{rna.get('strand_source', 'Unknown')}` |\n"
        "| RNA type | `Direct genomic transcription` |\n\n"
        "**Reference RNA**\n\n"
        f"```text\n{reference_rna}\n```\n\n"
        "**Mutated RNA**\n\n"
        f"```text\n{mutated_rna}\n```\n\n"
        "> This result represents direct genomic transcription. "
        "It is not necessarily mature processed mRNA."
    )


### RNA transcription

| Field | Result |
|---|---|
| RNA change | `{rna.get("reference_rna_base", "?")} → {rna.get("alternate_rna_base", "?")}` |
| RNA position | `{rna.get("position_1_based", "Unknown")}` |
| Gene strand | `{rna.get("strand", "Unknown")}` |
| Strand source | `{rna.get("strand_source", "Unknown")}` |
| RNA type | Direct genomic transcription |

**Reference RNA**

```text
{shorten_sequence(rna.get("reference_rna"))}