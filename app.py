"""Modern Gradio interface for the DNA Mutation Analyzer."""

from html import escape
from pathlib import Path
from typing import Any

import gradio as gr
import pandas as pd
from Bio import SeqIO

from src.ui.analysis_service import run_analysis_from_ui


APP_THEME = gr.themes.Soft(
    primary_hue="emerald",
    secondary_hue="teal",
    neutral_hue="slate",
)

APP_CSS = """
.gradio-container {
    min-height: 100vh;
    background:
        radial-gradient(circle at 10% 0%, rgba(16,185,129,.10), transparent 30rem),
        radial-gradient(circle at 95% 10%, rgba(6,182,212,.08), transparent 30rem);
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
    background: linear-gradient(135deg, #064e3b 0%, #0f766e 55%, #0e7490 100%);
    box-shadow: 0 22px 50px rgba(6,78,59,.20);
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
    opacity: .92;
}

.hero-tags {
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
    margin-top: 22px;
}

.hero-tag {
    padding: 7px 12px;
    border: 1px solid rgba(255,255,255,.25);
    border-radius: 999px;
    background: rgba(255,255,255,.11);
    font-size: .82rem;
    font-weight: 650;
}

.input-card {
    padding: 8px;
    border-radius: 20px !important;
    border: 1px solid var(--border-color-primary) !important;
    box-shadow: 0 12px 32px rgba(15,23,42,.05);
}

#run-button button {
    min-height: 52px;
    font-size: 1rem;
    font-weight: 750;
    box-shadow: 0 12px 26px rgba(5,150,105,.22);
}

.result-banner {
    display: flex;
    align-items: center;
    gap: 14px;
    margin: 22px 0 16px;
    padding: 17px 19px;
    border: 1px solid rgba(16,185,129,.30);
    border-radius: 18px;
    background: rgba(16,185,129,.08);
}

.result-banner.error {
    border-color: rgba(239,68,68,.35);
    background: rgba(239,68,68,.08);
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
    opacity: .76;
}

.metrics-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
    gap: 14px;
    margin-bottom: 18px;
}

.metric-card {
    position: relative;
    overflow: hidden;
    min-height: 130px;
    padding: 18px;
    border: 1px solid var(--border-color-primary);
    border-radius: 18px;
    background: var(--block-background-fill);
    box-shadow: 0 8px 23px rgba(15,23,42,.05);
}

.metric-card::before {
    content: "";
    position: absolute;
    inset: 0 0 auto;
    height: 4px;
    background: #94a3b8;
}

.metric-card.mutation::before { background: #059669; }
.metric-card.genomic::before { background: #0891b2; }
.metric-card.expression::before { background: #2563eb; }
.metric-card.splicing::before { background: #7c3aed; }
.metric-card.tf::before { background: #db2777; }

.metric-title {
    margin-bottom: 12px;
    color: var(--body-text-color-subdued);
    font-size: .78rem;
    font-weight: 750;
    letter-spacing: .07em;
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
    font-size: .86rem;
    line-height: 1.4;
}

.sequence-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
    gap: 16px;
}

.sequence-box {
    padding: 18px;
    border: 1px solid var(--border-color-primary);
    border-radius: 17px;
    background: var(--block-background-fill);
}

.sequence-title {
    margin-bottom: 12px;
    font-weight: 750;
}

.sequence-code {
    overflow-wrap: anywhere;
    font-family: "JetBrains Mono", "Fira Code", monospace;
    font-size: 1rem;
    line-height: 1.9;
    letter-spacing: .06em;
}

.reference-base,
.mutated-base {
    display: inline-block;
    padding: 1px 6px;
    border-radius: 6px;
    color: white;
    font-weight: 850;
}

.reference-base { background: #059669; }
.mutated-base { background: #e11d48; }

footer { display: none !important; }

@media (max-width: 720px) {
    #main-container { padding: 10px 8px 40px; }
    .hero { padding: 25px 20px; border-radius: 20px; }
    .hero h1 { font-size: 2.2rem; }
}
"""

INITIAL_SUMMARY = """
<section>
    <div class="result-banner">
        <span class="result-icon">→</span>
        <div>
            <strong>Ready for analysis</strong>
            <p>Enter the two DNA sequences, select the analyses, and press Run analysis.</p>
        </div>
    </div>
</section>
"""

INITIAL_FINDINGS = """
## Key findings

The scientific interpretation will appear here after running the analysis.
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


def get_section(
    result: dict[str, Any],
    name: str,
) -> dict[str, Any]:
    """Return one pipeline section safely."""
    section = result.get(name)
    return section if isinstance(section, dict) else {}


def is_completed(section: dict[str, Any]) -> bool:
    """Check whether an analysis section completed."""
    return section.get("status") == "completed"


def format_number(value: Any, digits: int = 4) -> str:
    """Format a numeric value for display."""
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
    """Create one summary card."""
    return f"""
    <div class="metric-card {escape(css_class)}">
        <div class="metric-title">{escape(title)}</div>
        <div class="metric-value">{escape(value)}</div>
        <div class="metric-subtitle">{escape(subtitle)}</div>
    </div>
    """


def build_summary_html(result: dict[str, Any]) -> str:
    """Create the result overview cards."""
    preprocessing = get_section(result, "preprocessing")
    mutation = preprocessing.get("mutation", {})

    reference_base = mutation.get("reference_base", "?")
    alternate_base = mutation.get("alternate_base", "?")
    position = mutation.get("position_1_based", "?")

    mutation_card = create_metric_card(
        "DNA mutation",
        f"{reference_base} → {alternate_base}",
        f"Sequence position {position}",
        "mutation",
    )

    annotation = get_section(result, "genomic_annotation")
    if is_completed(annotation):
        genomic_card = create_metric_card(
            "Genomic context",
            str(annotation.get("gene_name", "Unknown")),
            str(annotation.get("primary_region", "Unknown")).replace("_", " "),
            "genomic",
        )
    else:
        genomic_card = create_metric_card(
            "Genomic context",
            "Sequence only",
            "Genomic annotation was not run",
            "genomic",
        )

    expression = get_section(result, "expression_analysis")
    if is_completed(expression):
        summary = expression.get("summary", {})
        direction = str(summary.get("direction", "Unknown")).replace("_", " ").title()
        expression_card = create_metric_card(
            "Gene expression",
            direction,
            f"Raw score: {format_number(summary.get('raw_score'))}",
            "expression",
        )
    else:
        expression_card = create_metric_card(
            "Gene expression",
            "Not run",
            "Optional AlphaGenome analysis",
            "expression",
        )

    splicing = get_section(result, "splicing_analysis")
    if is_completed(splicing):
        summary = (
            splicing.get("target_gene_summary")
            or splicing.get("global_summary")
            or {}
        )
        splicing_card = create_metric_card(
            "Splicing score",
            format_number(summary.get("merged_splicing_score")),
            "Merged AlphaGenome score",
            "splicing",
        )
    else:
        splicing_card = create_metric_card(
            "Splicing score",
            "Not run",
            "Optional AlphaGenome analysis",
            "splicing",
        )

    tf_binding = get_section(result, "tf_binding_analysis")
    if is_completed(tf_binding):
        counts = tf_binding.get("counts", {})
        tf_card = create_metric_card(
            "TF motif changes",
            f"{counts.get('gained', 0)} gained",
            f"{counts.get('lost', 0)} lost motif matches",
            "tf",
        )
    else:
        tf_card = create_metric_card(
            "TF motif changes",
            "Not run",
            "Optional JASPAR/FIMO analysis",
            "tf",
        )

    return f"""
    <section>
        <div class="result-banner">
            <span class="result-icon">✓</span>
            <div>
                <strong>Analysis completed successfully</strong>
                <p>The most important predicted effects are summarized below.</p>
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


def highlight_nucleotide(
    sequence: str,
    index: int | None,
    css_class: str,
) -> str:
    """Highlight one nucleotide in a DNA sequence."""
    if not sequence:
        return "Sequence not available."

    if index is None or index < 0 or index >= len(sequence):
        return escape(sequence)

    return (
        escape(sequence[:index])
        + f'<span class="{escape(css_class)}">'
        + escape(sequence[index])
        + "</span>"
        + escape(sequence[index + 1:])
    )


def build_sequence_html(result: dict[str, Any]) -> str:
    """Show the reference and mutated sequences."""
    preprocessing = get_section(result, "preprocessing")
    window = preprocessing.get("window", {})

    reference_sequence = (
        window.get("reference_window")
        or preprocessing.get("reference_sequence")
        or ""
    )
    mutated_sequence = (
        window.get("mutated_window")
        or preprocessing.get("mutated_sequence")
        or ""
    )

    mutation_index = window.get("mutation_index_in_window")
    if mutation_index is None:
        mutation_index = preprocessing.get("mutation", {}).get("position_0_based")

    reference_html = highlight_nucleotide(
        reference_sequence,
        mutation_index,
        "reference-base",
    )
    mutated_html = highlight_nucleotide(
        mutated_sequence,
        mutation_index,
        "mutated-base",
    )

    return f"""
    <div class="sequence-grid">
        <div class="sequence-box">
            <div class="sequence-title">Normal / reference DNA</div>
            <div class="sequence-code">{reference_html}</div>
        </div>
        <div class="sequence-box">
            <div class="sequence-title">Mutated DNA</div>
            <div class="sequence-code">{mutated_html}</div>
        </div>
    </div>
    """


def build_findings_markdown(result: dict[str, Any]) -> str:
    """Create a short scientific interpretation."""
    findings: list[str] = []

    preprocessing = get_section(result, "preprocessing")
    mutation = preprocessing.get("mutation", {})
    findings.append(
        "A single-nucleotide variant "
        f"`{mutation.get('reference_base', '?')}>{mutation.get('alternate_base', '?')}` "
        f"was detected at sequence position `{mutation.get('position_1_based', '?')}`."
    )

    annotation = get_section(result, "genomic_annotation")
    if is_completed(annotation):
        findings.append(
            f"The variant is associated with gene `{annotation.get('gene_name', 'Unknown')}` "
            f"and is located in a predicted "
            f"`{annotation.get('primary_region', 'Unknown')}` region."
        )

    expression = get_section(result, "expression_analysis")
    if is_completed(expression):
        summary = expression.get("summary", {})
        direction = str(summary.get("direction", "Unknown")).replace("_", " ")
        findings.append(
            f"AlphaGenome predicted a `{direction}` in expression of "
            f"`{expression.get('gene_name', 'the target gene')}`."
        )

    splicing = get_section(result, "splicing_analysis")
    if is_completed(splicing):
        summary = (
            splicing.get("target_gene_summary")
            or splicing.get("global_summary")
            or {}
        )
        findings.append(
            "The predicted merged splicing score is "
            f"`{format_number(summary.get('merged_splicing_score'))}`."
        )

    tf_binding = get_section(result, "tf_binding_analysis")
    if is_completed(tf_binding):
        counts = tf_binding.get("counts", {})
        findings.append(
            "FIMO detected "
            f"`{counts.get('gained', 0)}` gained, "
            f"`{counts.get('lost', 0)}` lost, "
            f"`{counts.get('strengthened', 0)}` strengthened, and "
            f"`{counts.get('weakened', 0)}` weakened motif matches "
            "overlapping the mutation."
        )

    lines = ["## Key findings", ""]
    lines.extend(f"- {finding}" for finding in findings)

    warnings = result.get("warnings", [])
    if warnings:
        lines.extend(["", "### Warnings"])
        for warning in warnings:
            if isinstance(warning, dict):
                lines.append(f"- {warning.get('message', str(warning))}")
            else:
                lines.append(f"- {warning}")

    lines.extend(
        [
            "",
            "> These outputs are computational predictions for research and "
            "educational use. They are not a medical diagnosis.",
        ]
    )
    return "\n".join(lines)


def build_genomic_markdown(result: dict[str, Any]) -> str:
    """Display genomic annotation as a compact table."""
    annotation = get_section(result, "genomic_annotation")
    if not is_completed(annotation):
        reason = annotation.get(
            "reason",
            "Genomic annotation was not requested.",
        )
        return f"### Genomic context\n\n{reason}"

    regions = annotation.get("region_labels", ())
    region_text = ", ".join(
        str(region).replace("_", " ")
        for region in regions
    )
    transcripts = annotation.get("transcript_ids", ())

    return (
        "### Genomic context\n\n"
        "| Field | Result |\n"
        "|---|---|\n"
        f"| Gene | `{annotation.get('gene_name', 'Unknown')}` |\n"
        f"| Gene ID | `{annotation.get('gene_id', 'Unknown')}` |\n"
        f"| Gene type | `{annotation.get('gene_type', 'Unknown')}` |\n"
        f"| Chromosome | `{annotation.get('chromosome', 'Unknown')}` |\n"
        f"| Genomic position | `{annotation.get('position_1_based', 'Unknown')}` |\n"
        f"| Strand | `{annotation.get('strand', 'Unknown')}` |\n"
        f"| Primary region | "
        f"`{str(annotation.get('primary_region', 'Unknown')).replace('_', ' ')}` |\n"
        f"| Region labels | `{region_text or 'Unknown'}` |\n"
        f"| Distance to TSS | `{annotation.get('distance_to_tss', 'Unknown')}` |\n"
        f"| Matched transcripts | `{len(transcripts)}` |"
    )


def shorten_sequence(
    sequence: str | None,
    limit: int = 180,
) -> str:
    """Shorten a long RNA sequence for display."""
    if not sequence:
        return "Not available"
    if len(sequence) <= limit:
        return sequence

    half = limit // 2
    return sequence[:half] + "..." + sequence[-half:]


def build_rna_markdown(result: dict[str, Any]) -> str:
    """Display RNA transcription results."""
    rna = get_section(result, "rna_analysis")
    if not is_completed(rna):
        reason = rna.get(
            "reason",
            "RNA analysis was not run because the strand is unknown.",
        )
        return f"### RNA transcription\n\n{reason}"

    reference_rna = shorten_sequence(rna.get("reference_rna"))
    mutated_rna = shorten_sequence(rna.get("mutated_rna"))

    return (
        "### RNA transcription\n\n"
        "| Field | Result |\n"
        "|---|---|\n"
        f"| RNA change | "
        f"`{rna.get('reference_rna_base', '?')} → "
        f"{rna.get('alternate_rna_base', '?')}` |\n"
        f"| RNA position | `{rna.get('position_1_based', 'Unknown')}` |\n"
        f"| Gene strand | `{rna.get('strand', 'Unknown')}` |\n"
        f"| Strand source | `{rna.get('strand_source', 'Unknown')}` |\n"
        "| RNA type | `Direct genomic transcription` |\n\n"
        "**Reference RNA**\n\n"
        f"```text\n{reference_rna}\n```\n\n"
        "**Mutated RNA**\n\n"
        f"```text\n{mutated_rna}\n```\n\n"
        "> This result represents direct genomic transcription. "
        "It is not necessarily mature processed mRNA."
    )


def build_expression_dataframe(result: dict[str, Any]) -> pd.DataFrame:
    """Create the gene-expression table."""
    expression = get_section(result, "expression_analysis")
    if not is_completed(expression):
        return pd.DataFrame(columns=EXPRESSION_COLUMNS)

    rows = [
        {
            "Biosample": track.get("biosample_name"),
            "GTEx tissue": track.get("gtex_tissue"),
            "Track": track.get("track_name"),
            "Raw score": track.get("raw_score"),
            "Quantile score": track.get("quantile_score"),
        }
        for track in expression.get("tracks", [])
    ]
    return pd.DataFrame(rows, columns=EXPRESSION_COLUMNS)


def build_splicing_dataframe(result: dict[str, Any]) -> pd.DataFrame:
    """Create the RNA-splicing table."""
    splicing = get_section(result, "splicing_analysis")
    if not is_completed(splicing):
        return pd.DataFrame(columns=SPLICING_COLUMNS)

    rows = [
        {
            "Output": str(record.get("output_type", "")).replace("_", " ").title(),
            "Biosample": record.get("biosample_name"),
            "GTEx tissue": record.get("gtex_tissue"),
            "Track": record.get("track_name"),
            "Raw score": record.get("raw_score"),
            "Quantile score": record.get("quantile_score"),
        }
        for record in splicing.get("strongest_target_gene_records", [])
    ]
    return pd.DataFrame(rows, columns=SPLICING_COLUMNS)


def build_tf_dataframe(result: dict[str, Any]) -> pd.DataFrame:
    """Create the transcription-factor motif table."""
    tf_binding = get_section(result, "tf_binding_analysis")
    if not is_completed(tf_binding):
        return pd.DataFrame(columns=TF_COLUMNS)

    rows = [
        {
            "Change": str(change.get("change_type", "")).title(),
            "Factor": change.get("motif_alt_id"),
            "Motif": change.get("motif_id"),
            "Start": change.get("start"),
            "Stop": change.get("stop"),
            "Strand": change.get("strand"),
            "Reference score": change.get("reference_score"),
            "Mutated score": change.get("mutated_score"),
            "Score delta": change.get("score_delta"),
            "Reference p-value": change.get("reference_p_value"),
            "Mutated p-value": change.get("mutated_p_value"),
        }
        for change in tf_binding.get("changes", [])
    ]
    return pd.DataFrame(rows, columns=TF_COLUMNS)


def handle_analysis(
    reference_sequence: str,
    mutated_sequence: str,
    analysis_mode: str,
    chromosome: str | None,
    genomic_position: int | float | None,
    gene_name: str | None,
    tissue: str | None,
    strand: str | None,
    flank_size: int | float,
    selected_analyses: list[str] | None,
) -> tuple[Any, ...]:
    """Run the analysis and prepare UI-friendly outputs."""
    selected = set(selected_analyses or [])

    try:
        outputs = run_analysis_from_ui(
            reference_sequence=reference_sequence,
            mutated_sequence=mutated_sequence,
            analysis_mode=analysis_mode,
            chromosome=chromosome,
            genomic_position=genomic_position,
            gene_name=gene_name,
            tissue=tissue,
            strand=strand,
            flank_size=flank_size,
            run_expression="expression" in selected,
            run_splicing="splicing" in selected,
            run_tf_binding="tf_binding" in selected,
        )

        complete_result = outputs[-1]

        return (
            build_summary_html(complete_result),
            build_findings_markdown(complete_result),
            build_sequence_html(complete_result),
            build_genomic_markdown(complete_result),
            build_rna_markdown(complete_result),
            build_expression_dataframe(complete_result),
            build_splicing_dataframe(complete_result),
            build_tf_dataframe(complete_result),
            complete_result,
        )

    except Exception as error:
        error_result = {
            "status": "failed",
            "error_type": type(error).__name__,
            "message": str(error),
        }

        error_html = f"""
        <section>
            <div class="result-banner error">
                <span class="result-icon">!</span>
                <div>
                    <strong>Analysis failed</strong>
                    <p>{escape(type(error).__name__)}: {escape(str(error))}</p>
                </div>
            </div>
        </section>
        """

        error_markdown = (
            "## Unable to complete the analysis\n\n"
            f"**{type(error).__name__}:** {error}\n\n"
            "Check the DNA sequences, genomic fields, analysis mode, "
            "and selected analyses."
        )

        return (
            error_html,
            error_markdown,
            INITIAL_SEQUENCES,
            "### Genomic context\n\nNo result.",
            "### RNA transcription\n\nNo result.",
            pd.DataFrame(columns=EXPRESSION_COLUMNS),
            pd.DataFrame(columns=SPLICING_COLUMNS),
            pd.DataFrame(columns=TF_COLUMNS),
            error_result,
        )


def load_apol4_example() -> tuple[Any, ...]:
    """Load the verified APOL4 example."""
    fasta_path = Path("data/tf_binding/apol4_reference_mutated.fa")

    if fasta_path.exists():
        sequences = {
            record.id: str(record.seq)
            for record in SeqIO.parse(fasta_path, "fasta")
        }
        reference_sequence = sequences.get("reference", "GACTCACCCGA")
        mutated_sequence = sequences.get("mutated", "GACTCCCCCGA")
        flank_size = 50
        analyses = ["tf_binding", "expression", "splicing"]
    else:
        reference_sequence = "GACTCACCCGA"
        mutated_sequence = "GACTCCCCCGA"
        flank_size = 5
        analyses = ["expression", "splicing"]

    return (
        reference_sequence,
        mutated_sequence,
        "genomic",
        "chr22",
        36201698,
        "APOL4",
        "colon",
        "-",
        flank_size,
        analyses,
    )


def clear_interface() -> tuple[Any, ...]:
    """Reset all inputs and outputs."""
    return (
        "",
        "",
        "raw_sequence",
        "",
        None,
        "",
        "",
        "",
        100,
        ["tf_binding"],
        INITIAL_SUMMARY,
        INITIAL_FINDINGS,
        INITIAL_SEQUENCES,
        "### Genomic context\n\nNo analysis yet.",
        "### RNA transcription\n\nNo analysis yet.",
        pd.DataFrame(columns=EXPRESSION_COLUMNS),
        pd.DataFrame(columns=SPLICING_COLUMNS),
        pd.DataFrame(columns=TF_COLUMNS),
        {},
    )


with gr.Blocks(
    title="DNA Mutation Analyzer",
    fill_width=True,
) as demo:
    with gr.Column(elem_id="main-container"):
        gr.HTML(
            """
            <section class="hero">
                <h1>DNA Mutation Analyzer</h1>
                <p>
                    Compare normal and mutated DNA and study predicted effects
                    on RNA transcription, genomic context, gene expression,
                    splicing, and transcription-factor motif matching.
                </p>
                <div class="hero-tags">
                    <span class="hero-tag">AlphaGenome</span>
                    <span class="hero-tag">GENCODE</span>
                    <span class="hero-tag">GRCh38</span>
                    <span class="hero-tag">JASPAR</span>
                    <span class="hero-tag">FIMO</span>
                </div>
            </section>
            """
        )

        with gr.Row(equal_height=False):
            with gr.Column(scale=3, variant="panel", elem_classes="input-card"):
                gr.Markdown("## 1. Enter the DNA sequences")
                gr.Markdown(
                    "The sequences must have equal lengths and contain "
                    "exactly one nucleotide difference."
                )

                reference_sequence = gr.Textbox(
                    label="Normal / reference DNA",
                    placeholder="Example: GACTCACCCGA",
                    lines=7,
                    max_lines=20,
                    buttons=["copy"],
                )

                mutated_sequence = gr.Textbox(
                    label="Mutated DNA",
                    placeholder="Example: GACTCCCCCGA",
                    lines=7,
                    max_lines=20,
                    buttons=["copy"],
                )

            with gr.Column(scale=2, variant="panel", elem_classes="input-card"):
                gr.Markdown("## 2. Select the analysis")

                analysis_mode = gr.Radio(
                    choices=[
                        ("Sequence only", "raw_sequence"),
                        ("Genomic context", "genomic"),
                    ],
                    value="raw_sequence",
                    label="Analysis mode",
                    info="Expression and splicing require Genomic context.",
                )

                selected_analyses = gr.CheckboxGroup(
                    choices=[
                        ("TF motif changes", "tf_binding"),
                        ("Gene expression", "expression"),
                        ("RNA splicing", "splicing"),
                    ],
                    value=["tf_binding"],
                    label="Optional analyses",
                )

                with gr.Accordion("Genomic information", open=False):
                    chromosome = gr.Textbox(
                        label="Chromosome",
                        placeholder="Example: chr22",
                    )
                    genomic_position = gr.Number(
                        label="Genomic position",
                        precision=0,
                        placeholder="Example: 36201698",
                    )
                    gene_name = gr.Textbox(
                        label="Gene name",
                        placeholder="Example: APOL4",
                    )
                    tissue = gr.Textbox(
                        label="Tissue",
                        placeholder="Example: colon",
                    )

                with gr.Accordion("Advanced settings", open=False):
                    strand = gr.Dropdown(
                        choices=[
                            ("Unknown", ""),
                            ("Forward strand (+)", "+"),
                            ("Reverse strand (-)", "-"),
                        ],
                        value="",
                        label="Gene strand",
                    )
                    flank_size = gr.Slider(
                        minimum=1,
                        maximum=1000,
                        value=100,
                        step=1,
                        label="Mutation-window flank size",
                    )

                with gr.Row():
                    example_button = gr.Button(
                        "Load APOL4 example",
                        variant="secondary",
                    )
                    clear_button = gr.Button(
                        "Clear",
                        variant="secondary",
                    )

                analyze_button = gr.Button(
                    "Run analysis",
                    variant="primary",
                    elem_id="run-button",
                )

        summary_output = gr.HTML(value=INITIAL_SUMMARY)
        findings_output = gr.Markdown(value=INITIAL_FINDINGS)

        with gr.Tabs():
            with gr.Tab("DNA mutation"):
                gr.Markdown(
                    "### DNA sequence comparison\n\n"
                    "The normal nucleotide is green and the mutated "
                    "nucleotide is red."
                )
                sequences_output = gr.HTML(value=INITIAL_SEQUENCES)

            with gr.Tab("Genomic context and RNA"):
                with gr.Row(equal_height=False):
                    genomic_output = gr.Markdown(
                        value="### Genomic context\n\nNo analysis yet."
                    )
                    rna_output = gr.Markdown(
                        value="### RNA transcription\n\nNo analysis yet."
                    )

            with gr.Tab("Gene expression"):
                gr.Markdown(
                    "### AlphaGenome expression predictions\n\n"
                    "A negative raw score suggests a predicted decrease; "
                    "a positive score suggests a predicted increase."
                )
                expression_output = gr.Dataframe(
                    headers=EXPRESSION_COLUMNS,
                    datatype="auto",
                    interactive=False,
                    wrap=True,
                    show_search="filter",
                    buttons=["copy", "fullscreen"],
                    max_height=470,
                )

            with gr.Tab("RNA splicing"):
                gr.Markdown(
                    "### Strongest AlphaGenome splicing results\n\n"
                    "The table shows the strongest splice-site, usage, "
                    "and junction records."
                )
                splicing_output = gr.Dataframe(
                    headers=SPLICING_COLUMNS,
                    datatype="auto",
                    interactive=False,
                    wrap=True,
                    show_search="filter",
                    buttons=["copy", "fullscreen"],
                    max_height=470,
                )

            with gr.Tab("TF motifs"):
                gr.Markdown(
                    "### JASPAR/FIMO motif changes\n\n"
                    "Gained and lost results depend on the selected "
                    "FIMO statistical threshold."
                )
                tf_output = gr.Dataframe(
                    headers=TF_COLUMNS,
                    datatype="auto",
                    interactive=False,
                    wrap=True,
                    show_search="filter",
                    buttons=["copy", "fullscreen"],
                    max_height=520,
                )

        with gr.Accordion("Advanced technical result", open=False):
            complete_output = gr.JSON(
                label="Complete pipeline output",
                open=False,
                show_indices=False,
                max_height=650,
                buttons=["copy"],
            )

        gr.Markdown(
            "> **Research disclaimer:** This application provides "
            "computational predictions for educational and research "
            "purposes. It is not a medical diagnostic system."
        )

        analysis_inputs = [
            reference_sequence,
            mutated_sequence,
            analysis_mode,
            chromosome,
            genomic_position,
            gene_name,
            tissue,
            strand,
            flank_size,
            selected_analyses,
        ]

        analysis_outputs = [
            summary_output,
            findings_output,
            sequences_output,
            genomic_output,
            rna_output,
            expression_output,
            splicing_output,
            tf_output,
            complete_output,
        ]

        analyze_button.click(
            fn=handle_analysis,
            inputs=analysis_inputs,
            outputs=analysis_outputs,
            concurrency_limit=1,
            scroll_to_output=True,
            show_progress="full",
        )

        example_button.click(
            fn=load_apol4_example,
            inputs=[],
            outputs=analysis_inputs,
            queue=False,
        )

        clear_button.click(
            fn=clear_interface,
            inputs=[],
            outputs=[*analysis_inputs, *analysis_outputs],
            queue=False,
        )


if __name__ == "__main__":
    demo.queue(default_concurrency_limit=1).launch(
        inbrowser=True,
        show_error=True,
        theme=APP_THEME,
        css=APP_CSS,
    )
