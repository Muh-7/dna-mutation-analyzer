"""Gradio interface for the DNA Mutation Analyzer."""

from typing import Any

import gradio as gr

from src.ui.analysis_service import (
    run_analysis_from_ui,
)


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
    """Run analysis and convert errors to readable UI output."""
    try:
        return run_analysis_from_ui(
            reference_sequence=reference_sequence,
            mutated_sequence=mutated_sequence,
            analysis_mode=analysis_mode,
            chromosome=chromosome,
            genomic_position=genomic_position,
            gene_name=gene_name,
            tissue=tissue,
            strand=strand,
            flank_size=flank_size,
            run_expression=run_expression,
            run_splicing=run_splicing,
            run_tf_binding=run_tf_binding,
        )

    except Exception as error:
        error_result = {
            "status": "failed",
            "error_type": type(error).__name__,
            "message": str(error),
        }

        error_markdown = (
            "## Analysis failed\n\n"
            f"**{type(error).__name__}:** "
            f"{error}\n\n"
            "Please check the DNA sequences and the selected "
            "analysis options."
        )

        return (
            error_markdown,
            error_result,
            error_result,
            error_result,
            error_result,
            error_result,
            error_result,
            error_result,
        )


with gr.Blocks(
    title="DNA Mutation Analyzer",
    fill_width=True,
) as demo:
    gr.Markdown(
        """
# DNA Mutation Analyzer

Compare a normal DNA sequence with a mutated DNA sequence
and analyze the predicted effects on:

- RNA transcription
- Genomic context
- Gene expression
- RNA splicing
- Transcription-factor motif matching

> Educational and research use only.
"""
    )

    with gr.Row():
        reference_sequence = gr.Textbox(
            label="Normal / Reference DNA",
            placeholder=(
                "Enter the normal DNA sequence using A, C, G, T..."
            ),
            lines=8,
            max_lines=20,
        )

        mutated_sequence = gr.Textbox(
            label="Mutated DNA",
            placeholder=(
                "Enter the mutated DNA sequence using A, C, G, T..."
            ),
            lines=8,
            max_lines=20,
        )

    with gr.Row():
        analysis_mode = gr.Radio(
            choices=[
                (
                    "Sequence only",
                    "raw_sequence",
                ),
                (
                    "Genomic context",
                    "genomic",
                ),
            ],
            value="raw_sequence",
            label="Analysis mode",
        )

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

    with gr.Accordion(
        "Optional genomic information",
        open=False,
    ):
        gr.Markdown(
            """
These fields are required only when using
**Genomic context** mode.
"""
        )

        with gr.Row():
            chromosome = gr.Textbox(
                label="Chromosome",
                placeholder="Example: chr22",
            )

            genomic_position = gr.Number(
                label="Genomic position (1-based)",
                precision=0,
                placeholder="Example: 36201698",
            )

        with gr.Row():
            gene_name = gr.Textbox(
                label="Gene name",
                placeholder="Example: APOL4",
            )

            tissue = gr.Textbox(
                label="Tissue",
                placeholder="Example: colon",
            )

    with gr.Accordion(
        "Analysis options",
        open=True,
    ):
        gr.Markdown(
            """
Expression and splicing require **Genomic context** mode.
TF motif analysis works directly from the two DNA sequences.
"""
        )

        with gr.Row():
            run_tf_binding = gr.Checkbox(
                value=True,
                label="Run JASPAR/FIMO TF motif analysis",
            )

            run_expression = gr.Checkbox(
                value=False,
                label="Run AlphaGenome expression analysis",
            )

            run_splicing = gr.Checkbox(
                value=False,
                label="Run AlphaGenome splicing analysis",
            )

    with gr.Row():
        analyze_button = gr.Button(
            "Run analysis",
            variant="primary",
        )

        clear_button = gr.ClearButton(
            value="Clear",
        )

    status_output = gr.Markdown(
        value=(
            "Enter the two DNA sequences, choose the "
            "analysis options, and press **Run analysis**."
        )
    )

    with gr.Tabs():
        with gr.Tab("Mutation"):
            mutation_output = gr.JSON(
                label="Detected mutation",
                open=True,
            )

        with gr.Tab("Genomic annotation"):
            genomic_output = gr.JSON(
                label="GENCODE annotation",
                open=True,
            )

        with gr.Tab("RNA"):
            rna_output = gr.JSON(
                label="RNA change",
                open=True,
            )

        with gr.Tab("Expression"):
            expression_output = gr.JSON(
                label="AlphaGenome expression",
                open=False,
            )

        with gr.Tab("Splicing"):
            splicing_output = gr.JSON(
                label="AlphaGenome splicing",
                open=False,
            )

        with gr.Tab("TF motifs"):
            tf_binding_output = gr.JSON(
                label="JASPAR/FIMO analysis",
                open=False,
            )

        with gr.Tab("Complete result"):
            complete_output = gr.JSON(
                label="Complete pipeline output",
                open=False,
                max_height=700,
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
        run_expression,
        run_splicing,
        run_tf_binding,
    ]

    analysis_outputs = [
        status_output,
        mutation_output,
        genomic_output,
        rna_output,
        expression_output,
        splicing_output,
        tf_binding_output,
        complete_output,
    ]

    analyze_button.click(
        fn=handle_analysis,
        inputs=analysis_inputs,
        outputs=analysis_outputs,
        concurrency_limit=1,
    )

    clear_button.add(
        components=[
            reference_sequence,
            mutated_sequence,
            chromosome,
            genomic_position,
            gene_name,
            tissue,
            mutation_output,
            genomic_output,
            rna_output,
            expression_output,
            splicing_output,
            tf_binding_output,
            complete_output,
        ]
    )


if __name__ == "__main__":
    demo.queue().launch(
        inbrowser=True,
        show_error=True,
    )