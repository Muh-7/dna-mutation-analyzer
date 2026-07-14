"""Tests for the Gradio analysis service."""

from pathlib import Path

import pytest

from src.ui.analysis_service import (
    build_variant_input_from_ui,
    run_analysis_from_ui,
)


def test_build_raw_sequence_variant() -> None:
    variant = build_variant_input_from_ui(
        reference_sequence="AACCGAACTG",
        mutated_sequence="AACCGCACTG",
        analysis_mode="raw_sequence",
        chromosome=None,
        genomic_position=None,
        gene_name=None,
        tissue=None,
        strand="+",
        flank_size=2,
    )

    assert (
        variant.analysis_mode.value
        == "raw_sequence"
    )

    assert variant.reference_sequence == (
        "AACCGAACTG"
    )

    assert variant.mutated_sequence == (
        "AACCGCACTG"
    )

    assert variant.flank_size == 2
    assert variant.strand.value == "+"


def test_build_genomic_variant() -> None:
    variant = build_variant_input_from_ui(
        reference_sequence="GACTCACCCGA",
        mutated_sequence="GACTCCCCCGA",
        analysis_mode="genomic",
        chromosome="22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="-",
        flank_size=5,
    )

    assert variant.analysis_mode.value == "genomic"
    assert variant.chromosome == "chr22"
    assert variant.genomic_position == 36201698
    assert variant.gene_name == "APOL4"
    assert variant.tissue == "colon"


def test_fractional_position_is_rejected() -> None:
    with pytest.raises(
        ValueError,
        match="must be an integer",
    ):
        build_variant_input_from_ui(
            reference_sequence="GACTCACCCGA",
            mutated_sequence="GACTCCCCCGA",
            analysis_mode="genomic",
            chromosome="22",
            genomic_position=36201698.5,
            gene_name="APOL4",
            tissue="colon",
            strand="-",
            flank_size=5,
        )


def test_alphagenome_requires_genomic_mode() -> None:
    with pytest.raises(
        ValueError,
        match="require Genomic Context mode",
    ):
        run_analysis_from_ui(
            reference_sequence="AACCGAACTG",
            mutated_sequence="AACCGCACTG",
            analysis_mode="raw_sequence",
            chromosome=None,
            genomic_position=None,
            gene_name=None,
            tissue=None,
            strand="+",
            flank_size=2,
            run_expression=True,
            run_splicing=False,
            run_tf_binding=False,
        )


def test_ui_service_returns_sections(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake_pipeline_result = {
        "project_version": "0.7.0",
        "preprocessing": {
            "mutation": {
                "mutation_type": "SNV",
                "position_0_based": 5,
                "position_1_based": 6,
                "reference_base": "A",
                "alternate_base": "C",
                "sequence_length": 10,
            }
        },
        "genomic_annotation": {
            "status": "not_run",
        },
        "rna_analysis": {
            "status": "completed",
            "strand": "+",
        },
        "expression_analysis": {
            "status": "not_run",
        },
        "splicing_analysis": {
            "status": "not_run",
        },
        "tf_binding_analysis": {
            "status": "not_run",
        },
        "warnings": [],
    }

    monkeypatch.setattr(
        (
            "src.ui.analysis_service."
            "run_analysis_pipeline"
        ),
        lambda **kwargs: fake_pipeline_result,
    )

    outputs = run_analysis_from_ui(
        reference_sequence="AACCGAACTG",
        mutated_sequence="AACCGCACTG",
        analysis_mode="raw_sequence",
        chromosome=None,
        genomic_position=None,
        gene_name=None,
        tissue=None,
        strand="+",
        flank_size=2,
        run_expression=False,
        run_splicing=False,
        run_tf_binding=False,
    )

    assert len(outputs) == 8
    assert "Analysis completed" in outputs[0]
    assert outputs[1]["position_1_based"] == 6
    assert outputs[3]["status"] == "completed"