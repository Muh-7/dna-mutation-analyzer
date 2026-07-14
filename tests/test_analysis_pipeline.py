"""Tests for the main DNA mutation analysis pipeline."""

from src.pipeline.analysis_pipeline import run_analysis_pipeline
from src.schemas.variant_input import VariantInput


def test_genomic_pipeline_returns_complete_analysis() -> None:
    variant = VariantInput(
        analysis_mode="genomic",
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
        chromosome="22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="+",
        flank_size=2,
    )

    result = run_analysis_pipeline(variant)

    assert result["project_version"] == "0.2.0"
    assert result["analysis_status"] == "completed"

    preprocessing = result["preprocessing"]

    assert preprocessing["mutation"]["mutation_type"] == "SNV"
    assert preprocessing["mutation"]["position_1_based"] == 6
    assert preprocessing["mutation"]["reference_base"] == "T"
    assert preprocessing["mutation"]["alternate_base"] == "G"

    assert preprocessing["window"]["reference_window"] == "CGTAC"
    assert preprocessing["window"]["mutated_window"] == "CGGAC"

    assert preprocessing["metadata"]["chromosome"] == "chr22"
    assert preprocessing["metadata"]["gene_name"] == "APOL4"
    assert preprocessing["metadata"]["tissue"] == "colon"

    rna = result["rna_analysis"]

    assert rna["status"] == "completed"
    assert rna["analysis_type"] == "direct_genomic_transcription"
    assert rna["is_mature_mrna"] is False
    assert rna["reference_rna"] == "AACCGUACGU"
    assert rna["mutated_rna"] == "AACCGGACGU"
    assert rna["reference_rna_base"] == "U"
    assert rna["alternate_rna_base"] == "G"
    assert rna["position_1_based"] == 6


def test_raw_sequence_pipeline_skips_rna_without_strand() -> None:
    variant = VariantInput(
        reference_sequence="ACGT",
        mutated_sequence="AGGT",
        flank_size=1,
    )

    result = run_analysis_pipeline(variant)

    rna = result["rna_analysis"]

    assert rna["status"] == "not_run"
    assert "strand" in rna["reason"].lower()


def test_reverse_strand_pipeline() -> None:
    variant = VariantInput(
        analysis_mode="genomic",
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
        chromosome="22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="-",
        flank_size=2,
    )

    result = run_analysis_pipeline(variant)

    rna = result["rna_analysis"]

    assert rna["status"] == "completed"
    assert rna["strand"] == "-"
    assert rna["reference_rna"] == "ACGUACGGUU"
    assert rna["mutated_rna"] == "ACGUCCGGUU"
    assert rna["position_1_based"] == 5
    assert rna["reference_rna_base"] == "A"
    assert rna["alternate_rna_base"] == "C"


def test_pipeline_contains_scientific_notes() -> None:
    variant = VariantInput(
        reference_sequence="ACGT",
        mutated_sequence="AGGT",
    )

    result = run_analysis_pipeline(variant)

    notes = result["scientific_notes"]

    assert len(notes) >= 1
    assert any("mature mRNA" in note for note in notes)