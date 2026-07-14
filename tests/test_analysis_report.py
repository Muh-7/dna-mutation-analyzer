"""Tests for the unified mutation analysis report."""

import pytest

from src.reporting.analysis_report import (
    build_analysis_report,
)


def test_build_complete_analysis_report() -> None:
    pipeline_result = {
        "project_version": "0.7.0",
        "preprocessing": {
            "mutation": {
                "mutation_type": "SNV",
                "position_0_based": 50,
                "position_1_based": 51,
                "reference_base": "A",
                "alternate_base": "C",
                "sequence_length": 101,
            }
        },
        "genomic_annotation": {
            "status": "completed",
            "chromosome": "chr22",
            "position_1_based": 36201698,
            "gene_id": "ENSG00000100336.19",
            "gene_name": "APOL4",
            "gene_type": "protein_coding",
            "strand": "-",
            "primary_region": "splice_region",
            "region_labels": (
                "intron",
                "splice_region",
            ),
            "distance_to_tss": 3142,
            "transcription_start_site": 36204840,
        },
        "rna_analysis": {
            "status": "completed",
            "reference_rna_base": "U",
            "alternate_rna_base": "G",
            "position_1_based": 51,
            "strand": "-",
            "strand_source": "GENCODE",
            "reference_rna": "UCGGGUGAGUC",
            "mutated_rna": "UCGGGGGAGUC",
            "is_mature_mrna": False,
        },
        "expression_analysis": {
            "status": "completed",
            "model": "AlphaGenome",
            "gene_name": "APOL4",
            "tissue_query": "colon",
            "matched_tracks": 8,
            "summary": {
                "direction": "predicted_decrease",
                "raw_score": -1.9386,
                "quantile_score": -0.9999,
                "biosample_name": "sigmoid colon",
                "track_name": "total RNA-seq",
            },
        },
        "splicing_analysis": {
            "status": "completed",
            "model": "AlphaGenome",
            "target_gene": "APOL4",
            "target_gene_summary": {
                "splice_sites": 0.9726,
                "splice_site_usage": 0.8632,
                "splice_junctions": 7.46875,
                "merged_splicing_score": 3.3296,
            },
            "strongest_target_gene_records": [],
        },
        "tf_binding_analysis": {
            "status": "completed",
            "tool": "FIMO",
            "p_value_threshold": 1e-4,
            "total_fimo_hits": 49,
            "changed_sites": 5,
            "counts": {
                "gained": 4,
                "lost": 1,
                "strengthened": 0,
                "weakened": 0,
                "unchanged": 0,
            },
            "changes": [
                {
                    "change_type": "gained",
                    "motif_alt_id": "ZNF460",
                },
                {
                    "change_type": "lost",
                    "motif_alt_id": "SREBF1",
                },
            ],
        },
        "warnings": [],
        "scientific_notes": [
            "Computational predictions only."
        ],
    }

    report = build_analysis_report(
        pipeline_result
    )

    assert report["status"] == "completed"
    assert report["mutation"]["label"] == (
        "A>C at sequence position 51"
    )

    assert report[
        "genomic_context"
    ]["gene_name"] == "APOL4"

    assert report[
        "expression_effect"
    ]["direction"] == "predicted_decrease"

    assert report[
        "splicing_effect"
    ]["merged_splicing_score"] == pytest.approx(
        3.3296
    )

    assert report[
        "tf_binding_effect"
    ]["gained_factors"] == ["ZNF460"]

    assert report[
        "tf_binding_effect"
    ]["lost_factors"] == ["SREBF1"]

    assert len(report["key_findings"]) == 4


def test_report_supports_sequence_only_mode() -> None:
    pipeline_result = {
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
            "status": "not_run",
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
        "scientific_notes": [],
    }

    report = build_analysis_report(
        pipeline_result
    )

    assert report["mutation"]["reference_base"] == "A"
    assert report["mutation"]["alternate_base"] == "C"
    assert report["genomic_context"] is None
    assert report["expression_effect"] is None
    assert report["splicing_effect"] is None