"""Tests for the main DNA mutation analysis pipeline."""

from pathlib import Path

import pytest

from src.annotation.database import build_annotation_database

from src.pipeline.analysis_pipeline import run_analysis_pipeline
from src.schemas.variant_input import VariantInput
from src.annotation.genomic_annotator import GenomicAnnotation



PIPELINE_TEST_GTF = """\
chr1\tTEST\tgene\t100\t500\t.\t+\t.\tgene_id "GENE1"; gene_type "protein_coding"; gene_name "TEST1";
chr1\tTEST\ttranscript\t100\t500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_type "protein_coding"; gene_name "TEST1";
chr1\tTEST\texon\t100\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_number "1"; exon_id "EXON1"; gene_name "TEST1";
chr1\tTEST\tCDS\t120\t200\t.\t+\t0\tgene_id "GENE1"; transcript_id "TX1"; exon_number "1"; exon_id "EXON1"; gene_name "TEST1";
chr1\tTEST\texon\t300\t500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_number "2"; exon_id "EXON2"; gene_name "TEST1";
"""





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

    assert result["project_version"] == "0.7.0"
    assert result["analysis_status"] == "completed"
    assert result["final_report"]["status"] == "completed"
    
    
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
    
    
@pytest.fixture()
def pipeline_annotation_database(
    tmp_path: Path,
) -> Path:
    """Create a small temporary annotation database."""
    gtf_path = tmp_path / "pipeline_test.gtf"
    database_path = tmp_path / "pipeline_test.db"

    gtf_path.write_text(
        PIPELINE_TEST_GTF,
        encoding="utf-8",
    )

    return build_annotation_database(
        gtf_path=gtf_path,
        database_path=database_path,
        force=True,
    )
    


def test_pipeline_includes_genomic_annotation(
    pipeline_annotation_database: Path,
) -> None:
    variant = VariantInput(
        analysis_mode="genomic",
        reference_sequence="AACCGTACGT",
        mutated_sequence="AACCGGACGT",
        chromosome="1",
        genomic_position=150,
        gene_name="TEST1",
        tissue="test_tissue",
        strand="+",
        flank_size=2,
    )

    result = run_analysis_pipeline(
        variant=variant,
        annotation_database_path=pipeline_annotation_database,
    )

    annotation = result["genomic_annotation"]

    assert annotation["status"] == "completed"
    assert annotation["chromosome"] == "chr1"
    assert annotation["position_1_based"] == 150
    assert annotation["gene_name"] == "TEST1"
    assert annotation["primary_region"] == "CDS"
    assert "exon" in annotation["region_labels"]
    assert "CDS" in annotation["region_labels"]
    assert annotation["exon_numbers"] == (1,)





def test_raw_pipeline_skips_genomic_annotation() -> None:
    variant = VariantInput(
        reference_sequence="ACGT",
        mutated_sequence="AGGT",
    )

    result = run_analysis_pipeline(variant)

    annotation = result["genomic_annotation"]

    assert annotation["status"] == "not_run"
    assert "genomic" in annotation["reason"].lower()
    
    
    
    
def test_expression_analysis_is_disabled_by_default() -> None:
    variant = VariantInput(
        analysis_mode="genomic",
        reference_sequence="AACCGAACTG",
        mutated_sequence="AACCGCACTG",
        chromosome="22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="+",
    )

    result = run_analysis_pipeline(variant)

    expression = result["expression_analysis"]

    assert expression["status"] == "not_run"
    assert "not requested" in expression["reason"]
    



def test_pipeline_includes_expression_analysis(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    variant = VariantInput(
        analysis_mode="genomic",
        reference_sequence="AACCGAACTG",
        mutated_sequence="AACCGCACTG",
        chromosome="22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="+",
    )

    fake_model = object()

    fake_expression_result = {
        "status": "completed",
        "model": "AlphaGenome",
        "output_type": "RNA_SEQ",
        "gene_name": "APOL4",
        "tissue_query": "colon",
        "total_score_rows": 13727,
        "matched_tracks": 8,
        "summary": {
            "status": "completed",
            "direction": "predicted_decrease",
            "raw_score": -1.938638687133789,
            "quantile_score": -0.9999998807907104,
        },
        "tracks": [],
    }

    def fake_analyze_variant_expression(
        model: object,
        variant_input: VariantInput,
    ) -> dict[str, object]:
        assert model is fake_model
        assert variant_input.gene_name == "APOL4"

        return fake_expression_result

    monkeypatch.setattr(
        (
            "src.pipeline.analysis_pipeline."
            "analyze_variant_expression"
        ),
        fake_analyze_variant_expression,
    )

    result = run_analysis_pipeline(
        variant=variant,
        run_expression_analysis=True,
        alphagenome_model=fake_model,
    )

    expression = result["expression_analysis"]

    assert expression["status"] == "completed"
    assert expression["matched_tracks"] == 8
    assert (
        expression["summary"]["direction"]
        == "predicted_decrease"
    )
    assert (
        expression["summary"]["raw_score"]
        == pytest.approx(-1.938638687133789)
    )




def test_gencode_strand_overrides_incorrect_input_strand(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    variant = VariantInput(
        analysis_mode="genomic",
        reference_sequence="AACCGAACTG",
        mutated_sequence="AACCGCACTG",
        chromosome="22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="+",
    )

    fake_annotation = GenomicAnnotation(
        chromosome="chr22",
        position_1_based=36201698,
        primary_region="splice_region",
        region_labels=("intron", "splice_region"),
        gene_id="ENSG_TEST",
        gene_name="APOL4",
        gene_type="protein_coding",
        strand="-",
        gene_start=36189124,
        gene_end=36204840,
        transcription_start_site=36204840,
        distance_to_tss=3142,
        transcript_ids=("ENST_TEST",),
        exon_numbers=(),
        promoter_upstream_bp=2000,
        splice_window_bp=2,
    )

    monkeypatch.setattr(
        (
            "src.pipeline.analysis_pipeline."
            "annotate_genomic_position"
        ),
        lambda **kwargs: fake_annotation,
    )

    result = run_analysis_pipeline(
        variant=variant,
        annotation_database_path="fake.db",
    )

    rna = result["rna_analysis"]

    assert rna["strand"] == "-"
    assert rna["strand_source"] == "GENCODE"
    assert rna["reference_rna"] == "CAGUUCGGUU"
    assert rna["mutated_rna"] == "CAGUGCGGUU"
    assert rna["position_1_based"] == 5
    assert rna["reference_rna_base"] == "U"
    assert rna["alternate_rna_base"] == "G"

    assert result["warnings"][0]["type"] == "strand_mismatch"
    



def test_pipeline_includes_tf_binding_analysis(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    motif_path = (
        tmp_path
        / "motifs.meme"
    )

    motif_path.write_text(
        "MEME version 4\n",
        encoding="utf-8",
    )

    variant = VariantInput(
        reference_sequence="AACCGAACTG",
        mutated_sequence="AACCGCACTG",
        flank_size=2,
    )

    fake_result = {
        "status": "completed",
        "tool": "FIMO",
        "counts": {
            "gained": 1,
            "lost": 1,
            "strengthened": 0,
            "weakened": 0,
            "unchanged": 0,
        },
        "changed_sites": 2,
        "changes": [],
    }

    def fake_analysis(
        reference_sequence: str,
        mutated_sequence: str,
        motif_database_path: str | Path,
        fimo_executable: str,
        p_value_threshold: float,
    ) -> dict[str, object]:
        assert reference_sequence == "CGAAC"
        assert mutated_sequence == "CGCAC"
        assert motif_database_path == motif_path
        assert fimo_executable == "fimo"
        assert p_value_threshold == pytest.approx(
            1e-4
        )

        return fake_result

    monkeypatch.setattr(
        (
            "src.pipeline.analysis_pipeline."
            "analyze_tf_binding_sequences"
        ),
        fake_analysis,
    )

    result = run_analysis_pipeline(
        variant=variant,
        run_tf_binding_analysis=True,
        motif_database_path=motif_path,
    )

    tf_result = result[
        "tf_binding_analysis"
    ]

    assert tf_result["status"] == "completed"
    assert tf_result["counts"]["gained"] == 1
    assert tf_result["counts"]["lost"] == 1
    
    
    
def test_pipeline_includes_splicing_analysis(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    variant = VariantInput(
        analysis_mode="genomic",
        reference_sequence="GACTCACCCGA",
        mutated_sequence="GACTCCCCCGA",
        chromosome="chr22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="-",
    )

    fake_model = object()

    fake_splicing_result = {
        "status": "completed",
        "analysis_type": "merged_splicing",
        "target_gene": "APOL4",
        "target_gene_summary": {
            "splice_sites": 0.9726,
            "splice_site_usage": 0.8632,
            "splice_junctions": 7.46875,
            "merged_splicing_score": 3.3296,
        },
    }

    def fake_analyze_variant_splicing(
        model: object,
        variant_input: VariantInput,
    ) -> dict[str, object]:
        assert model is fake_model
        assert variant_input.gene_name == "APOL4"

        return fake_splicing_result

    monkeypatch.setattr(
        (
            "src.pipeline.analysis_pipeline."
            "analyze_variant_splicing"
        ),
        fake_analyze_variant_splicing,
    )

    result = run_analysis_pipeline(
        variant=variant,
        run_splicing_analysis=True,
        alphagenome_model=fake_model,
    )

    splicing = result["splicing_analysis"]

    assert splicing["status"] == "completed"

    assert splicing[
        "target_gene_summary"
    ]["merged_splicing_score"] == pytest.approx(
        3.3296
    )