"""Tests for genomic position annotation."""

from pathlib import Path

import pytest

from src.annotation.database import build_annotation_database
from src.annotation.genomic_annotator import (
    annotate_genomic_position,
)


TEST_GTF = """\
chr1\tTEST\tgene\t100\t500\t.\t+\t.\tgene_id "GENE1"; gene_type "protein_coding"; gene_name "TEST1";
chr1\tTEST\ttranscript\t100\t500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_type "protein_coding"; gene_name "TEST1"; transcript_type "protein_coding";
chr1\tTEST\texon\t100\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_number "1"; exon_id "EXON1"; gene_name "TEST1";
chr1\tTEST\tUTR\t100\t119\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_number "1"; exon_id "EXON1"; gene_name "TEST1";
chr1\tTEST\tCDS\t120\t200\t.\t+\t0\tgene_id "GENE1"; transcript_id "TX1"; exon_number "1"; exon_id "EXON1"; gene_name "TEST1";
chr1\tTEST\texon\t300\t500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_number "2"; exon_id "EXON2"; gene_name "TEST1";
chr1\tTEST\tCDS\t300\t450\t.\t+\t0\tgene_id "GENE1"; transcript_id "TX1"; exon_number "2"; exon_id "EXON2"; gene_name "TEST1";
chr1\tTEST\tUTR\t451\t500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_number "2"; exon_id "EXON2"; gene_name "TEST1";
chr2\tTEST\tgene\t1000\t1500\t.\t-\t.\tgene_id "GENE2"; gene_type "protein_coding"; gene_name "TEST2";
chr2\tTEST\ttranscript\t1000\t1500\t.\t-\t.\tgene_id "GENE2"; transcript_id "TX2"; gene_type "protein_coding"; gene_name "TEST2"; transcript_type "protein_coding";
chr2\tTEST\texon\t1000\t1200\t.\t-\t.\tgene_id "GENE2"; transcript_id "TX2"; exon_number "2"; exon_id "EXON3"; gene_name "TEST2";
chr2\tTEST\texon\t1400\t1500\t.\t-\t.\tgene_id "GENE2"; transcript_id "TX2"; exon_number "1"; exon_id "EXON4"; gene_name "TEST2";
chr3\tTEST\tgene\t100\t500\t.\t+\t.\tgene_id "GENE3"; gene_type "protein_coding"; gene_name "TEST3";
chr3\tTEST\ttranscript\t150\t450\t.\t+\t.\tgene_id "GENE3"; transcript_id "TX3"; gene_type "protein_coding"; gene_name "TEST3"; transcript_type "protein_coding";
chr3\tTEST\texon\t150\t200\t.\t+\t.\tgene_id "GENE3"; transcript_id "TX3"; exon_number "1"; exon_id "EXON5"; gene_name "TEST3";
chr3\tTEST\texon\t300\t450\t.\t+\t.\tgene_id "GENE3"; transcript_id "TX3"; exon_number "2"; exon_id "EXON6"; gene_name "TEST3";
"""


@pytest.fixture()
def annotation_database(
    tmp_path: Path,
) -> Path:
    """Create a tiny annotation database for each test."""
    gtf_path = tmp_path / "test.gtf"
    database_path = tmp_path / "test.db"

    gtf_path.write_text(
        TEST_GTF,
        encoding="utf-8",
    )

    return build_annotation_database(
        gtf_path=gtf_path,
        database_path=database_path,
        force=True,
    )


def test_cds_annotation(
    annotation_database: Path,
) -> None:
    result = annotate_genomic_position(
        database_path=annotation_database,
        chromosome="chr1",
        position_1_based=150,
        gene_name="TEST1",
    )

    assert result.primary_region == "CDS"
    assert "exon" in result.region_labels
    assert "CDS" in result.region_labels
    assert result.gene_name == "TEST1"
    assert result.exon_numbers == (1,)


def test_utr_annotation(
    annotation_database: Path,
) -> None:
    result = annotate_genomic_position(
        database_path=annotation_database,
        chromosome="chr1",
        position_1_based=110,
        gene_name="TEST1",
    )

    assert result.primary_region == "UTR"
    assert "UTR" in result.region_labels
    assert "exon" in result.region_labels


def test_intron_annotation(
    annotation_database: Path,
) -> None:
    result = annotate_genomic_position(
        database_path=annotation_database,
        chromosome="chr1",
        position_1_based=250,
        gene_name="TEST1",
    )

    assert result.primary_region == "intron"
    assert result.region_labels == ("intron",)


def test_splice_region_annotation(
    annotation_database: Path,
) -> None:
    result = annotate_genomic_position(
        database_path=annotation_database,
        chromosome="chr1",
        position_1_based=202,
        gene_name="TEST1",
        splice_window_bp=2,
    )

    assert result.primary_region == "splice_region"
    assert "splice_region" in result.region_labels
    assert "intron" in result.region_labels


def test_forward_strand_promoter(
    annotation_database: Path,
) -> None:
    result = annotate_genomic_position(
        database_path=annotation_database,
        chromosome="chr1",
        position_1_based=90,
        gene_name="TEST1",
        promoter_upstream_bp=20,
    )

    assert result.primary_region == "promoter"
    assert result.distance_to_tss == -10


def test_reverse_strand_promoter(
    annotation_database: Path,
) -> None:
    result = annotate_genomic_position(
        database_path=annotation_database,
        chromosome="chr2",
        position_1_based=1550,
        gene_name="TEST2",
        promoter_upstream_bp=100,
    )

    assert result.primary_region == "promoter"
    assert result.strand == "-"
    assert result.transcription_start_site == 1500
    assert result.distance_to_tss == -50


def test_missing_annotation_is_reported(
    annotation_database: Path,
) -> None:
    with pytest.raises(
        LookupError,
        match="No gene or promoter annotation",
    ):
        annotate_genomic_position(
            database_path=annotation_database,
            chromosome="chr1",
            position_1_based=10000,
            gene_name="TEST1",
        )


def test_position_must_be_positive(
    annotation_database: Path,
) -> None:
    with pytest.raises(
        ValueError,
        match="greater than or equal to 1",
    ):
        annotate_genomic_position(
            database_path=annotation_database,
            chromosome="chr1",
            position_1_based=0,
        )


def test_gene_body_outside_all_transcripts(
    annotation_database: Path,
) -> None:
    result = annotate_genomic_position(
        database_path=annotation_database,
        chromosome="chr3",
        position_1_based=110,
        gene_name="TEST3",
    )

    assert result.primary_region == "gene_body"
    assert result.region_labels == ("gene_body",)
    assert result.transcript_ids == ()
    assert result.exon_numbers == ()