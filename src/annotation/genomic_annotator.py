"""Annotate a genomic position using a GENCODE database."""

from dataclasses import dataclass
from pathlib import Path

from gffutils import Feature, FeatureDB

from src.annotation.database import open_annotation_database


@dataclass(frozen=True)
class GenomicAnnotation:
    """
    Store annotation information for one genomic position.
    """

    chromosome: str
    position_1_based: int

    primary_region: str
    region_labels: tuple[str, ...]

    gene_id: str
    gene_name: str
    gene_type: str | None
    strand: str

    gene_start: int
    gene_end: int
    transcription_start_site: int
    distance_to_tss: int

    transcript_ids: tuple[str, ...]
    exon_numbers: tuple[int, ...]

    promoter_upstream_bp: int
    splice_window_bp: int


def _first_attribute(
    feature: Feature,
    attribute_name: str,
) -> str | None:
    """Return the first value of a GTF feature attribute."""
    values = feature.attributes.get(attribute_name)

    if not values:
        return None

    return values[0]


def _transcription_start_site(gene: Feature) -> int:
    """Return the gene TSS according to its strand."""
    if gene.strand == "+":
        return gene.start

    return gene.end


def _distance_to_tss(
    position: int,
    tss: int,
    strand: str,
) -> int:
    """
    Calculate strand-aware distance from the transcription start site.

    Negative values are upstream of the TSS.
    Positive values are downstream of the TSS.
    """
    if strand == "+":
        return position - tss

    return tss - position


def _promoter_interval(
    gene: Feature,
    upstream_bp: int,
) -> tuple[int, int]:
    """Calculate a strand-aware promoter interval."""
    tss = _transcription_start_site(gene)

    if gene.strand == "+":
        return max(1, tss - upstream_bp), tss - 1

    return tss + 1, tss + upstream_bp


def _contains_position(
    feature: Feature,
    position: int,
) -> bool:
    """Check whether a 1-based position overlaps a feature."""
    return feature.start <= position <= feature.end


def _is_near_exon_boundary(
    exon: Feature,
    position: int,
    window_bp: int,
) -> bool:
    """Check whether a position is close to an exon boundary."""
    distance_to_start = abs(position - exon.start)
    distance_to_end = abs(position - exon.end)

    return min(distance_to_start, distance_to_end) <= window_bp


def _choose_primary_region(
    labels: set[str],
) -> str:
    """Select one primary region using project priority rules."""
    priority = (
        "splice_region",
        "UTR",
        "CDS",
        "exon",
        "intron",
        "gene_body",
        "promoter",
        "intergenic",
    )

    for region in priority:
        if region in labels:
            return region

    return "intergenic"


def _find_target_gene(
    database: FeatureDB,
    chromosome: str,
    position: int,
    gene_name: str | None,
    promoter_upstream_bp: int,
) -> Feature:
    """Find the gene overlapping the position or its promoter."""
    search_start = max(
        1,
        position - promoter_upstream_bp,
    )

    search_end = position + promoter_upstream_bp

    nearby_genes = list(
        database.region(
            region=(
                chromosome,
                search_start,
                search_end,
            ),
            featuretype="gene",
        )
    )

    candidates: list[Feature] = []

    for gene in nearby_genes:
        current_gene_name = _first_attribute(
            gene,
            "gene_name",
        )

        if (
            gene_name is not None
            and current_gene_name is not None
            and current_gene_name.upper() != gene_name.upper()
        ):
            continue

        promoter_start, promoter_end = _promoter_interval(
            gene,
            promoter_upstream_bp,
        )

        overlaps_gene = _contains_position(
            gene,
            position,
        )

        overlaps_promoter = (
            promoter_start
            <= position
            <= promoter_end
        )

        if overlaps_gene or overlaps_promoter:
            candidates.append(gene)

    if not candidates:
        requested_gene = (
            f" for gene {gene_name}"
            if gene_name is not None
            else ""
        )

        raise LookupError(
            f"No gene or promoter annotation was found at "
            f"{chromosome}:{position}{requested_gene}."
        )

    overlapping_genes = [
        gene
        for gene in candidates
        if _contains_position(gene, position)
    ]

    if len(overlapping_genes) == 1:
        return overlapping_genes[0]

    if len(overlapping_genes) > 1:
        candidates = overlapping_genes

    if len(candidates) > 1:
        names = [
            _first_attribute(gene, "gene_name")
            or gene.id
            for gene in candidates
        ]

        raise ValueError(
            "Multiple genes match this position: "
            f"{', '.join(names)}. Provide a gene name."
        )

    return candidates[0]


def annotate_genomic_position(
    database_path: str | Path,
    chromosome: str,
    position_1_based: int,
    gene_name: str | None = None,
    promoter_upstream_bp: int = 2000,
    splice_window_bp: int = 2,
) -> GenomicAnnotation:
    """
    Annotate one 1-based genomic position.

    The function determines whether the position is in a promoter,
    exon, CDS, UTR, intron, or splice-region context.
    """
    if position_1_based < 1:
        raise ValueError(
            "Genomic position must be greater than or equal to 1."
        )

    if promoter_upstream_bp < 0:
        raise ValueError(
            "Promoter upstream size cannot be negative."
        )

    if splice_window_bp < 0:
        raise ValueError(
            "Splice window size cannot be negative."
        )

    database = open_annotation_database(database_path)

    gene = _find_target_gene(
        database=database,
        chromosome=chromosome,
        position=position_1_based,
        gene_name=gene_name,
        promoter_upstream_bp=promoter_upstream_bp,
    )

    labels: set[str] = set()

    inside_gene = _contains_position(
        gene,
        position_1_based,
    )

    promoter_start, promoter_end = _promoter_interval(
        gene,
        promoter_upstream_bp,
    )

    if (
        not inside_gene
        and promoter_start
        <= position_1_based
        <= promoter_end
    ):
        labels.add("promoter")

    transcripts = list(
        database.children(
            gene,
            featuretype="transcript",
            order_by="start",
        )
    )

    affected_transcripts = [
        transcript
        for transcript in transcripts
        if _contains_position(
            transcript,
            position_1_based,
        )
    ]

    exons = list(
        database.children(
            gene,
            featuretype="exon",
            order_by="start",
        )
    )

    coding_regions = list(
        database.children(
            gene,
            featuretype="CDS",
            order_by="start",
        )
    )

    utr_regions = list(
        database.children(
            gene,
            featuretype="UTR",
            order_by="start",
        )
    )

    overlapping_exons = [
        exon
        for exon in exons
        if _contains_position(
            exon,
            position_1_based,
        )
    ]

    overlapping_cds = [
        cds
        for cds in coding_regions
        if _contains_position(
            cds,
            position_1_based,
        )
    ]

    overlapping_utrs = [
        utr
        for utr in utr_regions
        if _contains_position(
            utr,
            position_1_based,
        )
    ]

    near_splice_boundary = any(
        _is_near_exon_boundary(
            exon,
            position_1_based,
            splice_window_bp,
        )
        for exon in exons
    )

    if inside_gene:
        if overlapping_exons:
            labels.add("exon")

        elif affected_transcripts:
            # The position is inside at least one transcript,
            # but it is outside all of its exons.
            labels.add("intron")

        else:
            # The position is inside the broad gene boundaries,
            # but outside all annotated transcripts.
            labels.add("gene_body")

        if overlapping_cds:
            labels.add("CDS")

        if overlapping_utrs:
            labels.add("UTR")

        if affected_transcripts and near_splice_boundary:
            labels.add("splice_region")

    if not labels:
        labels.add("intergenic")

    exon_numbers: set[int] = set()

    for exon in overlapping_exons:
        for value in exon.attributes.get(
            "exon_number",
            [],
        ):
            try:
                exon_numbers.add(int(value))
            except ValueError:
                continue

    tss = _transcription_start_site(gene)

    return GenomicAnnotation(
        chromosome=chromosome,
        position_1_based=position_1_based,
        primary_region=_choose_primary_region(labels),
        region_labels=tuple(sorted(labels)),
        gene_id=(
            _first_attribute(gene, "gene_id")
            or gene.id
        ),
        gene_name=(
            _first_attribute(gene, "gene_name")
            or gene.id
        ),
        gene_type=_first_attribute(
            gene,
            "gene_type",
        ),
        strand=gene.strand,
        gene_start=gene.start,
        gene_end=gene.end,
        transcription_start_site=tss,
        distance_to_tss=_distance_to_tss(
            position=position_1_based,
            tss=tss,
            strand=gene.strand,
        ),
        transcript_ids=tuple(
            sorted(
                {
                    (
                        _first_attribute(
                            transcript,
                            "transcript_id",
                        )
                        or transcript.id
                    )
                    for transcript in affected_transcripts
                }
            )
        ),
        exon_numbers=tuple(
            sorted(exon_numbers)
        ),
        promoter_upstream_bp=promoter_upstream_bp,
        splice_window_bp=splice_window_bp,
    )