"""Build a unified scientific report from pipeline results."""

from typing import Any


def _completed_section(
    pipeline_result: dict[str, Any],
    section_name: str,
) -> dict[str, Any] | None:
    """Return a completed analysis section."""
    section = pipeline_result.get(section_name)

    if not isinstance(section, dict):
        return None

    if section.get("status") != "completed":
        return None

    return section


def _build_mutation_summary(
    pipeline_result: dict[str, Any],
) -> dict[str, Any]:
    """Build the basic DNA mutation summary."""
    preprocessing = pipeline_result.get(
        "preprocessing",
        {},
    )

    mutation = preprocessing.get(
        "mutation",
        {},
    )

    reference_base = mutation.get("reference_base")
    alternate_base = mutation.get("alternate_base")
    position = mutation.get("position_1_based")

    mutation_label = None

    if (
        reference_base is not None
        and alternate_base is not None
        and position is not None
    ):
        mutation_label = (
            f"{reference_base}>{alternate_base} "
            f"at sequence position {position}"
        )

    return {
        "mutation_type": mutation.get(
            "mutation_type"
        ),
        "position_0_based": mutation.get(
            "position_0_based"
        ),
        "position_1_based": position,
        "reference_base": reference_base,
        "alternate_base": alternate_base,
        "sequence_length": mutation.get(
            "sequence_length"
        ),
        "label": mutation_label,
    }


def _build_genomic_summary(
    pipeline_result: dict[str, Any],
) -> dict[str, Any] | None:
    """Build the genomic-context summary."""
    annotation = _completed_section(
        pipeline_result,
        "genomic_annotation",
    )

    if annotation is None:
        return None

    chromosome = annotation.get("chromosome")
    genomic_position = annotation.get(
        "position_1_based"
    )

    genomic_variant = None

    if (
        chromosome is not None
        and genomic_position is not None
    ):
        genomic_variant = (
            f"{chromosome}:{genomic_position}"
        )

    return {
        "genomic_variant": genomic_variant,
        "gene_id": annotation.get("gene_id"),
        "gene_name": annotation.get("gene_name"),
        "gene_type": annotation.get("gene_type"),
        "strand": annotation.get("strand"),
        "primary_region": annotation.get(
            "primary_region"
        ),
        "region_labels": list(
            annotation.get("region_labels", ())
        ),
        "distance_to_tss": annotation.get(
            "distance_to_tss"
        ),
        "transcription_start_site": annotation.get(
            "transcription_start_site"
        ),
    }


def _build_rna_summary(
    pipeline_result: dict[str, Any],
) -> dict[str, Any] | None:
    """Build the direct RNA-transcription summary."""
    rna = _completed_section(
        pipeline_result,
        "rna_analysis",
    )

    if rna is None:
        return None

    reference_base = rna.get(
        "reference_rna_base"
    )

    alternate_base = rna.get(
        "alternate_rna_base"
    )

    change = None

    if (
        reference_base is not None
        and alternate_base is not None
    ):
        change = (
            f"{reference_base}>{alternate_base}"
        )

    return {
        "change": change,
        "position_1_based": rna.get(
            "position_1_based"
        ),
        "strand": rna.get("strand"),
        "strand_source": rna.get(
            "strand_source"
        ),
        "reference_rna": rna.get(
            "reference_rna"
        ),
        "mutated_rna": rna.get(
            "mutated_rna"
        ),
        "is_mature_mrna": rna.get(
            "is_mature_mrna",
            False,
        ),
    }


def _build_expression_summary(
    pipeline_result: dict[str, Any],
) -> dict[str, Any] | None:
    """Build the AlphaGenome expression summary."""
    expression = _completed_section(
        pipeline_result,
        "expression_analysis",
    )

    if expression is None:
        return None

    strongest = expression.get(
        "summary",
        {},
    )

    return {
        "model": expression.get("model"),
        "gene_name": expression.get(
            "gene_name"
        ),
        "tissue_query": expression.get(
            "tissue_query"
        ),
        "direction": strongest.get(
            "direction"
        ),
        "raw_score": strongest.get(
            "raw_score"
        ),
        "quantile_score": strongest.get(
            "quantile_score"
        ),
        "strongest_biosample": strongest.get(
            "biosample_name"
        ),
        "strongest_track": strongest.get(
            "track_name"
        ),
        "matched_tracks": expression.get(
            "matched_tracks"
        ),
    }


def _build_splicing_summary(
    pipeline_result: dict[str, Any],
) -> dict[str, Any] | None:
    """Build the AlphaGenome splicing summary."""
    splicing = _completed_section(
        pipeline_result,
        "splicing_analysis",
    )

    if splicing is None:
        return None

    target_summary = splicing.get(
        "target_gene_summary"
    )

    if not isinstance(target_summary, dict):
        target_summary = splicing.get(
            "global_summary",
            {},
        )

    return {
        "model": splicing.get("model"),
        "target_gene": splicing.get(
            "target_gene"
        ),
        "splice_sites": target_summary.get(
            "splice_sites"
        ),
        "splice_site_usage": target_summary.get(
            "splice_site_usage"
        ),
        "splice_junctions": target_summary.get(
            "splice_junctions"
        ),
        "merged_splicing_score": (
            target_summary.get(
                "merged_splicing_score"
            )
        ),
        "strongest_records": splicing.get(
            "strongest_target_gene_records",
            [],
        ),
    }


def _build_tf_binding_summary(
    pipeline_result: dict[str, Any],
) -> dict[str, Any] | None:
    """Build the FIMO TF-motif summary."""
    tf_binding = _completed_section(
        pipeline_result,
        "tf_binding_analysis",
    )

    if tf_binding is None:
        return None

    counts = tf_binding.get(
        "counts",
        {},
    )

    changes = tf_binding.get(
        "changes",
        [],
    )

    gained_factors = [
        change.get("motif_alt_id")
        for change in changes
        if change.get("change_type") == "gained"
    ]

    lost_factors = [
        change.get("motif_alt_id")
        for change in changes
        if change.get("change_type") == "lost"
    ]

    strengthened_factors = [
        change.get("motif_alt_id")
        for change in changes
        if change.get("change_type")
        == "strengthened"
    ]

    weakened_factors = [
        change.get("motif_alt_id")
        for change in changes
        if change.get("change_type")
        == "weakened"
    ]

    return {
        "tool": tf_binding.get("tool"),
        "p_value_threshold": tf_binding.get(
            "p_value_threshold"
        ),
        "total_fimo_hits": tf_binding.get(
            "total_fimo_hits"
        ),
        "changed_sites": tf_binding.get(
            "changed_sites"
        ),
        "counts": counts,
        "gained_factors": gained_factors,
        "lost_factors": lost_factors,
        "strengthened_factors": (
            strengthened_factors
        ),
        "weakened_factors": weakened_factors,
    }


def _build_key_findings(
    genomic: dict[str, Any] | None,
    expression: dict[str, Any] | None,
    splicing: dict[str, Any] | None,
    tf_binding: dict[str, Any] | None,
) -> list[str]:
    """Create human-readable key findings."""
    findings: list[str] = []

    if genomic is not None:
        gene_name = genomic.get("gene_name")
        region = genomic.get("primary_region")

        if gene_name and region:
            findings.append(
                f"The variant is associated with "
                f"{gene_name} and is located in "
                f"a predicted {region} region."
            )

    if expression is not None:
        direction = expression.get("direction")
        gene_name = expression.get("gene_name")
        biosample = expression.get(
            "strongest_biosample"
        )

        if direction and gene_name:
            readable_direction = direction.replace(
                "_",
                " ",
            )

            finding = (
                f"AlphaGenome predicted a "
                f"{readable_direction} in "
                f"{gene_name} expression"
            )

            if biosample:
                finding += (
                    f", with the strongest matched "
                    f"signal in {biosample}"
                )

            findings.append(f"{finding}.")

    if splicing is not None:
        merged_score = splicing.get(
            "merged_splicing_score"
        )

        if merged_score is not None:
            findings.append(
                "AlphaGenome detected a predicted "
                "splicing change with a merged score "
                f"of {merged_score:.4f}."
            )

    if tf_binding is not None:
        counts = tf_binding.get("counts", {})

        gained = counts.get("gained", 0)
        lost = counts.get("lost", 0)
        strengthened = counts.get(
            "strengthened",
            0,
        )
        weakened = counts.get(
            "weakened",
            0,
        )

        findings.append(
            "FIMO detected "
            f"{gained} gained, "
            f"{lost} lost, "
            f"{strengthened} strengthened, and "
            f"{weakened} weakened motif matches "
            "overlapping the mutation."
        )

    if not findings:
        findings.append(
            "No optional contextual analyses were completed."
        )

    return findings


def build_analysis_report(
    pipeline_result: dict[str, Any],
) -> dict[str, Any]:
    """Build one unified report from all pipeline sections."""
    mutation = _build_mutation_summary(
        pipeline_result
    )

    genomic = _build_genomic_summary(
        pipeline_result
    )

    rna = _build_rna_summary(
        pipeline_result
    )

    expression = _build_expression_summary(
        pipeline_result
    )

    splicing = _build_splicing_summary(
        pipeline_result
    )

    tf_binding = _build_tf_binding_summary(
        pipeline_result
    )

    key_findings = _build_key_findings(
        genomic=genomic,
        expression=expression,
        splicing=splicing,
        tf_binding=tf_binding,
    )

    limitations = list(
        pipeline_result.get(
            "scientific_notes",
            [],
        )
    )

    limitations.extend(
        [
            (
                "The report contains computational predictions "
                "and does not constitute a medical diagnosis."
            ),
            (
                "Predicted effects should be validated with "
                "experimental and clinical evidence before "
                "biological conclusions are made."
            ),
        ]
    )

    return {
        "status": "completed",
        "report_type": (
            "integrated_dna_mutation_analysis"
        ),
        "project_version": pipeline_result.get(
            "project_version"
        ),
        "mutation": mutation,
        "genomic_context": genomic,
        "rna_change": rna,
        "expression_effect": expression,
        "splicing_effect": splicing,
        "tf_binding_effect": tf_binding,
        "key_findings": key_findings,
        "warnings": pipeline_result.get(
            "warnings",
            [],
        ),
        "limitations": limitations,
    }
