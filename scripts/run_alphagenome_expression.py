"""Run one real AlphaGenome gene-expression analysis."""

from pprint import pprint

from src.models.alphagenome_client import (
    build_alphagenome_variant,
    create_alphagenome_client,
    filter_expression_scores,
    score_gene_expression,
    summarize_expression_scores,
)


def main() -> None:
    """Run the official APOL4 example variant."""
    model = create_alphagenome_client()

    variant = build_alphagenome_variant(
        chromosome="chr22",
        position_1_based=36201698,
        reference_bases="A",
        alternate_bases="C",
    )

    print("Sending AlphaGenome expression request...")

    all_scores = score_gene_expression(
        model=model,
        variant=variant,
    )

    print(
        f"Total RNA expression score rows: "
        f"{len(all_scores)}"
    )

    filtered_scores = filter_expression_scores(
        scores=all_scores,
        gene_name="APOL4",
        tissue_query="colon",
    )

    print(
        f"APOL4 colon-matched tracks: "
        f"{len(filtered_scores)}"
    )

    summary = summarize_expression_scores(
        filtered_scores
    )

    print("\nExpression summary:")
    pprint(summary)

    if not filtered_scores.empty:
        columns_to_show = [
            column
            for column in (
                "gene_name",
                "biosample_name",
                "gtex_tissue",
                "ontology_curie",
                "track_name",
                "raw_score",
                "quantile_score",
            )
            if column in filtered_scores.columns
        ]

        print("\nMatched tracks:")
        print(
            filtered_scores[columns_to_show]
            .sort_values(
                by="raw_score",
                key=lambda column: column.abs(),
                ascending=False,
            )
            .head(10)
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()