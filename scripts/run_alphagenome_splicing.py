"""Run AlphaGenome splicing analysis for the APOL4 variant."""

from pprint import pprint

from Bio import SeqIO

from src.models.alphagenome_client import (
    analyze_variant_splicing,
    create_alphagenome_client,
)
from src.schemas.variant_input import (
    VariantInput,
)


FASTA_PATH = (
    "data/tf_binding/"
    "apol4_reference_mutated.fa"
)


def main() -> None:
    """Run the APOL4 splicing example."""
    sequences = {
        record.id: str(record.seq)
        for record in SeqIO.parse(
            FASTA_PATH,
            "fasta",
        )
    }

    variant = VariantInput(
        analysis_mode="genomic",
        reference_sequence=(
            sequences["reference"]
        ),
        mutated_sequence=(
            sequences["mutated"]
        ),
        chromosome="chr22",
        genomic_position=36201698,
        gene_name="APOL4",
        tissue="colon",
        strand="-",
        flank_size=50,
    )

    model = create_alphagenome_client()

    result = analyze_variant_splicing(
        model=model,
        variant_input=variant,
    )

    print("\nAlphaGenome splicing analysis:\n")
    pprint(result)


if __name__ == "__main__":
    main()