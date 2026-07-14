"""Parse and filter FIMO motif-search results."""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class FimoHit:
    """Represent one motif match reported by FIMO."""

    motif_id: str
    motif_alt_id: str
    sequence_name: str

    start: int
    stop: int
    strand: str

    score: float
    p_value: float
    q_value: float | None

    matched_sequence: str

    @property
    def length(self) -> int:
        """Return the motif-match length."""
        return self.stop - self.start + 1

    def overlaps_position(
        self,
        position_1_based: int,
    ) -> bool:
        """Check whether this hit covers a 1-based sequence position."""
        return self.start <= position_1_based <= self.stop


def _optional_float(value: object) -> float | None:
    """Convert a possibly missing dataframe value to float."""
    if pd.isna(value):
        return None

    return float(value)


def parse_fimo_tsv(
    fimo_tsv_path: str | Path,
) -> list[FimoHit]:
    """Read FIMO's tab-separated output file."""
    path = Path(fimo_tsv_path)

    if not path.exists():
        raise FileNotFoundError(
            f"FIMO result file was not found: {path}"
        )

    dataframe = pd.read_csv(
        path,
        sep="\t",
        comment="#",
    )

    required_columns = {
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "start",
        "stop",
        "strand",
        "score",
        "p-value",
        "q-value",
        "matched_sequence",
    }

    missing_columns = required_columns - set(dataframe.columns)

    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))

        raise ValueError(
            f"FIMO output is missing required columns: {missing_text}"
        )

    hits: list[FimoHit] = []

    for _, row in dataframe.iterrows():
        hits.append(
            FimoHit(
                motif_id=str(row["motif_id"]),
                motif_alt_id=str(row["motif_alt_id"]),
                sequence_name=str(row["sequence_name"]),
                start=int(row["start"]),
                stop=int(row["stop"]),
                strand=str(row["strand"]),
                score=float(row["score"]),
                p_value=float(row["p-value"]),
                q_value=_optional_float(row["q-value"]),
                matched_sequence=str(row["matched_sequence"]),
            )
        )

    return hits


def filter_hits_by_sequence(
    hits: list[FimoHit],
    sequence_name: str,
) -> list[FimoHit]:
    """Return hits belonging to one FASTA sequence."""
    return [
        hit
        for hit in hits
        if hit.sequence_name == sequence_name
    ]


def filter_hits_overlapping_position(
    hits: list[FimoHit],
    position_1_based: int,
) -> list[FimoHit]:
    """Return motif hits covering a selected sequence position."""
    if isinstance(position_1_based, bool) or not isinstance(
        position_1_based,
        int,
    ):
        raise TypeError("Sequence position must be an integer.")

    if position_1_based < 1:
        raise ValueError(
            "Sequence position must be greater than or equal to 1."
        )

    return [
        hit
        for hit in hits
        if hit.overlaps_position(position_1_based)
    ]
