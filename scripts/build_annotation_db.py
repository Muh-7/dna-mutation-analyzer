"""Build the local GENCODE annotation database."""

from pathlib import Path

from src.annotation.database import build_annotation_database


PROJECT_ROOT = Path(__file__).resolve().parents[1]

GTF_PATH = (
    PROJECT_ROOT
    / "data"
    / "annotations"
    / "gencode.v50.basic.annotation.gtf"
)

DATABASE_PATH = (
    PROJECT_ROOT
    / "data"
    / "annotations"
    / "gencode.v50.basic.annotation.db"
)


def main() -> None:
    """Build the GENCODE SQLite database."""
    print("Building genomic annotation database...")
    print(f"GTF source: {GTF_PATH}")
    print(f"Database:   {DATABASE_PATH}")

    result = build_annotation_database(
        gtf_path=GTF_PATH,
        database_path=DATABASE_PATH,
        force=False,
    )

    print(f"Database ready: {result}")


if __name__ == "__main__":
    main()