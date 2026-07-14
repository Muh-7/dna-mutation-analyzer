"""Create and access a genomic annotation database from a GTF file."""

from pathlib import Path

import gffutils


def build_annotation_database(
    gtf_path: str | Path,
    database_path: str | Path,
    force: bool = False,
) -> Path:
    """
    Convert a GTF annotation file into a SQLite database.

    Parameters
    ----------
    gtf_path:
        Path to the source GTF file.

    database_path:
        Path where the SQLite database will be stored.

    force:
        Whether an existing database should be overwritten.

    Returns
    -------
    Path
        Path to the created or existing database.

    Raises
    ------
    FileNotFoundError
        If the GTF file does not exist.
    """
    gtf = Path(gtf_path)
    database = Path(database_path)

    if not gtf.exists():
        raise FileNotFoundError(
            f"GTF annotation file was not found: {gtf}"
        )

    if database.exists() and not force:
        return database

    database.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    gffutils.create_db(
        data=str(gtf),
        dbfn=str(database),
        force=force,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )

    return database


def open_annotation_database(
    database_path: str | Path,
) -> gffutils.FeatureDB:
    """
    Open an existing annotation SQLite database.

    Raises
    ------
    FileNotFoundError
        If the database does not exist.
    """
    database = Path(database_path)

    if not database.exists():
        raise FileNotFoundError(
            f"Annotation database was not found: {database}"
        )

    return gffutils.FeatureDB(
        str(database),
        keep_order=True,
    )