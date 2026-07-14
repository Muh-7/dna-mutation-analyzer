"""Run FIMO on reference and mutated DNA sequences."""

from pathlib import Path
import shutil
import subprocess
import tempfile

from src.preprocessing.sequence_validator import (
    validate_sequence,
)
from src.tf_binding.fimo_parser import (
    FimoHit,
    parse_fimo_tsv,
)


def resolve_fimo_executable(
    fimo_executable: str = "fimo",
) -> str:
    """Find the FIMO executable on the current system."""
    resolved_path = shutil.which(fimo_executable)

    if resolved_path is not None:
        return resolved_path

    supplied_path = Path(
        fimo_executable
    ).expanduser()

    if supplied_path.is_file():
        return str(supplied_path.resolve())

    raise FileNotFoundError(
        "FIMO executable was not found. "
        f"Requested executable: {fimo_executable}"
    )


def validate_p_value_threshold(
    threshold: float,
) -> float:
    """Validate a FIMO p-value threshold."""
    if isinstance(threshold, bool) or not isinstance(
        threshold,
        (int, float),
    ):
        raise TypeError(
            "FIMO p-value threshold must be numeric."
        )

    normalized_threshold = float(threshold)

    if not 0 < normalized_threshold <= 1:
        raise ValueError(
            "FIMO p-value threshold must be greater than 0 "
            "and less than or equal to 1."
        )

    return normalized_threshold


def _wrap_fasta_sequence(
    sequence: str,
    line_width: int = 80,
) -> str:
    """Wrap a DNA sequence into FASTA-style lines."""
    return "\n".join(
        sequence[index:index + line_width]
        for index in range(
            0,
            len(sequence),
            line_width,
        )
    )


def write_sequence_pair_fasta(
    reference_sequence: str,
    mutated_sequence: str,
    output_path: str | Path,
) -> Path:
    """Write reference and mutated DNA to one FASTA file."""
    reference = validate_sequence(
        reference_sequence,
        allow_n=False,
    )

    mutated = validate_sequence(
        mutated_sequence,
        allow_n=False,
    )

    if len(reference) != len(mutated):
        raise ValueError(
            "Reference and mutated sequences must have "
            "the same length for SNV TF-binding analysis."
        )

    path = Path(output_path)
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fasta_content = (
        ">reference\n"
        f"{_wrap_fasta_sequence(reference)}\n"
        ">mutated\n"
        f"{_wrap_fasta_sequence(mutated)}\n"
    )

    path.write_text(
        fasta_content,
        encoding="utf-8",
    )

    return path


def run_fimo(
    reference_sequence: str,
    mutated_sequence: str,
    motif_database_path: str | Path,
    fimo_executable: str = "fimo",
    p_value_threshold: float = 1e-4,
) -> list[FimoHit]:
    """
    Run FIMO and return its parsed motif hits.

    A temporary working directory is removed automatically
    after the results have been parsed.
    """
    motif_path = Path(
        motif_database_path
    ).expanduser()

    if not motif_path.exists():
        raise FileNotFoundError(
            "Motif database was not found: "
            f"{motif_path}"
        )

    if not motif_path.is_file():
        raise ValueError(
            "Motif database path must point to a file."
        )

    executable = resolve_fimo_executable(
        fimo_executable
    )

    threshold = validate_p_value_threshold(
        p_value_threshold
    )

    with tempfile.TemporaryDirectory(
        prefix="dna_mutation_fimo_",
    ) as temporary_directory:
        working_directory = Path(
            temporary_directory
        )

        sequence_path = write_sequence_pair_fasta(
            reference_sequence=reference_sequence,
            mutated_sequence=mutated_sequence,
            output_path=(
                working_directory
                / "sequence_pair.fa"
            ),
        )

        output_directory = (
            working_directory
            / "fimo_output"
        )

        command = [
            executable,
            "--oc",
            str(output_directory),
            "--thresh",
            f"{threshold:g}",
            str(motif_path.resolve()),
            str(sequence_path),
        ]

        try:
            subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
            )

        except subprocess.CalledProcessError as error:
            error_message = (
                error.stderr.strip()
                or error.stdout.strip()
                or "FIMO returned no diagnostic message."
            )

            raise RuntimeError(
                "FIMO execution failed:\n"
                f"{error_message}"
            ) from error

        result_path = (
            output_directory
            / "fimo.tsv"
        )

        if not result_path.exists():
            raise RuntimeError(
                "FIMO finished without producing fimo.tsv."
            )

        return parse_fimo_tsv(
            result_path
        )
