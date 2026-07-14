"""Tests for automated FIMO execution."""

from pathlib import Path
import subprocess

import pytest

from src.tf_binding.fimo_runner import (
    resolve_fimo_executable,
    run_fimo,
    validate_p_value_threshold,
    write_sequence_pair_fasta,
)


def test_write_sequence_pair_fasta(
    tmp_path: Path,
) -> None:
    output_path = (
        tmp_path
        / "sequences.fa"
    )

    write_sequence_pair_fasta(
        reference_sequence="AACCGAACTG",
        mutated_sequence="AACCGCACTG",
        output_path=output_path,
    )

    content = output_path.read_text(
        encoding="utf-8"
    )

    assert ">reference" in content
    assert "AACCGAACTG" in content
    assert ">mutated" in content
    assert "AACCGCACTG" in content


def test_validate_p_value_threshold() -> None:
    assert validate_p_value_threshold(
        1e-4
    ) == pytest.approx(1e-4)

    with pytest.raises(
        ValueError,
        match="greater than 0",
    ):
        validate_p_value_threshold(0)


def test_missing_fimo_is_rejected(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "src.tf_binding.fimo_runner.shutil.which",
        lambda executable: None,
    )

    with pytest.raises(
        FileNotFoundError,
        match="was not found",
    ):
        resolve_fimo_executable(
            "missing-fimo"
        )


def test_run_fimo_parses_generated_output(
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

    monkeypatch.setattr(
        "src.tf_binding.fimo_runner.shutil.which",
        lambda executable: "/fake/fimo",
    )

    def fake_subprocess_run(
        command: list[str],
        check: bool,
        capture_output: bool,
        text: bool,
    ) -> subprocess.CompletedProcess[str]:
        assert check is True
        assert capture_output is True
        assert text is True

        output_index = (
            command.index("--oc") + 1
        )

        output_directory = Path(
            command[output_index]
        )

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        result_path = (
            output_directory
            / "fimo.tsv"
        )

        result_path.write_text(
            "motif_id\tmotif_alt_id\t"
            "sequence_name\tstart\tstop\t"
            "strand\tscore\tp-value\t"
            "q-value\tmatched_sequence\n"
            "MA0595.1\tSREBF1\treference\t"
            "3\t8\t+\t13.8\t1.7e-05\t"
            "0.005\tCCGAAC\n"
            "MA1596.1\tZNF460\tmutated\t"
            "3\t8\t+\t12.9\t4.4e-06\t"
            "0.001\tCCGCAC\n",
            encoding="utf-8",
        )

        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout="",
            stderr="",
        )

    monkeypatch.setattr(
        "src.tf_binding.fimo_runner.subprocess.run",
        fake_subprocess_run,
    )

    hits = run_fimo(
        reference_sequence="AACCGAACTG",
        mutated_sequence="AACCGCACTG",
        motif_database_path=motif_path,
    )

    assert len(hits) == 2
    assert hits[0].motif_alt_id == "SREBF1"
    assert hits[1].motif_alt_id == "ZNF460"
