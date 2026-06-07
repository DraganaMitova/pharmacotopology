from __future__ import annotations

"""MSA-free single-sequence structure-source adapter.

This module does not implement a neural structure predictor.  It provides the
safe boundary that was missing in the project: a non-AlphaFold, MSA-free
structure source can be plugged in as an independent evidence layer without
letting the run hang, without persisting raw sequence by default, and without
pretending that sequence-only hand rules are an independent structure source.
"""

import shlex
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping

from pharmacotopology.folding_real_coordinate_visual_benchmark import RealCoordinateVisualRow


SINGLE_SEQUENCE_STRUCTURE_SOURCE_KIND = "single_sequence_structure_source_adapter_v0"
SINGLE_SEQUENCE_STRUCTURE_CONTACT_SOURCE_KIND = "single_sequence_structure_model_contacts_v0"
SINGLE_SEQUENCE_STRUCTURE_REPORT_KIND = "msa_free_single_sequence_structure_probe_v0"
ALPHAFOLD_SOURCE_TOKENS = ("alphafold", "alpha_fold", "afdb", "af-")


@dataclass(frozen=True)
class SingleSequenceStructureSourceRun:
    kind: str
    row_id: str
    source_accession: str
    sequence_hash: str
    predictor_source_id: str
    predicted_pdb_path: str
    output_pdb_exists: bool
    prediction_command_provided: bool
    prediction_command_executed: bool
    prediction_command_returncode: int | None
    prediction_command_timed_out: bool
    timeout_seconds: float
    query_fasta_persisted: bool
    raw_sequence_exposed_in_persisted_artifacts: bool
    alphafold_source_allowed: bool
    alphafold_like_source_id: bool
    source_rejected: bool
    status: str
    stdout_tail: str = ""
    stderr_tail: str = ""

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def is_alphafold_like_source_id(source_id: str) -> bool:
    normalized = source_id.lower().replace("/", "_").replace("-", "-")
    return any(token in normalized for token in ALPHAFOLD_SOURCE_TOKENS)


def fasta_text_for_row(row: RealCoordinateVisualRow) -> str:
    header = f">{row.source_accession}|{row.row_id}|sha256:{row.sequence_sha256}\n"
    chunks = [row.sequence[index : index + 80] for index in range(0, len(row.sequence), 80)]
    return header + "\n".join(chunks) + "\n"


def _tail(text: str, *, limit: int = 1200) -> str:
    if len(text) <= limit:
        return text
    return text[-limit:]


def _format_command(command_template: str, replacements: Mapping[str, str]) -> list[str]:
    try:
        rendered = command_template.format(**replacements)
    except KeyError as exc:
        known = ", ".join(sorted(replacements))
        raise ValueError(f"unknown prediction command placeholder {exc.args[0]!r}; known: {known}") from exc
    return shlex.split(rendered)


def run_single_sequence_prediction_source(
    *,
    row: RealCoordinateVisualRow,
    predicted_pdb_path: Path,
    predictor_source_id: str,
    prediction_command: str | None = None,
    timeout_seconds: float = 120.0,
    keep_query_fasta: bool = False,
    allow_alphafold_source: bool = False,
    working_directory: Path | None = None,
) -> SingleSequenceStructureSourceRun:
    """Optionally run an MSA-free predictor command and return a safety manifest.

    The command is opt-in and bounded by ``timeout_seconds``.  It receives only a
    FASTA file and an output PDB path through placeholders:

    ``{fasta}``, ``{output_pdb}``, ``{source_accession}``, ``{row_id}``,
    ``{sequence_length}``.
    """

    predicted_pdb_path = predicted_pdb_path.resolve()
    predicted_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    timeout = max(1.0, float(timeout_seconds))
    alphafold_like = is_alphafold_like_source_id(
        f"{predictor_source_id} {predicted_pdb_path.name} {predicted_pdb_path.parent.name}"
    )
    if alphafold_like and not allow_alphafold_source:
        return SingleSequenceStructureSourceRun(
            kind=SINGLE_SEQUENCE_STRUCTURE_SOURCE_KIND,
            row_id=row.row_id,
            source_accession=row.source_accession,
            sequence_hash=row.sequence_sha256,
            predictor_source_id=predictor_source_id,
            predicted_pdb_path=str(predicted_pdb_path),
            output_pdb_exists=predicted_pdb_path.exists(),
            prediction_command_provided=prediction_command is not None,
            prediction_command_executed=False,
            prediction_command_returncode=None,
            prediction_command_timed_out=False,
            timeout_seconds=timeout,
            query_fasta_persisted=False,
            raw_sequence_exposed_in_persisted_artifacts=False,
            alphafold_source_allowed=allow_alphafold_source,
            alphafold_like_source_id=True,
            source_rejected=True,
            status="rejected_alphafold_like_source_id_for_msa_free_probe",
        )

    stdout_tail = ""
    stderr_tail = ""
    returncode: int | None = None
    timed_out = False
    executed = False
    query_fasta_persisted = False
    raw_sequence_persisted = False

    if prediction_command:
        executed = True
        if keep_query_fasta:
            query_fasta_path = (working_directory or predicted_pdb_path.parent) / "query.fasta"
            query_fasta_path.parent.mkdir(parents=True, exist_ok=True)
            query_fasta_path.write_text(fasta_text_for_row(row), encoding="utf-8")
            query_fasta_persisted = True
            raw_sequence_persisted = True
            cleanup_context = None
        else:
            cleanup_context = tempfile.TemporaryDirectory(prefix="pharm_msa_free_")
            query_fasta_path = Path(cleanup_context.name) / "query.fasta"
            query_fasta_path.write_text(fasta_text_for_row(row), encoding="utf-8")
        replacements = {
            "fasta": str(query_fasta_path),
            "output_pdb": str(predicted_pdb_path),
            "source_accession": row.source_accession,
            "row_id": row.row_id,
            "sequence_length": str(row.sequence_length),
        }
        try:
            command = _format_command(prediction_command, replacements)
            completed = subprocess.run(
                command,
                cwd=str(working_directory) if working_directory is not None else None,
                check=False,
                timeout=timeout,
                capture_output=True,
                text=True,
            )
            returncode = completed.returncode
            stdout_tail = _tail(completed.stdout or "")
            stderr_tail = _tail(completed.stderr or "")
        except subprocess.TimeoutExpired as exc:
            timed_out = True
            stdout_tail = _tail((exc.stdout or "") if isinstance(exc.stdout, str) else "")
            stderr_tail = _tail((exc.stderr or "") if isinstance(exc.stderr, str) else "")
        finally:
            if cleanup_context is not None:
                cleanup_context.cleanup()

    exists = predicted_pdb_path.exists() and predicted_pdb_path.stat().st_size > 0
    if timed_out:
        status = "prediction_command_timed_out"
    elif prediction_command and returncode not in (0, None):
        status = "prediction_command_failed"
    elif not exists:
        status = "missing_single_sequence_predictor_output"
    elif raw_sequence_persisted:
        status = "prediction_available_but_raw_sequence_persisted"
    else:
        status = "prediction_available"

    return SingleSequenceStructureSourceRun(
        kind=SINGLE_SEQUENCE_STRUCTURE_SOURCE_KIND,
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_hash=row.sequence_sha256,
        predictor_source_id=predictor_source_id,
        predicted_pdb_path=str(predicted_pdb_path),
        output_pdb_exists=exists,
        prediction_command_provided=prediction_command is not None,
        prediction_command_executed=executed,
        prediction_command_returncode=returncode,
        prediction_command_timed_out=timed_out,
        timeout_seconds=timeout,
        query_fasta_persisted=query_fasta_persisted,
        raw_sequence_exposed_in_persisted_artifacts=raw_sequence_persisted,
        alphafold_source_allowed=allow_alphafold_source,
        alphafold_like_source_id=alphafold_like,
        source_rejected=False,
        status=status,
        stdout_tail=stdout_tail,
        stderr_tail=stderr_tail,
    )
