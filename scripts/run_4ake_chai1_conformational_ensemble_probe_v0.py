#!/usr/bin/env python3
"""Bounded 4AKE conformational-ensemble probe for Chai-1/Boltz-style outputs.

This is the path that actually matches the 4AKE failure mode: AdK has open and
closed conformations, so a single static contact predictor is the wrong source.
The script can consume already-generated Chai-1/Boltz PDB/mmCIF samples or run a
user-supplied local command with a hard timeout.  No GIFs are generated.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_chai1_conformational_ensemble import build_conformational_ensemble_packet  # noqa: E402
from pharmacotopology.folding_native_contact_eval import evaluate_contact_prediction  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows  # noqa: E402


DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_OUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "4ake_chai1_conformational_ensemble_v0"


def _row_by_accession(rows: Sequence[object], accession: str):
    matches = [row for row in rows if getattr(row, "source_accession") == accession]
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one row for {accession}, found {len(matches)}")
    return matches[0]


def _tail(text: str, limit: int = 2000) -> str:
    return text if len(text) <= limit else text[-limit:]


def _command_is_unsafe_for_no_msa_claim(command: str) -> bool:
    lowered = command.lower()
    return any(token in lowered for token in ("--use-msa", "use_msa", "msa-server", "mmseqs", "colabfold"))


def _run_prediction_command(
    *,
    command_template: str,
    sequence: str,
    output_dir: Path,
    timeout_seconds: float,
    allow_msa_command: bool,
) -> dict[str, object]:
    if _command_is_unsafe_for_no_msa_claim(command_template) and not allow_msa_command:
        return {
            "attempted": False,
            "status": "rejected_msa_like_prediction_command_for_no_msa_probe",
            "timed_out": False,
            "returncode": None,
            "stdout_tail": "",
            "stderr_tail": "",
            "raw_sequence_persisted": False,
        }
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pharm_chai1_4ake_") as tmp:
        fasta_path = Path(tmp) / "query.fasta"
        fasta_path.write_text(f">4AKE_A_sequence_only\n{sequence}\n", encoding="utf-8")
        command = command_template.format(fasta=str(fasta_path), out_dir=str(output_dir), output_dir=str(output_dir))
        try:
            completed = subprocess.run(
                shlex.split(command),
                cwd=str(REPO_ROOT),
                check=False,
                timeout=max(5.0, float(timeout_seconds)),
                capture_output=True,
                text=True,
                env={**os.environ, "PYTHONPATH": str(SRC_ROOT)},
            )
            return {
                "attempted": True,
                "status": "completed" if completed.returncode == 0 else "prediction_command_failed",
                "timed_out": False,
                "returncode": completed.returncode,
                "stdout_tail": _tail(completed.stdout or ""),
                "stderr_tail": _tail(completed.stderr or ""),
                "raw_sequence_persisted": False,
                "command_template": command_template,
            }
        except subprocess.TimeoutExpired as exc:
            return {
                "attempted": True,
                "status": "prediction_command_timed_out",
                "timed_out": True,
                "returncode": None,
                "stdout_tail": _tail((exc.stdout or "") if isinstance(exc.stdout, str) else ""),
                "stderr_tail": _tail((exc.stderr or "") if isinstance(exc.stderr, str) else ""),
                "raw_sequence_persisted": False,
                "command_template": command_template,
            }


def _write_contacts_csv(path: Path, contacts: Sequence[tuple[int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["i", "j", "sequence_separation", "selected"])
        writer.writeheader()
        for i, j in contacts:
            writer.writerow({"i": i, "j": j, "sequence_separation": j - i, "selected": True})


def main() -> None:
    parser = argparse.ArgumentParser(description="4AKE Chai-1/Boltz conformational ensemble probe")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--predicted-structure-dir", default=None)
    parser.add_argument("--prediction-command", default=None, help="Example: 'chai-lab fold {fasta} {out_dir}'")
    parser.add_argument("--predicted-source-id", default="chai1_single_sequence_conformational_ensemble")
    parser.add_argument("--predicted-chain", default="A")
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument("--allow-msa-command", action="store_true")
    parser.add_argument("--allow-alphafold-like-input", action="store_true")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = tuple(load_real_coordinate_visual_rows(Path(args.benchmark_file)))
    row = _row_by_accession(rows, args.source_accession)

    structure_dir = Path(args.predicted_structure_dir) if args.predicted_structure_dir else out_dir / "predicted_structures"
    command_report = {"attempted": False, "status": "not_requested", "raw_sequence_persisted": False}
    if args.prediction_command:
        command_report = _run_prediction_command(
            command_template=args.prediction_command,
            sequence=row.sequence,
            output_dir=structure_dir,
            timeout_seconds=args.timeout_seconds,
            allow_msa_command=args.allow_msa_command,
        )

    packet = build_conformational_ensemble_packet(
        structure_dir=structure_dir,
        source_id=args.predicted_source_id,
        chain_id=args.predicted_chain or None,
        reject_alphafold_like_inputs=not args.allow_alphafold_like_input,
        minimum_ca_count=max(100, int(row.sequence_length * 0.70)),
    )
    metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=packet.contacts)
    solved = (
        metric.native_contact_precision >= 0.70
        and metric.native_contact_recall >= 0.70
        and metric.long_range_contact_recall >= 0.70
    )
    contacts_csv = out_dir / "4ake_chai1_conformational_ensemble_selected_contacts.csv"
    _write_contacts_csv(contacts_csv, packet.contacts)
    report = {
        "kind": "4ake_chai1_conformational_ensemble_probe_v0",
        "source_accession": row.source_accession,
        "sequence_length": row.sequence_length,
        "sequence_hash": row.sequence_sha256,
        "raw_sequence_persisted": False,
        "gif_generation_used": False,
        "alphafold_used_as_evidence": False,
        "msa_used_as_evidence": bool(args.allow_msa_command),
        "prediction_command": command_report,
        "ensemble_packet": packet.to_report_dict(),
        "selected_contacts_csv": str(contacts_csv),
        "metric_after_native_audit": metric.to_dict(),
        "folding_problem_solved_after_native_audit": solved,
        "benchmark_claim_allowed_after_native_audit": solved,
        "claim_rejection_reason": "none" if solved else "missing_or_insufficient_chai1_boltz_conformational_ensemble_output",
    }
    report_path = out_dir / "4ake_chai1_conformational_ensemble_probe_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "report_path": str(report_path),
        "parsed_conformer_count": packet.parsed_conformer_count,
        "selected_contact_count": packet.selected_contact_count,
        "precision": metric.native_contact_precision,
        "recall": metric.native_contact_recall,
        "long_range_recall": metric.long_range_contact_recall,
        "folding_problem_solved_after_native_audit": solved,
        "claim_rejection_reason": report["claim_rejection_reason"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
