#!/usr/bin/env python3
"""Try internet MSA-free predictors for 4AKE and immediately run the probe.

Default behavior is deliberately bounded:
- no AlphaFold endpoint is called;
- no raw FASTA/sequence file is persisted;
- each internet request has a timeout;
- each probe subprocess has a timeout;
- GIF generation is never invoked.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_msa_free_internet_predictors import (  # noqa: E402
    post_biolm_esmfold,
    post_esm_atlas_fold_sequence,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_OUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "4ake_msa_free_internet_sweep_v0"


def _row_by_accession(rows: Sequence[RealCoordinateVisualRow], accession: str) -> RealCoordinateVisualRow:
    matches = [row for row in rows if row.source_accession == accession]
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one row for {accession}, found {len(matches)}")
    return matches[0]


def _tail(text: str, *, limit: int = 1600) -> str:
    return text if len(text) <= limit else text[-limit:]


def _run_probe(
    *,
    pdb_path: Path,
    predictor_source_id: str,
    chain_id: str | None,
    out_dir: Path,
    source_accession: str,
    timeout_seconds: float,
    include_physical_priors: bool,
) -> dict[str, object]:
    probe_dir = out_dir / f"probe_{predictor_source_id}"
    command = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "run_msa_free_single_sequence_structure_probe_v0.py"),
        "--source-accession",
        source_accession,
        "--predicted-pdb",
        str(pdb_path),
        "--predicted-source-id",
        predictor_source_id,
        "--out-dir",
        str(probe_dir),
        "--timeout-seconds",
        str(max(1.0, timeout_seconds)),
    ]
    if chain_id:
        command.extend(["--predicted-pdb-chain", chain_id])
    if include_physical_priors:
        command.append("--include-sequence-physical-priors")
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    try:
        completed = subprocess.run(
            command,
            cwd=str(REPO_ROOT),
            check=False,
            timeout=max(5.0, timeout_seconds + 20.0),
            capture_output=True,
            text=True,
            env=env,
        )
        report_path = probe_dir / "msa_free_single_sequence_structure_probe_report.json"
        parsed_report: Mapping[str, object] | None = None
        if report_path.exists():
            parsed_report = json.loads(report_path.read_text(encoding="utf-8"))
        return {
            "command": command,
            "returncode": completed.returncode,
            "timed_out": False,
            "stdout_tail": _tail(completed.stdout or ""),
            "stderr_tail": _tail(completed.stderr or ""),
            "report_path": str(report_path),
            "report": parsed_report,
        }
    except subprocess.TimeoutExpired as exc:
        return {
            "command": command,
            "returncode": None,
            "timed_out": True,
            "stdout_tail": _tail((exc.stdout or "") if isinstance(exc.stdout, str) else ""),
            "stderr_tail": _tail((exc.stderr or "") if isinstance(exc.stderr, str) else ""),
            "report_path": str(probe_dir / "msa_free_single_sequence_structure_probe_report.json"),
            "report": None,
        }


def main() -> None:
    parser = argparse.ArgumentParser(description="Bounded 4AKE MSA-free internet predictor sweep.")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--timeout-seconds", type=float, default=120.0)
    parser.add_argument("--predicted-pdb-chain", default="A")
    parser.add_argument("--skip-esm-atlas", action="store_true")
    parser.add_argument("--try-biolm", action="store_true")
    parser.add_argument("--biolm-api-key-env", default="BIOLM_API_KEY")
    parser.add_argument("--include-sequence-physical-priors", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    row = _row_by_accession(tuple(load_real_coordinate_visual_rows(Path(args.benchmark_file))), args.source_accession)

    attempts = []
    if not args.skip_esm_atlas:
        attempts.append(
            post_esm_atlas_fold_sequence(
                sequence=row.sequence,
                output_pdb_path=out_dir / "4ake_esm_atlas_esmfold.pdb",
                timeout_seconds=args.timeout_seconds,
            )
        )
    if args.try_biolm:
        attempts.append(
            post_biolm_esmfold(
                sequence=row.sequence,
                output_pdb_path=out_dir / "4ake_biolm_esmfold.pdb",
                api_key=os.environ.get(args.biolm_api_key_env),
                timeout_seconds=args.timeout_seconds,
            )
        )

    probe_runs = []
    for attempt in attempts:
        if not attempt.success:
            continue
        probe_runs.append(
            _run_probe(
                pdb_path=Path(attempt.output_pdb_path),
                predictor_source_id=attempt.predictor_id,
                chain_id=args.predicted_pdb_chain,
                out_dir=out_dir,
                source_accession=args.source_accession,
                timeout_seconds=args.timeout_seconds,
                include_physical_priors=args.include_sequence_physical_priors,
            )
        )

    solved_reports = []
    for run in probe_runs:
        report = run.get("report")
        if not isinstance(report, Mapping):
            continue
        safety = report.get("safety")
        if isinstance(safety, Mapping) and safety.get("folding_problem_solved") is True:
            solved_reports.append(str(run.get("report_path")))

    summary = {
        "kind": "4ake_msa_free_internet_predictor_sweep_v0",
        "source_accession": row.source_accession,
        "sequence_hash": row.sequence_sha256,
        "sequence_length": row.sequence_length,
        "raw_sequence_persisted": False,
        "alphafold_endpoint_used": False,
        "gif_generation_used": False,
        "attempts": [attempt.to_dict() for attempt in attempts],
        "probe_runs": probe_runs,
        "folding_problem_solved": bool(solved_reports),
        "solved_report_paths": solved_reports,
    }
    summary_path = out_dir / "4ake_msa_free_internet_sweep_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "summary_path": str(summary_path),
        "attempt_statuses": [attempt.status for attempt in attempts],
        "successful_prediction_count": sum(1 for attempt in attempts if attempt.success),
        "folding_problem_solved": bool(solved_reports),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
