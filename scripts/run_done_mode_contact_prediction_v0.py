#!/usr/bin/env python3
"""One-command DONE-mode contact prediction wrapper.

DONE-mode means every target can be processed through the same interface, but a
benchmark/scientific claim is only allowed when an independent source is present
and the leakage gates stay closed.  Without independent evidence the runner
returns a structured abstention instead of pretending to solve the target.
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PROBE = REPO_ROOT / "scripts" / "run_independent_contact_ensemble_probe_v0.py"
VISUAL_PROOF = REPO_ROOT / "scripts" / "render_4ake_contact_visual_proof_gif_v0.py"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "done_mode_v0"


def _slug(source_accession: str) -> str:
    return source_accession.lower().replace(":", "_").replace("/", "_")


def _safe_rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def _kill_process_group(process: subprocess.Popen[str]) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    except PermissionError:
        process.terminate()
    deadline = time.monotonic() + 2.0
    while process.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)
    if process.poll() is None:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            return
        except PermissionError:
            process.kill()


def _run(command: list[str], *, timeout: int) -> subprocess.CompletedProcess[str]:
    process = subprocess.Popen(
        command,
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        _kill_process_group(process)
        stdout, stderr = process.communicate()
        raise subprocess.TimeoutExpired(
            cmd=command,
            timeout=timeout,
            output=stdout or exc.output,
            stderr=stderr or exc.stderr,
        )
    completed = subprocess.CompletedProcess(
        command,
        process.returncode,
        stdout,
        stderr,
    )
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(
            completed.returncode,
            command,
            output=stdout,
            stderr=stderr,
        )
    return completed


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the honest DONE-mode predictor for any locked protein target. "
            "A claim is allowed only when an independent source exists and leakage gates pass."
        )
    )
    parser.add_argument("--source-accession", required=True, help="Locked benchmark target, e.g. 4AKE:A")
    parser.add_argument("--predicted-pdb", default=None, help="Independent predicted/AF/RoseTTAFold PDB")
    parser.add_argument("--predicted-pdb-chain", default=None)
    parser.add_argument("--predicted-source-id", default="external_predicted_structure")
    parser.add_argument("--independent-contact-json", action="append", default=[])
    parser.add_argument(
        "--event-source",
        default="aggressive_relative_frontier",
        choices=(
            "candidate_pool",
            "competitive_pool",
            "aggressive_relative_frontier",
            "self_deciding_generated_frontier",
        ),
    )
    parser.add_argument("--min-votes", type=int, default=2)
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--timeout-seconds", type=int, default=60)
    parser.add_argument(
        "--visual-proof-timeout-seconds",
        type=int,
        default=None,
        help="Optional separate timeout for GIF rendering. Defaults to --timeout-seconds.",
    )
    parser.add_argument(
        "--render-visual-proof",
        action="store_true",
        help="Render the 4AKE visual proof GIF when source-accession is 4AKE:A.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir) if args.out_dir else DEFAULT_OUTPUT_ROOT / _slug(args.source_accession)
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "contact_prediction_report.json"
    decisions_path = out_dir / "contact_prediction_decisions.csv"
    selected_path = out_dir / "contact_prediction_selected_contacts.csv"
    evidence_path = out_dir / "contact_prediction_evidence.json"

    command = [
        sys.executable,
        str(PROBE),
        "--source-accession",
        args.source_accession,
        "--event-source",
        args.event_source,
        "--min-votes",
        str(args.min_votes),
        "--report-output",
        str(report_path),
        "--decisions-output",
        str(decisions_path),
        "--selected-output",
        str(selected_path),
        "--evidence-output",
        str(evidence_path),
    ]
    for item in args.independent_contact_json:
        command.extend(["--independent-contact-json", item])
    if args.predicted_pdb:
        command.extend(["--predicted-pdb", args.predicted_pdb])
        if args.predicted_pdb_chain:
            command.extend(["--predicted-pdb-chain", args.predicted_pdb_chain])
        command.extend(["--predicted-source-id", args.predicted_source_id])

    result = _run(command, timeout=args.timeout_seconds)
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    report = payload["report"]

    manifest: dict[str, object] = {
        "kind": "done_mode_contact_prediction_v0",
        "source_accession": args.source_accession,
        "done_contract": {
            "predictor_runs_for_every_locked_target": True,
            "claim_requires_independent_evidence": True,
            "abstains_instead_of_guessing_without_independent_evidence": True,
            "native_or_coordinate_truth_before_selection_allowed": False,
            "raw_sequence_exposure_allowed": False,
        },
        "status": "claim_allowed" if report["benchmark_claim_allowed"] else "abstained_or_rejected",
        "claim_rejection_reason": report["claim_rejection_reason"],
        "metrics": {
            "long_range_precision": report["long_range_precision"],
            "long_range_recall": report["long_range_recall"],
            "contact_precision": report["contact_precision"],
            "contact_recall": report["contact_recall"],
        },
        "outputs": {
            "report": _safe_rel(report_path),
            "decisions": _safe_rel(decisions_path),
            "selected_contacts": _safe_rel(selected_path),
            "evidence": _safe_rel(evidence_path),
        },
        "visual_proof_enabled": False,
        "visual_proof_policy": "off_by_default_explicit_render_visual_proof_flag_required",
        "probe_stdout_tail": result.stdout.splitlines()[-4:],
    }

    if args.render_visual_proof and args.source_accession == "4AKE:A":
        gif_path = out_dir / "4ake_visual_proof.gif"
        gif_manifest_path = out_dir / "4ake_visual_proof_manifest.json"
        frames_dir = out_dir / "frames"
        _run(
            [
                sys.executable,
                str(VISUAL_PROOF),
                "--report",
                str(report_path),
                "--evidence",
                str(evidence_path),
                "--selected",
                str(selected_path),
                "--output-gif",
                str(gif_path),
                "--output-manifest",
                str(gif_manifest_path),
                "--output-frames-dir",
                str(frames_dir),
            ],
            timeout=args.visual_proof_timeout_seconds or args.timeout_seconds,
        )
        manifest["visual_proof_enabled"] = True
        manifest["outputs"] = dict(manifest["outputs"], visual_proof_gif=_safe_rel(gif_path), visual_proof_manifest=_safe_rel(gif_manifest_path))
    elif args.render_visual_proof:
        manifest["visual_proof_policy"] = "render_visual_proof_requested_but_only_4ake_a_has_visual_renderer"

    manifest_path = out_dir / "done_mode_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    print(manifest_path)


if __name__ == "__main__":
    main()
