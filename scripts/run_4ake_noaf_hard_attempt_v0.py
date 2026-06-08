#!/usr/bin/env python3
"""One-command hard attempt for 4AKE without AlphaFold as evidence.

This orchestrator tries the internet MSA-free ESMFold path first.  It then runs a
locked AlphaFold fixture only as a positive-control safety test: the direct
metrics should be high, but the claim must stay disallowed because the source is
AlphaFold-like.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
DEFAULT_OUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "4ake_noaf_hard_attempt_v0"
AF_FIXTURE = REPO_ROOT / "data" / "independent_contact_sources" / "AF-P69441-F1-model_v4.pdb"


def _tail(text: str, *, limit: int = 2400) -> str:
    return text if len(text) <= limit else text[-limit:]


def _run(command: list[str], *, timeout_seconds: float, cwd: Path = REPO_ROOT) -> dict[str, object]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    env.setdefault("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")
    try:
        completed = subprocess.run(
            command,
            cwd=str(cwd),
            check=False,
            timeout=max(1.0, float(timeout_seconds)),
            capture_output=True,
            text=True,
            env=env,
        )
        return {
            "command": command,
            "returncode": completed.returncode,
            "timed_out": False,
            "stdout_tail": _tail(completed.stdout or ""),
            "stderr_tail": _tail(completed.stderr or ""),
        }
    except subprocess.TimeoutExpired as exc:
        return {
            "command": command,
            "returncode": None,
            "timed_out": True,
            "stdout_tail": _tail((exc.stdout or "") if isinstance(exc.stdout, str) else ""),
            "stderr_tail": _tail((exc.stderr or "") if isinstance(exc.stderr, str) else ""),
        }


def _load_json(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Hard bounded 4AKE no-AlphaFold attempt.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--internet-timeout-seconds", type=float, default=60.0)
    parser.add_argument("--probe-timeout-seconds", type=float, default=90.0)
    parser.add_argument("--test-timeout-seconds", type=float, default=120.0)
    parser.add_argument("--try-biolm", action="store_true")
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument("--skip-positive-control", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    internet_dir = out_dir / "internet_msa_free_sweep"
    internet_cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "run_4ake_msa_free_internet_sweep_v0.py"),
        "--out-dir",
        str(internet_dir),
        "--timeout-seconds",
        str(args.internet_timeout_seconds),
        "--include-sequence-physical-priors",
    ]
    if args.try_biolm:
        internet_cmd.append("--try-biolm")
    internet_run = _run(internet_cmd, timeout_seconds=args.internet_timeout_seconds + args.probe_timeout_seconds + 30)
    internet_summary = _load_json(internet_dir / "4ake_msa_free_internet_sweep_summary.json")

    positive_control_run = None
    positive_control_report = None
    if not args.skip_positive_control and AF_FIXTURE.exists():
        positive_dir = out_dir / "alphafold_positive_control_not_claimed"
        positive_cmd = [
            sys.executable,
            str(REPO_ROOT / "scripts" / "run_msa_free_single_sequence_structure_probe_v0.py"),
            "--source-accession",
            "4AKE:A",
            "--predicted-pdb",
            str(AF_FIXTURE),
            "--predicted-source-id",
            "alphafold_positive_control_not_claimed",
            "--allow-alphafold-source",
            "--predicted-pdb-chain",
            "A",
            "--out-dir",
            str(positive_dir),
            "--include-sequence-physical-priors",
        ]
        positive_control_run = _run(positive_cmd, timeout_seconds=args.probe_timeout_seconds)
        positive_control_report = _load_json(positive_dir / "msa_free_single_sequence_structure_probe_report.json")

    test_run = None
    if not args.skip_tests:
        test_run = _run(
            [sys.executable, "-m", "pytest", "-q", "-vv"],
            timeout_seconds=args.test_timeout_seconds,
        )

    folding_problem_solved = False
    if internet_summary and internet_summary.get("folding_problem_solved") is True:
        folding_problem_solved = True

    positive_control_claim_blocked = None
    positive_control_direct_metric = None
    if isinstance(positive_control_report, dict):
        safety = positive_control_report.get("safety", {})
        positive_control_claim_blocked = (
            isinstance(safety, dict)
            and safety.get("alphafold_used_by_this_script") is True
            and safety.get("folding_problem_solved") is False
            and safety.get("benchmark_claim_allowed_by_direct_structure") is False
        )
        positive_control_direct_metric = positive_control_report.get("direct_structure_metric")

    summary = {
        "kind": "4ake_noaf_hard_attempt_v0",
        "gif_generation_used": False,
        "alphafold_used_as_claim_evidence": False,
        "internet_msa_free_sweep": internet_summary,
        "internet_msa_free_sweep_run": internet_run,
        "alphafold_positive_control_not_claimed_run": positive_control_run,
        "alphafold_positive_control_claim_blocked": positive_control_claim_blocked,
        "alphafold_positive_control_direct_metric": positive_control_direct_metric,
        "fast_pytest_run": test_run,
        "folding_problem_solved": folding_problem_solved,
    }
    summary_path = out_dir / "4ake_noaf_hard_attempt_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    md = [
        "# 4AKE no-AlphaFold hard attempt v0",
        "",
        f"folding_problem_solved: `{folding_problem_solved}`",
        "gif_generation_used: `False`",
        "alphafold_used_as_claim_evidence: `False`",
        "",
        "## Internet MSA-free sweep",
        f"status: `{[a.get('status') for a in (internet_summary or {}).get('attempts', [])] if internet_summary else 'missing_summary'}`",
        f"successful_prediction_count: `{sum(1 for a in (internet_summary or {}).get('attempts', []) if a.get('success')) if internet_summary else 0}`",
        "",
        "## AlphaFold positive control (blocked as claim)",
        f"claim_blocked: `{positive_control_claim_blocked}`",
        f"direct_metric: `{positive_control_direct_metric}`",
        "",
        "## Test run",
        f"returncode: `{test_run.get('returncode') if test_run else None}`",
        f"timed_out: `{test_run.get('timed_out') if test_run else None}`",
        "",
        "This artifact treats the AlphaFold fixture only as a safety positive-control. It is not used as no-AlphaFold evidence.",
    ]
    (out_dir / "4ake_noaf_hard_attempt_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps({"summary_path": str(summary_path), "folding_problem_solved": folding_problem_solved}, indent=2))


if __name__ == "__main__":
    main()
