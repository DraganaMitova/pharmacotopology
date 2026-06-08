#!/usr/bin/env python3
"""Run the locked 4AKE AlphaFold independent-source ensemble probe.

This script is intentionally only an orchestration wrapper around
run_independent_contact_ensemble_probe_v0.py.  It does not read native/contact
truth before selection; truth is attached only inside the probe evaluator after
selected contacts are produced.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PDB_PATH = REPO_ROOT / "data" / "independent_contact_sources" / "AF-P69441-F1-model_v4.pdb"
RUNNER = REPO_ROOT / "scripts" / "run_independent_contact_ensemble_probe_v0.py"
OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
SOURCE_ID = "alphafold_db_AF-P69441-F1-model_v4"


PROBES = (
    {
        "label": "aggressive_relative_frontier_min2",
        "event_source": "aggressive_relative_frontier",
        "min_votes": "2",
        "prefix": "independent_contact_ensemble_4ake_alphafold",
    },
    {
        "label": "aggressive_relative_frontier_min3",
        "event_source": "aggressive_relative_frontier",
        "min_votes": "3",
        "prefix": "independent_contact_ensemble_4ake_alphafold_min3",
    },
    {
        "label": "competitive_pool_min2",
        "event_source": "competitive_pool",
        "min_votes": "2",
        "prefix": "independent_contact_ensemble_4ake_alphafold_competitive_pool",
    },
    {
        "label": "candidate_pool_min2",
        "event_source": "candidate_pool",
        "min_votes": "2",
        "prefix": "independent_contact_ensemble_4ake_alphafold_candidate_pool",
    },
)


def _run_probe(probe: dict[str, str], *, timeout_seconds: int, fail_on_probe_error: bool) -> dict[str, object]:
    if not PDB_PATH.exists():
        raise FileNotFoundError(f"missing AlphaFold PDB source: {PDB_PATH}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    report_path = OUTPUT_DIR / f"{probe['prefix']}_v0.json"
    command = [
        sys.executable,
        str(RUNNER),
        "--predicted-pdb",
        str(PDB_PATH.relative_to(REPO_ROOT)),
        "--predicted-pdb-chain",
        "A",
        "--predicted-source-id",
        SOURCE_ID,
        "--event-source",
        probe["event_source"],
        "--min-votes",
        probe["min_votes"],
        "--report-output",
        str(report_path.relative_to(REPO_ROOT)),
        "--decisions-output",
        str((OUTPUT_DIR / f"{probe['prefix']}_decisions_v0.csv").relative_to(REPO_ROOT)),
        "--selected-output",
        str((OUTPUT_DIR / f"{probe['prefix']}_selected_contacts_v0.csv").relative_to(REPO_ROOT)),
        "--evidence-output",
        str((OUTPUT_DIR / f"{probe['prefix']}_evidence_v0.json").relative_to(REPO_ROOT)),
    ]
    try:
        subprocess.run(
            command,
            check=True,
            cwd=REPO_ROOT,
            timeout=timeout_seconds,
            capture_output=True,
            text=True,
        )
    except subprocess.TimeoutExpired as exc:
        if fail_on_probe_error:
            raise
        return {
            "label": probe["label"],
            "probe_status": "timeout",
            "timeout_seconds": timeout_seconds,
            "report_path": str(report_path.relative_to(REPO_ROOT)),
            "stdout_tail": (exc.stdout or "").splitlines()[-4:] if isinstance(exc.stdout, str) else [],
            "stderr_tail": (exc.stderr or "").splitlines()[-4:] if isinstance(exc.stderr, str) else [],
        }
    except subprocess.CalledProcessError as exc:
        if fail_on_probe_error:
            raise
        return {
            "label": probe["label"],
            "probe_status": "failed",
            "returncode": exc.returncode,
            "report_path": str(report_path.relative_to(REPO_ROOT)),
            "stdout_tail": (exc.stdout or "").splitlines()[-4:],
            "stderr_tail": (exc.stderr or "").splitlines()[-4:],
        }
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    report = payload["report"]
    return {
        "label": probe["label"],
        "probe_status": "ok",
        "report_path": str(report_path.relative_to(REPO_ROOT)),
        "benchmark_claim_allowed": report["benchmark_claim_allowed"],
        "claim_rejection_reason": report["claim_rejection_reason"],
        "final_pair_count": report["final_pair_count"],
        "final_long_range_pair_count": report["final_long_range_pair_count"],
        "long_range_precision": report["long_range_precision"],
        "long_range_recall": report["long_range_recall"],
        "true_positive_long_range_contacts": report["true_positive_long_range_contacts"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the locked 4AKE AlphaFold independent-source probe suite.")
    parser.add_argument("--timeout-seconds", type=int, default=30)
    parser.add_argument(
        "--suite-mode",
        choices=("fast", "full"),
        default="fast",
        help="fast runs the primary aggressive frontier probe; full runs all historical probes.",
    )
    parser.add_argument("--fail-on-probe-error", action="store_true")
    args = parser.parse_args()

    selected_probes = PROBES if args.suite_mode == "full" else PROBES[:1]
    summaries = [
        _run_probe(
            probe,
            timeout_seconds=args.timeout_seconds,
            fail_on_probe_error=args.fail_on_probe_error,
        )
        for probe in selected_probes
    ]
    successful = [item for item in summaries if item.get("probe_status") == "ok"]
    summary_path = OUTPUT_DIR / "independent_contact_ensemble_4ake_alphafold_real_source_summary_v0.json"
    summary = {
        "kind": "real_alphafold_independent_source_probe_summary_v0",
        "suite_mode": args.suite_mode,
        "timeout_seconds_per_probe": args.timeout_seconds,
        "source_id": SOURCE_ID,
        "source_pdb": str(PDB_PATH.relative_to(REPO_ROOT)),
        "source_is_native_coordinate_leakage": False,
        "native_truth_attached_after_selection_for_evaluation": True,
        "probes": summaries,
        "best_recall_safe_probe": (
            max(
                successful,
                key=lambda item: (item["long_range_recall"], item["long_range_precision"]),
            )
            if successful
            else None
        ),
        "highest_precision_safe_probe": (
            max(
                successful,
                key=lambda item: (item["long_range_precision"], item["long_range_recall"]),
            )
            if successful
            else None
        ),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(summary_path)


if __name__ == "__main__":
    main()
