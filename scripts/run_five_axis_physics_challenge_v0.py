#!/usr/bin/env python3
from __future__ import annotations

"""Run the sequence-only five-axis physics challenge with an honest claim gate."""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_five_axis_physics import (  # noqa: E402
    FIVE_AXIS_PHYSICS_KIND,
    run_five_axis_challenge,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)

DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "five_axis_physics_challenge_v0"


def _write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _flatten_row_report(row: dict[str, object]) -> dict[str, object]:
    flattened = dict(row)
    boundary = flattened.pop("boundary", {})
    metric = flattened.pop("metric_after_native_audit", {})
    if isinstance(boundary, dict):
        for key, value in boundary.items():
            flattened[f"boundary_{key}"] = value
    if isinstance(metric, dict):
        for key, value in metric.items():
            flattened[f"metric_{key}"] = value
    return flattened


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--max-sequence-separation", type=int, default=160)
    args = parser.parse_args()

    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    packet = run_five_axis_challenge(
        rows,
        max_sequence_separation=args.max_sequence_separation,
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "five_axis_physics_challenge_report.json"
    rows_path = out_dir / "five_axis_physics_challenge_rows.csv"
    controls_path = out_dir / "five_axis_physics_challenge_controls.csv"
    decisions_path = out_dir / "five_axis_physics_challenge_selected_and_top_rejected_decisions.csv"
    certificate_path = out_dir / "five_axis_physics_challenge_certificate.json"

    payload = packet.to_dict()
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(rows_path, [_flatten_row_report(row.to_dict()) for row in packet.rows])
    _write_csv(controls_path, [control.to_dict() for control in packet.controls])
    _write_csv(decisions_path, [decision.to_dict() for decision in packet.decisions])
    certificate = {
        "kind": FIVE_AXIS_PHYSICS_KIND,
        "report": str(report_path),
        "rows": str(rows_path),
        "controls": str(controls_path),
        "decisions": str(decisions_path),
        "coordinate_truth_used_before_selection": packet.coordinate_truth_used_before_selection,
        "native_truth_used_before_selection": packet.native_truth_used_before_selection,
        "learned_prior_used_before_selection": packet.learned_prior_used_before_selection,
        "msa_used_before_selection": packet.msa_used_before_selection,
        "template_used_before_selection": packet.template_used_before_selection,
        "raw_sequence_exposed": packet.raw_sequence_exposed,
        "universal_physical_law_claim_allowed": packet.universal_physical_law_claim_allowed,
        "folding_problem_solved": packet.folding_problem_solved,
        "claim_rejection_reason": packet.claim_rejection_reason,
    }
    certificate_path.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")

    print(json.dumps({
        "kind": packet.kind,
        "row_count": packet.row_count,
        "mean_native_contact_precision_after_audit": packet.mean_native_contact_precision_after_audit,
        "mean_native_contact_recall_after_audit": packet.mean_native_contact_recall_after_audit,
        "mean_long_range_contact_recall_after_audit": packet.mean_long_range_contact_recall_after_audit,
        "mean_contact_map_f1_after_audit": packet.mean_contact_map_f1_after_audit,
        "mean_f1_margin_vs_best_control": packet.mean_f1_margin_vs_best_control,
        "universal_physical_law_claim_allowed": packet.universal_physical_law_claim_allowed,
        "folding_problem_solved": packet.folding_problem_solved,
        "claim_rejection_reason": packet.claim_rejection_reason,
        "report_path": str(report_path),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
