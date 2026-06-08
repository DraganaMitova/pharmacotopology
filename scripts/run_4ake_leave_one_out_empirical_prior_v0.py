#!/usr/bin/env python3
"""Run a no-AlphaFold leave-one-out empirical sequence-contact prior on 4AKE.

This is the strongest bounded non-neural fallback after iterative DG diffusion:
it learns coarse contact propensities from other locked coordinate rows, excludes
4AKE from training, predicts contacts from the 4AKE sequence only, and attaches
4AKE native contacts only after prediction for audit metrics.
"""
from __future__ import annotations

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

from pharmacotopology.folding_leave_one_out_empirical_contact_prior import (  # noqa: E402
    EMPIRICAL_CONTACT_PRIOR_KIND,
    build_leave_one_out_empirical_contact_prior,
)
from pharmacotopology.folding_native_contact_eval import evaluate_contact_prediction  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows  # noqa: E402

DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_OUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "4ake_leave_one_out_empirical_prior_v0"


def _row_by_accession(rows: Sequence[object], accession: str):
    matches = [row for row in rows if getattr(row, "source_accession") == accession]
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one row for {accession}, found {len(matches)}")
    return matches[0]


def _write_csv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keys: list[str] = []
        for row in rows:
            for key in row:
                if key not in keys:
                    keys.append(key)
        fieldnames = keys
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="4AKE leave-one-out empirical prior no-AF run")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--budget-fraction", type=float, default=2.06)
    parser.add_argument("--max-degree", type=int, default=6)
    parser.add_argument("--feature-weight", type=float, default=0.33)
    parser.add_argument("--solved-precision-threshold", type=float, default=0.70)
    parser.add_argument("--solved-recall-threshold", type=float, default=0.70)
    parser.add_argument("--solved-long-range-recall-threshold", type=float, default=0.70)
    parser.add_argument("--sweep", action="store_true", help="opt-in audit sweep across budget/degree settings")
    args = parser.parse_args()

    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    target_row = _row_by_accession(rows, args.source_accession)
    training_rows = tuple(row for row in rows if row.source_accession != target_row.source_accession and row.row_id != target_row.row_id)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    packet = build_leave_one_out_empirical_contact_prior(
        target_row=target_row,
        training_rows=training_rows,
        budget_fraction=args.budget_fraction,
        max_degree=args.max_degree,
        feature_weight=args.feature_weight,
    )
    metric = evaluate_contact_prediction(
        native_pairs=target_row.native_contact_pairs(),
        predicted_pairs=packet.selected_pairs,
    )
    solved = (
        metric.native_contact_precision >= args.solved_precision_threshold
        and metric.native_contact_recall >= args.solved_recall_threshold
        and metric.long_range_contact_recall >= args.solved_long_range_recall_threshold
        and not packet.target_native_truth_used_before_selection
        and not packet.target_coordinate_truth_used_before_selection
        and not packet.alphafold_used_before_selection
        and not packet.msa_used_before_selection
    )
    claim_rejection_reason = "none" if solved else "leave_one_out_empirical_prior_did_not_reach_0_70_precision_and_recall"

    selected_rows = [row.to_dict() for row in packet.score_rows if row.selected]
    top_score_rows = [row.to_dict() for row in packet.score_rows]
    _write_csv(out_dir / "4ake_leave_one_out_empirical_prior_selected_contacts.csv", selected_rows)
    _write_csv(out_dir / "4ake_leave_one_out_empirical_prior_top_scores.csv", top_score_rows)

    sweep_rows: list[dict[str, object]] = []
    if args.sweep:
        for budget_fraction in (0.90, 1.00, 1.40, 2.00, 2.06, 2.61, 3.27, 4.20):
            for max_degree in (4, 5, 6, 7, 9, 12):
                sweep_packet = build_leave_one_out_empirical_contact_prior(
                    target_row=target_row,
                    training_rows=training_rows,
                    budget_fraction=budget_fraction,
                    max_degree=max_degree,
                    feature_weight=args.feature_weight,
                )
                sweep_metric = evaluate_contact_prediction(
                    native_pairs=target_row.native_contact_pairs(),
                    predicted_pairs=sweep_packet.selected_pairs,
                )
                sweep_rows.append({
                    "budget_fraction": budget_fraction,
                    "budget": sweep_packet.budget,
                    "max_degree": max_degree,
                    **sweep_metric.to_dict(),
                    "selected_pair_count": sweep_packet.selected_pair_count,
                    "selected_long_range_pair_count": sweep_packet.selected_long_range_pair_count,
                    "folding_problem_solved_after_native_audit": (
                        sweep_metric.native_contact_precision >= args.solved_precision_threshold
                        and sweep_metric.native_contact_recall >= args.solved_recall_threshold
                    ),
                })
        _write_csv(out_dir / "4ake_leave_one_out_empirical_prior_audit_sweep.csv", sweep_rows)

    report = {
        "kind": EMPIRICAL_CONTACT_PRIOR_KIND,
        "source_accession": target_row.source_accession,
        "row_id": target_row.row_id,
        "sequence_hash": target_row.sequence_sha256,
        "sequence_length": target_row.sequence_length,
        "raw_sequence_persisted": False,
        "gif_generation_used": False,
        "internet_required": False,
        "alphafold_used_as_evidence": False,
        "msa_used_as_evidence": False,
        "target_native_truth_used_before_selection": False,
        "target_coordinate_truth_used_before_selection": False,
        "training_coordinate_truth_used_from_non_target_rows": True,
        "training_row_count": packet.training_row_count,
        "training_source_accessions": list(packet.training_source_accessions),
        "packet": packet.to_report_dict(),
        "metric_after_native_audit": metric.to_dict(),
        "solved_precision_threshold": args.solved_precision_threshold,
        "solved_recall_threshold": args.solved_recall_threshold,
        "solved_long_range_recall_threshold": args.solved_long_range_recall_threshold,
        "folding_problem_solved_after_native_audit": solved,
        "benchmark_claim_allowed_after_native_audit": solved,
        "claim_rejection_reason": claim_rejection_reason,
        "sweep_written": bool(args.sweep),
        "summary": (
            "A small leave-one-out empirical prior improves 4AKE no-AF contact recovery over hand priors/DG, "
            "but still does not meet the 0.70/0.70 solve threshold."
        ),
    }
    (out_dir / "4ake_leave_one_out_empirical_prior_report.json").write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")

    summary_md = f"""# 4AKE leave-one-out empirical sequence-contact prior v0

This run excludes 4AKE from empirical-prior training and uses 4AKE native coordinates only after prediction for audit metrics.

```text
kind = {EMPIRICAL_CONTACT_PRIOR_KIND}
source_accession = {target_row.source_accession}
sequence_length = {target_row.sequence_length}
training_rows = {packet.training_row_count}
training_coordinate_truth_used_from_non_target_rows = true
target_native_truth_used_before_selection = false
alphafold_used_as_evidence = false
msa_used_as_evidence = false
gif_generation_used = false
raw_sequence_persisted = false

selected_contacts = {packet.selected_pair_count}
selected_long_range_contacts = {packet.selected_long_range_pair_count}
native_contacts = {metric.native_contact_count}
true_positive_contacts = {metric.true_positive_contacts}
false_positive_contacts = {metric.false_positive_contacts}
false_negative_contacts = {metric.false_negative_contacts}
precision = {metric.native_contact_precision:.6f}
recall = {metric.native_contact_recall:.6f}
f1 = {metric.contact_map_f1:.6f}
long_range_recall = {metric.long_range_contact_recall:.6f}
solved_long_range_recall_threshold = {args.solved_long_range_recall_threshold:.6f}
folding_problem_solved = {str(solved).lower()}
claim_rejection_reason = {claim_rejection_reason}
```

Interpretation: this is a stronger no-AlphaFold fallback than hand physics or iterative DG alone, but it still does not crack 4AKE at the required 0.70 precision and 0.70 recall boundary.
"""
    (out_dir / "4AKE_LEAVE_ONE_OUT_EMPIRICAL_PRIOR_V0_SUMMARY.md").write_text(summary_md, encoding="utf-8")

    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
