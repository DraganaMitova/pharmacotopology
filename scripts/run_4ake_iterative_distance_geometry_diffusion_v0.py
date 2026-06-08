#!/usr/bin/env python3
"""Run the strongest native-free iterative distance-geometry diffusion attempt on 4AKE.

No AlphaFold, no templates, no internet dependency, no GIF generation. Native
4AKE coordinates are attached only after prediction for audit metrics.
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

from pharmacotopology.folding_contact_law_features import contact_law_feature_rows_for_row  # noqa: E402
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset  # noqa: E402
from pharmacotopology.folding_iterative_distance_geometry_diffusion import (  # noqa: E402
    ITERATIVE_DISTANCE_GEOMETRY_DIFFUSION_KIND,
    run_iterative_distance_geometry_diffusion,
)
from pharmacotopology.folding_native_contact_eval import contact_map_hash, evaluate_contact_prediction  # noqa: E402
from pharmacotopology.folding_nucleus_closure_search import nucleus_closure_events_for_row  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows  # noqa: E402

DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
DEFAULT_OUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "4ake_iterative_distance_geometry_diffusion_v0"


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


def _metric_dict(packet, row, *, solved_precision_threshold: float, solved_recall_threshold: float) -> dict[str, object]:
    native_pairs = set(row.native_contact_pairs())
    selected_pairs = tuple(packet.selected_pairs)
    selected_long = {pair for pair in selected_pairs if pair[1] - pair[0] >= 24}
    native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
    long_tp = selected_long & native_long
    metric = evaluate_contact_prediction(native_pairs=native_pairs, predicted_pairs=selected_pairs)
    long_precision = round(len(long_tp) / len(selected_long), 6) if selected_long else 0.0
    long_recall = round(len(long_tp) / len(native_long), 6) if native_long else 0.0
    solved = (
        metric.native_contact_precision >= solved_precision_threshold
        and metric.native_contact_recall >= solved_recall_threshold
        and not packet.coordinate_truth_used_before_selection
        and not packet.native_truth_used_before_selection
        and not packet.alphafold_used_before_selection
        and not packet.structure_template_used_before_selection
    )
    return {
        "predicted_contact_count": metric.predicted_contact_count,
        "native_contact_count": metric.native_contact_count,
        "true_positive_contacts": metric.true_positive_contacts,
        "false_positive_contacts": metric.false_positive_contacts,
        "false_negative_contacts": metric.false_negative_contacts,
        "native_contact_precision": metric.native_contact_precision,
        "native_contact_recall": metric.native_contact_recall,
        "contact_map_f1": metric.contact_map_f1,
        "long_range_precision": long_precision,
        "long_range_recall": long_recall,
        "true_positive_long_range_contacts": len(long_tp),
        "native_long_range_contact_count": len(native_long),
        "predicted_long_range_contact_count": len(selected_long),
        "solved_precision_threshold": solved_precision_threshold,
        "solved_recall_threshold": solved_recall_threshold,
        "folding_problem_solved_after_native_audit": bool(solved),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="4AKE iterative distance-geometry diffusion no-AF run")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--coupling-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--iteration-count", type=int, default=None)
    parser.add_argument("--solved-precision-threshold", type=float, default=0.70)
    parser.add_argument("--solved-recall-threshold", type=float, default=0.70)
    parser.add_argument("--mode-sweep", action="store_true", help="opt-in: run several heavier native-free diffusion configurations")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    row = _row_by_accession(load_real_coordinate_visual_rows(Path(args.benchmark_file)), args.source_accession)
    coupling_dataset = load_coupling_dataset(Path(args.coupling_file))
    constraints = tuple(
        item for item in coupling_dataset.constraints
        if item.source_accession == row.source_accession or item.row_id == row.row_id
    )
    features = contact_law_feature_rows_for_row(row)
    events = nucleus_closure_events_for_row(row, features)

    modes = [
        {"mode_id": "balanced_dg_diffusion", "attraction_budget_fraction": 1.72, "final_budget_fraction": 2.22, "max_degree": 9},
    ]
    if args.mode_sweep:
        modes.extend([
            {"mode_id": "tight_low_degree", "attraction_budget_fraction": 1.25, "final_budget_fraction": 1.70, "max_degree": 6},
            {"mode_id": "wide_recall", "attraction_budget_fraction": 2.20, "final_budget_fraction": 3.05, "max_degree": 12},
            {"mode_id": "seed_conservative", "attraction_budget_fraction": 1.00, "final_budget_fraction": 1.38, "max_degree": 6},
            {"mode_id": "geometry_recall_degree9", "attraction_budget_fraction": 2.35, "final_budget_fraction": 3.20, "max_degree": 9},
        ])

    mode_results: list[dict[str, object]] = []
    packets = []
    for mode in modes:
        packet = run_iterative_distance_geometry_diffusion(
            row=row,
            constraints=constraints,
            events=events,
            iteration_count=args.iteration_count,
            attraction_budget_fraction=float(mode["attraction_budget_fraction"]),
            final_budget_fraction=float(mode["final_budget_fraction"]),
            max_degree=int(mode["max_degree"]),
        )
        metric = _metric_dict(
            packet,
            row,
            solved_precision_threshold=args.solved_precision_threshold,
            solved_recall_threshold=args.solved_recall_threshold,
        )
        internal_quality = round(
            0.34 * (packet.iterations[-1].mean_selected_score if packet.iterations else 0.0)
            + 0.22 * (packet.iterations[-1].mean_selected_geometry_score if packet.iterations else 0.0)
            + 0.22 * packet.mean_final_physical_prior_score
            + 0.22 * min(1.0, packet.seed_anchor_pair_count / max(1, packet.final_pair_count)),
            6,
        )
        row_result = {
            "mode_id": mode["mode_id"],
            "internal_native_free_quality_score": internal_quality,
            "attraction_budget_fraction": mode["attraction_budget_fraction"],
            "final_budget_fraction": mode["final_budget_fraction"],
            "max_degree": mode["max_degree"],
            **packet.to_report_dict(),
            **metric,
        }
        mode_results.append(row_result)
        packets.append((mode["mode_id"], packet, metric, internal_quality))

    # The selected run is chosen only by native-free internal quality. The audit
    # best is reported separately but never used to select the claim run.
    selected_mode_id, selected_packet, selected_metric, selected_quality = max(
        packets, key=lambda item: (item[3], -item[1].final_pair_count, str(item[0]))
    )
    audit_best = max(mode_results, key=lambda item: (float(item["contact_map_f1"]), float(item["native_contact_recall"])))
    selected_solved = bool(selected_metric["folding_problem_solved_after_native_audit"])
    claim_rejection = "none" if selected_solved else "native_free_iterative_diffusion_did_not_reach_0_70_precision_and_recall"
    if selected_packet.seed_anchor_pair_count <= 0:
        claim_rejection = "no_safe_external_dca_anchors_available"

    summary = {
        "kind": ITERATIVE_DISTANCE_GEOMETRY_DIFFUSION_KIND,
        "source_accession": row.source_accession,
        "row_id": row.row_id,
        "sequence_hash": row.sequence_sha256,
        "sequence_length": row.sequence_length,
        "raw_sequence_persisted": False,
        "gif_generation_used": False,
        "alphafold_used_as_evidence": False,
        "structure_template_used_as_evidence": False,
        "internet_required": False,
        "native_truth_used_before_selection": False,
        "coordinate_truth_used_before_selection": False,
        "selected_mode_id_by_native_free_quality": selected_mode_id,
        "selected_internal_native_free_quality_score": selected_quality,
        "selected_packet": selected_packet.to_report_dict(),
        "selected_metric_after_native_audit": selected_metric,
        "benchmark_claim_allowed_after_native_audit": selected_solved,
        "folding_problem_solved": selected_solved,
        "claim_rejection_reason": claim_rejection,
        "audit_best_mode_not_used_for_selection": audit_best,
        "mode_results": mode_results,
    }

    (out_dir / "4ake_iterative_distance_geometry_diffusion_report.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    _write_csv(out_dir / "4ake_iterative_distance_geometry_diffusion_modes.csv", mode_results)
    _write_csv(
        out_dir / "4ake_iterative_distance_geometry_diffusion_iterations.csv",
        [
            {"mode_id": selected_mode_id, **iteration.to_dict()}
            for iteration in selected_packet.iterations
        ],
    )
    _write_csv(
        out_dir / "4ake_iterative_distance_geometry_diffusion_selected_contacts.csv",
        [score.to_dict() for score in selected_packet.scores if score.selected],
    )
    _write_csv(
        out_dir / "4ake_iterative_distance_geometry_diffusion_top_scores.csv",
        [score.to_dict() for score in selected_packet.scores[:700]],
    )
    md = [
        "# 4AKE iterative distance-geometry diffusion v0",
        "",
        f"selected_mode_id_by_native_free_quality: `{selected_mode_id}`",
        f"folding_problem_solved: `{selected_solved}`",
        f"claim_rejection_reason: `{claim_rejection}`",
        "gif_generation_used: `False`",
        "alphafold_used_as_evidence: `False`",
        "structure_template_used_as_evidence: `False`",
        "native_truth_used_before_selection: `False`",
        "",
        "## Selected metric after native audit",
        f"precision: `{selected_metric['native_contact_precision']}`",
        f"recall: `{selected_metric['native_contact_recall']}`",
        f"long_range_recall: `{selected_metric['long_range_recall']}`",
        f"predicted_contact_count: `{selected_metric['predicted_contact_count']}`",
        f"true_positive_contacts: `{selected_metric['true_positive_contacts']}`",
        "",
        "## Audit best mode, not used for selection",
        f"mode: `{audit_best['mode_id']}`",
        f"precision: `{audit_best['native_contact_precision']}`",
        f"recall: `{audit_best['native_contact_recall']}`",
        f"f1: `{audit_best['contact_map_f1']}`",
    ]
    (out_dir / "4AKE_ITERATIVE_DG_DIFFUSION_V0_SUMMARY.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps({
        "report_path": str(out_dir / "4ake_iterative_distance_geometry_diffusion_report.json"),
        "selected_mode_id": selected_mode_id,
        "folding_problem_solved": selected_solved,
        "selected_precision": selected_metric["native_contact_precision"],
        "selected_recall": selected_metric["native_contact_recall"],
        "audit_best_mode": audit_best["mode_id"],
        "audit_best_f1": audit_best["contact_map_f1"],
        "selected_contact_map_hash": contact_map_hash(selected_packet.selected_pairs),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
