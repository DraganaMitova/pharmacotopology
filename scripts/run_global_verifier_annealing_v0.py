#!/usr/bin/env python3
from __future__ import annotations

"""Run the native-free global verifier / structure-level annealing challenge."""

import argparse
import csv
import json
from pathlib import Path
from typing import Sequence

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_global_verifier_annealing import (
    GlobalVerifierAnnealingPacket,
    run_global_verifier_annealing_packet,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_EXTERNAL_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "global_verifier_annealing_v0"


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


def _flatten_metric(prefix: str, payload: dict[str, object], metric_key: str) -> None:
    metric = payload.pop(metric_key, {})
    if isinstance(metric, dict):
        for key, value in metric.items():
            payload[f"{prefix}_{key}"] = value


def _flatten_row(row: dict[str, object]) -> dict[str, object]:
    flat = dict(row)
    _flatten_metric("metric", flat, "metric_after_native_audit")
    return flat


def _write_packet(packet: GlobalVerifierAnnealingPacket, out_dir: Path) -> dict[str, str]:
    payload = packet.to_dict()
    report_path = out_dir / "global_verifier_annealing_report.json"
    rows_path = out_dir / "global_verifier_annealing_rows.csv"
    decisions_path = out_dir / "global_verifier_annealing_decisions.csv"
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(rows_path, [_flatten_row(row) for row in payload["rows"]])
    _write_csv(decisions_path, payload["decisions"])
    return {"report": str(report_path), "rows": str(rows_path), "decisions": str(decisions_path)}


def _summary(packet: GlobalVerifierAnnealingPacket) -> dict[str, object]:
    return {
        "kind": packet.kind,
        "source_mode": packet.source_mode,
        "row_count": packet.row_count,
        "global_verifier_included": packet.global_verifier_included,
        "statistical_potential_surrogate_included": packet.statistical_potential_surrogate_included,
        "degree_distribution_coherence_included": packet.degree_distribution_coherence_included,
        "loop_patch_coherence_included": packet.loop_patch_coherence_included,
        "monte_carlo_global_candidate_selection_included": packet.monte_carlo_global_candidate_selection_included,
        "contacts_predicted_before_structure": packet.contacts_predicted_before_structure,
        "mean_selected_global_score": packet.mean_selected_global_score,
        "mean_statistical_potential_score": packet.mean_statistical_potential_score,
        "mean_degree_coherence_score": packet.mean_degree_coherence_score,
        "mean_restraint_satisfaction_score": packet.mean_restraint_satisfaction_score,
        "mean_native_contact_precision_after_audit": packet.mean_native_contact_precision_after_audit,
        "mean_native_contact_recall_after_audit": packet.mean_native_contact_recall_after_audit,
        "mean_long_range_contact_recall_after_audit": packet.mean_long_range_contact_recall_after_audit,
        "mean_contact_map_f1_after_audit": packet.mean_contact_map_f1_after_audit,
        "min_row_precision_after_audit": packet.min_row_precision_after_audit,
        "min_row_recall_after_audit": packet.min_row_recall_after_audit,
        "min_row_long_range_recall_after_audit": packet.min_row_long_range_recall_after_audit,
        "mean_f1_margin_vs_best_control": packet.mean_f1_margin_vs_best_control,
        "mean_long_range_recall_margin_vs_best_control": packet.mean_long_range_recall_margin_vs_best_control,
        "global_verifier_claim_allowed": packet.global_verifier_claim_allowed,
        "universal_physical_law_claim_allowed": packet.universal_physical_law_claim_allowed,
        "folding_problem_solved": packet.folding_problem_solved,
        "claim_rejection_reason": packet.claim_rejection_reason,
        "coordinate_truth_used_before_selection": packet.coordinate_truth_used_before_selection,
        "native_truth_used_before_selection": packet.native_truth_used_before_selection,
        "structure_model_used_before_selection": packet.structure_model_used_before_selection,
        "learned_geometry_prior_used_before_selection": packet.learned_geometry_prior_used_before_selection,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_EXTERNAL_COUPLING_FILE))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--md-steps", type=int, default=84)
    parser.add_argument("--candidate-count", type=int, default=5)
    parser.add_argument("--max-direct", type=int, default=0)
    parser.add_argument("--max-sequence-closure", type=int, default=0)
    parser.add_argument("--max-geodesic", type=int, default=0)
    args = parser.parse_args()
    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    constraints = load_coupling_dataset(Path(args.external_coupling_file)).constraints
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    packet = run_global_verifier_annealing_packet(
        rows,
        constraints=constraints,
        md_steps=args.md_steps,
        candidate_count=args.candidate_count,
        max_direct=args.max_direct or None,
        max_sequence_closure=args.max_sequence_closure or None,
        max_geodesic=args.max_geodesic or None,
    )
    artifacts = _write_packet(packet, out_dir)
    summary = _summary(packet)
    summary.update({
        "out_dir": str(out_dir),
        "artifacts": artifacts,
        "full_suite_run": False,
        "gif_generation_used": False,
        "internet_required": False,
        "predictor_subprocesses_used": False,
        "atomistic_md_engine_used": False,
        "dependency_free_md_style_relaxation_used": True,
    })
    certificate_path = out_dir / "global_verifier_annealing_certificate.json"
    certificate_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    summary["certificate"] = str(certificate_path)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
