#!/usr/bin/env python3
from __future__ import annotations

"""Run the consensus-filtered geodesic restrained MD-style geometry challenge."""

import argparse
import csv
import json
from pathlib import Path
from typing import Sequence

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_geodesic_consensus_filter_md import (
    GeodesicConsensusFilterPacket,
    run_geodesic_consensus_filter_packet,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_EXTERNAL_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "geodesic_consensus_filter_md_v0"


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


def _flatten_md_row(row: dict[str, object]) -> dict[str, object]:
    flat = dict(row)
    _flatten_metric("metric", flat, "metric_after_native_audit")
    return flat


def _write_packet(packet: GeodesicConsensusFilterPacket, out_dir: Path) -> dict[str, str]:
    payload = packet.to_dict()
    report_path = out_dir / "geodesic_consensus_filter_report.json"
    rows_path = out_dir / "geodesic_consensus_filter_rows.csv"
    restraints_path = out_dir / "geodesic_consensus_filter_restraints.csv"
    md_rows_path = out_dir / "geodesic_consensus_filter_md_rows.csv"
    md_decisions_path = out_dir / "geodesic_consensus_filter_md_decisions.csv"
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(rows_path, payload["filter_rows"])
    _write_csv(restraints_path, payload["restraints"])
    md_payload = payload["md_packet"]
    _write_csv(md_rows_path, [_flatten_md_row(row) for row in md_payload["rows"]])
    _write_csv(md_decisions_path, md_payload["decisions"])
    return {
        "report": str(report_path),
        "filter_rows": str(rows_path),
        "restraints": str(restraints_path),
        "md_rows": str(md_rows_path),
        "md_decisions": str(md_decisions_path),
    }


def _summary(packet: GeodesicConsensusFilterPacket) -> dict[str, object]:
    md = packet.md_packet
    return {
        "kind": packet.kind,
        "source_mode": packet.source_mode,
        "row_count": packet.row_count,
        "mean_total_restraint_count": packet.mean_total_restraint_count,
        "mean_accepted_geodesic_contact_count": packet.mean_accepted_geodesic_contact_count,
        "mean_rejected_geodesic_contact_count": packet.mean_rejected_geodesic_contact_count,
        "geodesic_distance_filter_included": packet.geodesic_distance_filter_included,
        "bending_energy_filter_included": packet.bending_energy_filter_included,
        "consensus_geodesic_filter_included": packet.consensus_geodesic_filter_included,
        "direct_global_structure_generation_included": packet.direct_global_structure_generation_included,
        "contacts_predicted_before_structure": packet.contacts_predicted_before_structure,
        "mean_native_contact_precision_after_audit": md.mean_native_contact_precision_after_audit,
        "mean_native_contact_recall_after_audit": md.mean_native_contact_recall_after_audit,
        "mean_long_range_contact_recall_after_audit": md.mean_long_range_contact_recall_after_audit,
        "mean_contact_map_f1_after_audit": md.mean_contact_map_f1_after_audit,
        "min_row_precision_after_audit": md.min_row_precision_after_audit,
        "min_row_recall_after_audit": md.min_row_recall_after_audit,
        "min_row_long_range_recall_after_audit": md.min_row_long_range_recall_after_audit,
        "geodesic_consensus_filter_claim_allowed": packet.geodesic_consensus_filter_claim_allowed,
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
    parser.add_argument("--md-steps", type=int, default=72)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--max-direct", type=int, default=0)
    parser.add_argument("--max-sequence-closure", type=int, default=0)
    parser.add_argument("--max-geodesic", type=int, default=0)
    parser.add_argument("--max-geodesic-distance", type=int, default=30)
    parser.add_argument("--max-bending-energy", type=float, default=0.50)
    parser.add_argument("--consensus-vote-threshold", type=int, default=3)
    args = parser.parse_args()

    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    constraints = load_coupling_dataset(Path(args.external_coupling_file)).constraints
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    packet = run_geodesic_consensus_filter_packet(
        rows,
        constraints=constraints,
        md_steps=args.md_steps,
        restarts=args.restarts,
        max_direct=args.max_direct or None,
        max_sequence_closure=args.max_sequence_closure or None,
        max_geodesic=args.max_geodesic or None,
        max_geodesic_distance=args.max_geodesic_distance,
        max_bending_energy=args.max_bending_energy,
        consensus_vote_threshold=args.consensus_vote_threshold,
    )
    artifacts = _write_packet(packet, out_dir)
    summary = _summary(packet)
    summary.update(
        {
            "out_dir": str(out_dir),
            "artifacts": artifacts,
            "full_suite_run": False,
            "gif_generation_used": False,
            "internet_required": False,
            "predictor_subprocesses_used": False,
            "atomistic_md_engine_used": False,
            "dependency_free_md_style_relaxation_used": True,
        }
    )
    certificate_path = out_dir / "geodesic_consensus_filter_certificate.json"
    certificate_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    summary["certificate"] = str(certificate_path)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
