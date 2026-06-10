#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Sequence

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_local_coevolution_expansion import LocalCoevolutionExpansionPacket, run_local_coevolution_expansion_packet
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_EXTERNAL_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "local_coevolution_expansion_v0"


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


def _flatten_row(row: dict[str, object]) -> dict[str, object]:
    flat = dict(row)
    metric = flat.pop("metric_after_native_audit", {})
    if isinstance(metric, dict):
        for key, value in metric.items():
            flat[f"metric_{key}"] = value
    return flat


def _write_packet(packet: LocalCoevolutionExpansionPacket, out_dir: Path) -> dict[str, str]:
    payload = packet.to_dict()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "local_coevolution_expansion_report.json"
    rows_path = out_dir / "local_coevolution_expansion_rows.csv"
    contacts_path = out_dir / "local_coevolution_expansion_contacts.csv"
    decisions_path = out_dir / "local_coevolution_expansion_decisions.csv"
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(rows_path, [_flatten_row(row) for row in payload["rows"]])
    _write_csv(contacts_path, payload["contacts"])
    _write_csv(decisions_path, payload["decisions"])
    return {"report": str(report_path), "rows": str(rows_path), "contacts": str(contacts_path), "decisions": str(decisions_path)}


def _summary(packet: LocalCoevolutionExpansionPacket) -> dict[str, object]:
    return {
        "kind": packet.kind,
        "source_mode": packet.source_mode,
        "row_count": packet.row_count,
        "raw_msa_available_for_true_local_mi": packet.raw_msa_available_for_true_local_mi,
        "local_mi_channel_is_proxy_not_new_msa_calculation": packet.local_mi_channel_is_proxy_not_new_msa_calculation,
        "top_safe_dca_anchor_count_requested": packet.top_safe_dca_anchor_count_requested,
        "local_window": packet.local_window,
        "threshold": packet.threshold,
        "mean_safe_anchor_count": packet.mean_safe_anchor_count,
        "mean_accepted_local_pair_count": packet.mean_accepted_local_pair_count,
        "mean_native_contact_precision_after_audit": packet.mean_native_contact_precision_after_audit,
        "mean_native_contact_recall_after_audit": packet.mean_native_contact_recall_after_audit,
        "mean_long_range_contact_recall_after_audit": packet.mean_long_range_contact_recall_after_audit,
        "mean_contact_map_f1_after_audit": packet.mean_contact_map_f1_after_audit,
        "min_row_precision_after_audit": packet.min_row_precision_after_audit,
        "min_row_recall_after_audit": packet.min_row_recall_after_audit,
        "min_row_long_range_recall_after_audit": packet.min_row_long_range_recall_after_audit,
        "local_coevolution_claim_allowed": packet.local_coevolution_claim_allowed,
        "universal_physical_law_claim_allowed": packet.universal_physical_law_claim_allowed,
        "folding_problem_solved": packet.folding_problem_solved,
        "claim_rejection_reason": packet.claim_rejection_reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_EXTERNAL_COUPLING_FILE))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--top-anchor-count", type=int, default=50)
    parser.add_argument("--window", type=int, default=5)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--degree-cap", type=int, default=5)
    parser.add_argument("--source-accession", default="", help="Optional single accession, e.g. 4AKE:A")
    args = parser.parse_args()
    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    if args.source_accession and not any(row.source_accession == args.source_accession for row in rows):
        raise SystemExit(f"No rows matched --source-accession={args.source_accession!r}")
    constraints = load_coupling_dataset(Path(args.external_coupling_file)).constraints
    packet = run_local_coevolution_expansion_packet(
        rows,
        constraints=constraints,
        evaluation_source_accessions=tuple([args.source_accession]) if args.source_accession else None,
        top_anchor_count=args.top_anchor_count,
        window=args.window,
        threshold=args.threshold,
        degree_cap=args.degree_cap,
    )
    out_dir = Path(args.out_dir)
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
        "native_truth_used_before_selection": False,
        "coordinate_truth_used_before_selection": False,
    })
    cert_path = out_dir / "local_coevolution_expansion_certificate.json"
    cert_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    summary["certificate"] = str(cert_path)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
