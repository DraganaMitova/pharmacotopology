#!/usr/bin/env python3
from __future__ import annotations

"""Run the native-free global geometry-field challenge.

Bounded by design: no internet, no predictor subprocesses, no GIF generation,
no full pytest suite. Native coordinates are used only after selection for audit.
"""

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

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset  # noqa: E402
from pharmacotopology.folding_native_free_geometry_field import (  # noqa: E402
    EXTERNAL_DCA_GEOMETRY_FIELD_MODE,
    NATIVE_FREE_GEOMETRY_FIELD_KIND,
    ORACLE_ANCHOR_GEOMETRY_FIELD_CONTROL_MODE,
    PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE,
    NativeFreeGeometryFieldPacket,
    run_native_free_geometry_field_packet,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows  # noqa: E402

DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_EXTERNAL_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
DEFAULT_ORACLE_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_couplings.locked.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "native_free_geometry_field_v0"


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


def _flatten_metric_row(row: dict[str, object]) -> dict[str, object]:
    flattened = dict(row)
    metric = flattened.pop("metric_after_native_audit", {})
    best_metric = flattened.pop("best_solution_metric_after_native_audit", {})
    if isinstance(metric, dict):
        for key, value in metric.items():
            flattened[f"metric_{key}"] = value
    if isinstance(best_metric, dict):
        for key, value in best_metric.items():
            flattened[f"best_solution_metric_{key}"] = value
    return flattened


def _flatten_solution(solution: dict[str, object]) -> dict[str, object]:
    flattened = dict(solution)
    metric = flattened.pop("metric_after_native_audit", {})
    if isinstance(metric, dict):
        for key, value in metric.items():
            flattened[f"metric_{key}"] = value
    factor_ids = flattened.get("selected_factor_ids")
    if isinstance(factor_ids, list):
        flattened["selected_factor_ids"] = ";".join(str(item) for item in factor_ids[:40])
    return flattened


def _compact_factors(factors: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    compact: list[dict[str, object]] = []
    for factor in factors:
        row = dict(factor)
        contacts = row.pop("contacts", [])
        if isinstance(contacts, list):
            row["contacts_preview"] = ";".join(
                f"{pair[0]}-{pair[1]}" for pair in contacts[:12] if isinstance(pair, list) and len(pair) == 2
            )
            row["contacts_preview_count"] = min(len(contacts), 12)
        compact.append(row)
    return compact


def _write_packet(prefix: str, packet: NativeFreeGeometryFieldPacket, out_dir: Path) -> dict[str, str]:
    payload = packet.to_dict()
    report_path = out_dir / f"{prefix}_report.json"
    rows_path = out_dir / f"{prefix}_rows.csv"
    factors_path = out_dir / f"{prefix}_factors.csv"
    solutions_path = out_dir / f"{prefix}_solutions.csv"
    decisions_path = out_dir / f"{prefix}_decisions.csv"
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(rows_path, [_flatten_metric_row(row) for row in payload["rows"]])
    _write_csv(factors_path, _compact_factors(payload["factors"]))
    _write_csv(solutions_path, [_flatten_solution(row) for row in payload["solutions"]])
    _write_csv(decisions_path, payload["decisions"])
    return {
        "report": str(report_path),
        "rows": str(rows_path),
        "factors": str(factors_path),
        "solutions": str(solutions_path),
        "decisions": str(decisions_path),
    }


def _summary_for_packet(packet: NativeFreeGeometryFieldPacket) -> dict[str, object]:
    return {
        "row_count": packet.row_count,
        "mean_native_contact_precision_after_audit": packet.mean_native_contact_precision_after_audit,
        "mean_native_contact_recall_after_audit": packet.mean_native_contact_recall_after_audit,
        "mean_long_range_contact_recall_after_audit": packet.mean_long_range_contact_recall_after_audit,
        "mean_contact_map_f1_after_audit": packet.mean_contact_map_f1_after_audit,
        "min_row_precision_after_audit": packet.min_row_precision_after_audit,
        "min_row_recall_after_audit": packet.min_row_recall_after_audit,
        "min_row_long_range_recall_after_audit": packet.min_row_long_range_recall_after_audit,
        "native_free_geometry_field_claim_allowed": packet.native_free_geometry_field_claim_allowed,
        "universal_physical_law_claim_allowed": packet.universal_physical_law_claim_allowed,
        "folding_problem_solved": packet.folding_problem_solved,
        "claim_rejection_reason": packet.claim_rejection_reason,
        "coordinate_truth_used_before_selection": packet.coordinate_truth_used_before_selection,
        "native_truth_used_before_selection": packet.native_truth_used_before_selection,
        "structure_model_used_before_selection": packet.structure_model_used_before_selection,
        "learned_geometry_prior_used_before_selection": packet.learned_geometry_prior_used_before_selection,
        "msa_dca_used_before_selection": packet.msa_dca_used_before_selection,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_EXTERNAL_COUPLING_FILE))
    parser.add_argument("--oracle-coupling-file", default=str(DEFAULT_ORACLE_COUPLING_FILE))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--mode", choices=["pure", "external", "both"], default="both")
    parser.add_argument("--include-oracle-control", action="store_true")
    parser.add_argument("--pair-pool-multiplier", type=float, default=6.0)
    parser.add_argument("--max-selected-contact-multiplier", type=float, default=2.5)
    parser.add_argument("--ensemble-size", type=int, default=6)
    args = parser.parse_args()

    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    external_constraints = load_coupling_dataset(Path(args.external_coupling_file)).constraints
    oracle_constraints = load_coupling_dataset(Path(args.oracle_coupling_file)).constraints
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    requested_modes: list[tuple[str, str, Sequence[object]]] = []
    if args.mode in {"pure", "both"}:
        requested_modes.append(("pure_sequence_symbolic_geometry", PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE, ()))
    if args.mode in {"external", "both"}:
        requested_modes.append(("external_dca_geometry_field", EXTERNAL_DCA_GEOMETRY_FIELD_MODE, external_constraints))
    if args.include_oracle_control:
        requested_modes.append(("oracle_anchor_geometry_field_control", ORACLE_ANCHOR_GEOMETRY_FIELD_CONTROL_MODE, oracle_constraints))

    summary: dict[str, object] = {
        "kind": NATIVE_FREE_GEOMETRY_FIELD_KIND,
        "out_dir": str(out_dir),
        "gif_generation_used": False,
        "internet_required": False,
        "predictor_subprocesses_used": False,
        "full_suite_run": False,
        "artifacts": {},
        "modes": {},
    }

    for prefix, source_mode, constraints in requested_modes:
        packet = run_native_free_geometry_field_packet(
            rows,
            constraints=constraints,  # type: ignore[arg-type]
            source_mode=source_mode,
            pair_pool_multiplier=args.pair_pool_multiplier,
            max_selected_contact_multiplier=args.max_selected_contact_multiplier,
            ensemble_size=args.ensemble_size,
        )
        paths = _write_packet(prefix, packet, out_dir)
        summary["modes"][prefix] = _summary_for_packet(packet)  # type: ignore[index]
        summary["artifacts"][prefix] = paths  # type: ignore[index]

    certificate_path = out_dir / "native_free_geometry_field_certificate.json"
    certificate_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    summary["certificate"] = str(certificate_path)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
