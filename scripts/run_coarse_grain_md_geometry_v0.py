#!/usr/bin/env python3
from __future__ import annotations

"""Run the native-free coarse-grain MD global-geometry challenge.

Default mode is bounded multistart DCA-restrained global geometry. This is not a
full suite and does not call OpenMM/GROMACS or predictor subprocesses. It tests
the structure-first pivot: restraints -> global C-alpha geometry -> extracted
contacts -> native audit only after selection.
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

from pharmacotopology.folding_coarse_grain_md_geometry import (  # noqa: E402
    COARSE_GRAIN_MD_GEOMETRY_KIND,
    EXTERNAL_DCA_MD_MODE,
    MULTISTART_DCA_MD_MODE,
    PURE_SEQUENCE_MD_MODE,
    CoarseGrainMDGeometryPacket,
    run_coarse_grain_md_geometry_packet,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows  # noqa: E402

DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_EXTERNAL_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "coarse_grain_md_geometry_v0"


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
    flattened = dict(row)
    _flatten_metric("metric", flattened, "metric_after_native_audit")
    return flattened


def _write_packet(prefix: str, packet: CoarseGrainMDGeometryPacket, out_dir: Path) -> dict[str, str]:
    payload = packet.to_dict()
    report_path = out_dir / f"{prefix}_report.json"
    rows_path = out_dir / f"{prefix}_rows.csv"
    decisions_path = out_dir / f"{prefix}_decisions.csv"
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(rows_path, [_flatten_row(row) for row in payload["rows"]])
    _write_csv(decisions_path, payload["decisions"])
    return {"report": str(report_path), "rows": str(rows_path), "decisions": str(decisions_path)}


def _summary_for_packet(packet: CoarseGrainMDGeometryPacket) -> dict[str, object]:
    return {
        "row_count": packet.row_count,
        "source_mode": packet.source_mode,
        "direct_global_structure_generation_included": packet.direct_global_structure_generation_included,
        "contacts_predicted_before_structure": packet.contacts_predicted_before_structure,
        "mean_native_contact_precision_after_audit": packet.mean_native_contact_precision_after_audit,
        "mean_native_contact_recall_after_audit": packet.mean_native_contact_recall_after_audit,
        "mean_long_range_contact_recall_after_audit": packet.mean_long_range_contact_recall_after_audit,
        "mean_contact_map_f1_after_audit": packet.mean_contact_map_f1_after_audit,
        "min_row_precision_after_audit": packet.min_row_precision_after_audit,
        "min_row_recall_after_audit": packet.min_row_recall_after_audit,
        "min_row_long_range_recall_after_audit": packet.min_row_long_range_recall_after_audit,
        "mean_f1_margin_vs_best_control": packet.mean_f1_margin_vs_best_control,
        "mean_long_range_recall_margin_vs_best_control": packet.mean_long_range_recall_margin_vs_best_control,
        "md_geometry_claim_allowed": packet.md_geometry_claim_allowed,
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
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--mode", choices=["pure", "dca", "multistart", "all"], default="multistart")
    parser.add_argument("--md-steps", type=int, default=72)
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--max-restraints", type=int, default=96)
    args = parser.parse_args()

    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    constraints = load_coupling_dataset(Path(args.external_coupling_file)).constraints
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    requested: list[tuple[str, str, Sequence[object]]] = []
    if args.mode in {"pure", "all"}:
        requested.append(("pure_sequence_md_geometry", PURE_SEQUENCE_MD_MODE, ()))
    if args.mode in {"dca", "all"}:
        requested.append(("external_dca_md_geometry", EXTERNAL_DCA_MD_MODE, constraints))
    if args.mode in {"multistart", "all"}:
        requested.append(("external_dca_multistart_md_geometry", MULTISTART_DCA_MD_MODE, constraints))

    summary: dict[str, object] = {
        "kind": COARSE_GRAIN_MD_GEOMETRY_KIND,
        "out_dir": str(out_dir),
        "gif_generation_used": False,
        "internet_required": False,
        "predictor_subprocesses_used": False,
        "full_suite_run": False,
        "atomistic_md_engine_used": False,
        "dependency_free_md_style_relaxation_used": True,
        "artifacts": {},
        "modes": {},
    }
    for prefix, source_mode, mode_constraints in requested:
        packet = run_coarse_grain_md_geometry_packet(
            rows,
            constraints=mode_constraints,  # type: ignore[arg-type]
            source_mode=source_mode,
            md_steps=args.md_steps,
            restarts=args.restarts,
            max_restraints=args.max_restraints,
        )
        summary["modes"][prefix] = _summary_for_packet(packet)  # type: ignore[index]
        summary["artifacts"][prefix] = _write_packet(prefix, packet, out_dir)  # type: ignore[index]

    certificate_path = out_dir / "coarse_grain_md_geometry_certificate.json"
    certificate_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    summary["certificate"] = str(certificate_path)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
