#!/usr/bin/env python3
from __future__ import annotations

"""Run the global contact factor-graph + ensemble challenge.

Default behavior is intentionally bounded:

* sequence-only factor graph over the locked 8-row coordinate benchmark;
* optional 4AKE learned-geometry factor graph if the included ESMFold PDB exists;
* no internet calls, no predictor subprocesses, no GIF generation, no full suite.
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

from pharmacotopology.folding_global_factor_graph_ensemble import (  # noqa: E402
    GLOBAL_FACTOR_GRAPH_KIND,
    LEARNED_GEOMETRY_MODE,
    SEQUENCE_ONLY_MODE,
    GlobalFactorGraphModelSpec,
    run_global_factor_graph_packet,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)

DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "global_factor_graph_ensemble_v0"
DEFAULT_ESMFOLD_4AKE_PDB = REPO_ROOT / "external_msa_free_predictors" / "tryhard_runs_live" / "4ake_esmfold_api.pdb"


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


def _write_packet(prefix: str, packet, out_dir: Path) -> dict[str, str]:
    payload = packet.to_dict()
    report_path = out_dir / f"{prefix}_report.json"
    rows_path = out_dir / f"{prefix}_rows.csv"
    factors_path = out_dir / f"{prefix}_factors.csv"
    mutex_path = out_dir / f"{prefix}_mutex_groups.csv"
    solutions_path = out_dir / f"{prefix}_solutions.csv"
    decisions_path = out_dir / f"{prefix}_decisions.csv"
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(rows_path, [_flatten_metric_row(row) for row in payload["rows"]])
    _write_csv(factors_path, _compact_factors(payload["factors"]))
    _write_csv(mutex_path, payload["mutex_groups"])
    _write_csv(solutions_path, [_flatten_solution(row) for row in payload["solutions"]])
    _write_csv(decisions_path, payload["decisions"])
    return {
        "report": str(report_path),
        "rows": str(rows_path),
        "factors": str(factors_path),
        "mutex_groups": str(mutex_path),
        "solutions": str(solutions_path),
        "decisions": str(decisions_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--max-sequence-separation", type=int, default=160)
    parser.add_argument("--ensemble-size", type=int, default=6)
    parser.add_argument("--skip-sequence-only", action="store_true")
    parser.add_argument("--skip-4ake-esmfold", action="store_true")
    parser.add_argument("--esmfold-4ake-pdb", default=str(DEFAULT_ESMFOLD_4AKE_PDB))
    args = parser.parse_args()

    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary: dict[str, object] = {
        "kind": GLOBAL_FACTOR_GRAPH_KIND,
        "out_dir": str(out_dir),
        "sequence_only": None,
        "learned_geometry_4ake": None,
        "artifacts": {},
    }

    if not args.skip_sequence_only:
        sequence_packet = run_global_factor_graph_packet(
            rows,
            source_mode=SEQUENCE_ONLY_MODE,
            max_sequence_separation=args.max_sequence_separation,
            ensemble_size=args.ensemble_size,
        )
        paths = _write_packet("sequence_only_global_factor_graph", sequence_packet, out_dir)
        summary["sequence_only"] = {
            "row_count": sequence_packet.row_count,
            "mean_native_contact_precision_after_audit": sequence_packet.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": sequence_packet.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": sequence_packet.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": sequence_packet.mean_contact_map_f1_after_audit,
            "factor_graph_ensemble_claim_allowed": sequence_packet.factor_graph_ensemble_claim_allowed,
            "universal_physical_law_claim_allowed": sequence_packet.universal_physical_law_claim_allowed,
            "folding_problem_solved": sequence_packet.folding_problem_solved,
            "claim_rejection_reason": sequence_packet.claim_rejection_reason,
        }
        summary["artifacts"]["sequence_only"] = paths

    esmfold_path = Path(args.esmfold_4ake_pdb)
    if not args.skip_4ake_esmfold and esmfold_path.exists():
        row_4ake = next(row for row in rows if row.source_accession == "4AKE:A")
        learned_packet = run_global_factor_graph_packet(
            [row_4ake],
            source_mode=LEARNED_GEOMETRY_MODE,
            model_specs_by_accession={
                "4AKE:A": [
                    GlobalFactorGraphModelSpec(
                        source_id="esmfold_single_sequence_4ake",
                        pdb_path=str(esmfold_path),
                        chain_id="A",
                    )
                ]
            },
            max_sequence_separation=None,
            ensemble_size=args.ensemble_size,
        )
        paths = _write_packet("learned_geometry_4ake_global_factor_graph", learned_packet, out_dir)
        summary["learned_geometry_4ake"] = {
            "row_count": learned_packet.row_count,
            "mean_native_contact_precision_after_audit": learned_packet.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": learned_packet.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": learned_packet.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": learned_packet.mean_contact_map_f1_after_audit,
            "factor_graph_ensemble_claim_allowed": learned_packet.factor_graph_ensemble_claim_allowed,
            "universal_physical_law_claim_allowed": learned_packet.universal_physical_law_claim_allowed,
            "folding_problem_solved": learned_packet.folding_problem_solved,
            "claim_rejection_reason": learned_packet.claim_rejection_reason,
        }
        summary["artifacts"]["learned_geometry_4ake"] = paths

    certificate_path = out_dir / "global_factor_graph_ensemble_certificate.json"
    certificate_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    summary["certificate"] = str(certificate_path)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
