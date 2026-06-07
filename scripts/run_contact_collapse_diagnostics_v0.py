from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_contact_law_features import (
    contact_law_feature_rows,
    feature_rows_by_row_id,
)
from pharmacotopology.folding_coupling_nucleus_selector import (
    build_coupling_nucleus_context,
    select_coupling_trace_loop_self_deciding_frontier_expanded_events,
    self_deciding_frontier_expansion_rows,
)
from pharmacotopology.folding_event_region_contact_collapse import (
    DEFAULT_BALANCED_PAIRS_PER_EVENT,
    DEFAULT_PRECISION_MAX_PAIRS_PER_EVENT,
    DEFAULT_RECALL_MAX_PAIRS_PER_EVENT,
    EVENT_REGION_CONTACT_COLLAPSE_BOUNDARY,
    EVENT_REGION_CONTACT_COLLAPSE_KIND,
    SELF_DECIDING_STRATEGY_NAME,
    collapse_row_event_regions,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_native_contact_eval import contact_map_hash
from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
FRONTIER_FILE = OUTPUT_DIR / "external_coupling_trace_loop_frontier.csv"


def _csv_write(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _frontier_event_ids_by_row() -> dict[str, list[str]]:
    if not FRONTIER_FILE.exists():
        raise FileNotFoundError(
            f"missing frontier file: {FRONTIER_FILE}. Run the external trace-loop benchmark first."
        )
    output: dict[str, list[str]] = {}
    with FRONTIER_FILE.open(newline="", encoding="utf-8") as handle:
        for item in csv.DictReader(handle):
            row_id = str(item.get("row_id", ""))
            event_id = str(item.get("event_id", ""))
            if row_id and event_id:
                output.setdefault(row_id, []).append(event_id)
    return output


def _selected_frontier_events(
    *,
    row: RealCoordinateVisualRow,
    event_ids_by_row: Mapping[str, Sequence[str]],
    event_by_id: Mapping[str, NucleusClosureEvent],
) -> tuple[NucleusClosureEvent, ...]:
    event_ids = tuple(event_ids_by_row.get(row.row_id, ()))
    return tuple(event_by_id[event_id] for event_id in event_ids if event_id in event_by_id)


def _row_result_to_report(
    result,
) -> dict[str, object]:
    collapsed_pairs = tuple(pair.pair() for pair in result.collapsed_pairs)
    return {
        **result.evaluation.to_dict(),
        "collapsed_contact_map_hash": contact_map_hash(collapsed_pairs),
        "uncollapsed_region_contact_map_hash": contact_map_hash(result.uncollapsed_pairs),
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = tuple(load_real_coordinate_visual_rows(BENCHMARK_FILE))
    row_by_id = {row.row_id: row for row in rows}
    dataset = load_coupling_dataset(COUPLING_FILE)
    constraints_by_row = dataset.constraints_by_row_id()
    context = build_coupling_nucleus_context(rows=rows, coupling_dataset=dataset)
    features_by_row = feature_rows_by_row_id(contact_law_feature_rows(rows))
    event_ids_by_row = _frontier_event_ids_by_row()

    pair_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    strategy_reports: dict[str, list[dict[str, object]]] = {
        SELF_DECIDING_STRATEGY_NAME: [],
        "frontier_precision": [],
        "frontier_balanced": [],
        "frontier_internal_gap_balanced": [],
        "frontier_recall": [],
        "ridge_coupling": [],
    }
    one_cll_reports: dict[str, dict[str, object]] = {}

    for strategy in strategy_reports:
        for row_id, event_ids in sorted(event_ids_by_row.items()):
            row = row_by_id.get(row_id)
            if row is None:
                continue
            events = _selected_frontier_events(
                row=row,
                event_ids_by_row={row_id: event_ids},
                event_by_id=context.event_by_id,
            )
            if not events:
                continue
            if strategy == SELF_DECIDING_STRATEGY_NAME:
                max_pairs = 0
                min_pairs = 0
            elif strategy == "frontier_recall":
                max_pairs = DEFAULT_RECALL_MAX_PAIRS_PER_EVENT
                min_pairs = 8
            elif strategy == "frontier_balanced":
                max_pairs = DEFAULT_BALANCED_PAIRS_PER_EVENT
                min_pairs = DEFAULT_BALANCED_PAIRS_PER_EVENT
            elif strategy == "frontier_internal_gap_balanced":
                max_pairs = DEFAULT_BALANCED_PAIRS_PER_EVENT
                min_pairs = 1
            else:
                max_pairs = DEFAULT_PRECISION_MAX_PAIRS_PER_EVENT
                min_pairs = 1
            result = collapse_row_event_regions(
                row=row,
                events=events,
                row_features=features_by_row.get(row.row_id, ()),
                row_constraints=constraints_by_row.get(row.row_id, ()),
                collapse_strategy=strategy,
                min_pairs_per_event=min_pairs,
                max_pairs_per_event=max_pairs,
            )
            row_report = _row_result_to_report(result)
            strategy_reports[strategy].append(row_report)
            if row.source_accession == "1CLL:A":
                one_cll_reports[strategy] = row_report
            for pair in result.collapsed_pairs:
                pair_rows.append(pair.to_dict())
            for summary in result.event_summaries:
                event_rows.append(summary.to_dict())

    one_cll_budget_probe: list[dict[str, object]] = []
    one_cll_row = next((row for row in rows if row.source_accession == "1CLL:A"), None)
    if one_cll_row is not None:
        one_cll_events = _selected_frontier_events(
            row=one_cll_row,
            event_ids_by_row=event_ids_by_row,
            event_by_id=context.event_by_id,
        )
        for budget in range(DEFAULT_BALANCED_PAIRS_PER_EVENT, DEFAULT_RECALL_MAX_PAIRS_PER_EVENT + 1):
            budget_result = collapse_row_event_regions(
                row=one_cll_row,
                events=one_cll_events,
                row_features=features_by_row.get(one_cll_row.row_id, ()),
                row_constraints=constraints_by_row.get(one_cll_row.row_id, ()),
                collapse_strategy="frontier_balanced",
                min_pairs_per_event=budget,
                max_pairs_per_event=budget,
            )
            row_report = _row_result_to_report(budget_result)
            row_report["pairs_per_event_budget"] = budget
            one_cll_budget_probe.append(row_report)

    seed_frontier_events_list = []
    for row_id, event_ids in sorted(event_ids_by_row.items()):
        row = row_by_id.get(row_id)
        if row is None:
            continue
        seed_frontier_events_list.extend(
            _selected_frontier_events(
                row=row,
                event_ids_by_row={row_id: event_ids},
                event_by_id=context.event_by_id,
            )
        )
    seed_frontier_events = tuple(seed_frontier_events_list)
    self_expanded_frontier_events = (
        select_coupling_trace_loop_self_deciding_frontier_expanded_events(
            context,
            seed_events=(),
        )
    )
    merged_self_expanded_frontier_events = (
        select_coupling_trace_loop_self_deciding_frontier_expanded_events(
            context,
            seed_events=seed_frontier_events,
        )
    )
    self_expansion_event_ids_by_row: dict[str, list[str]] = {}
    for event in self_expanded_frontier_events:
        self_expansion_event_ids_by_row.setdefault(event.row_id, []).append(event.event_id)
    merged_expansion_event_ids_by_row: dict[str, list[str]] = {}
    for event in merged_self_expanded_frontier_events:
        merged_expansion_event_ids_by_row.setdefault(event.row_id, []).append(event.event_id)

    hard_target_rescue_probe: dict[str, dict[str, object]] = {}
    for target_accession in ("4AKE:A", "1MBN:A"):
        target_row = next((row for row in rows if row.source_accession == target_accession), None)
        if target_row is None:
            continue
        target_events = _selected_frontier_events(
            row=target_row,
            event_ids_by_row=event_ids_by_row,
            event_by_id=context.event_by_id,
        )
        target_reports: dict[str, object] = {}
        for target_strategy, min_pairs, max_pairs in (
            (SELF_DECIDING_STRATEGY_NAME, 0, 0),
            ("frontier_internal_gap_balanced", 1, DEFAULT_BALANCED_PAIRS_PER_EVENT),
            ("ridge_coupling", 1, DEFAULT_PRECISION_MAX_PAIRS_PER_EVENT),
            ("frontier_recall", 8, DEFAULT_RECALL_MAX_PAIRS_PER_EVENT),
        ):
            target_result = collapse_row_event_regions(
                row=target_row,
                events=target_events,
                row_features=features_by_row.get(target_row.row_id, ()),
                row_constraints=constraints_by_row.get(target_row.row_id, ()),
                collapse_strategy=target_strategy,
                min_pairs_per_event=min_pairs,
                max_pairs_per_event=max_pairs,
            )
            target_reports[target_strategy] = _row_result_to_report(target_result)
        for expansion_name, expansion_event_ids_by_row in (
            ("self_deciding_frontier_expansion_only", self_expansion_event_ids_by_row),
            ("self_deciding_frontier_expansion_merged", merged_expansion_event_ids_by_row),
        ):
            expansion_events = _selected_frontier_events(
                row=target_row,
                event_ids_by_row=expansion_event_ids_by_row,
                event_by_id=context.event_by_id,
            )
            expansion_result = collapse_row_event_regions(
                row=target_row,
                events=expansion_events,
                row_features=features_by_row.get(target_row.row_id, ()),
                row_constraints=constraints_by_row.get(target_row.row_id, ()),
                collapse_strategy=SELF_DECIDING_STRATEGY_NAME,
                min_pairs_per_event=0,
                max_pairs_per_event=0,
            )
            target_reports[expansion_name] = _row_result_to_report(expansion_result)
        hard_target_rescue_probe[target_accession] = target_reports

    expansion_rows_for_report = self_deciding_frontier_expansion_rows(
        context,
        seed_events=seed_frontier_events,
    )
    report = {
        "report_kind": "event_region_contact_collapse_diagnostics_v0",
        "collapse_kind": EVENT_REGION_CONTACT_COLLAPSE_KIND,
        "collapse_boundary": EVENT_REGION_CONTACT_COLLAPSE_BOUNDARY,
        "benchmark_file": str(BENCHMARK_FILE.relative_to(REPO_ROOT)),
        "external_coupling_file": str(COUPLING_FILE.relative_to(REPO_ROOT)),
        "frontier_file": str(FRONTIER_FILE.relative_to(REPO_ROOT)),
        "coordinate_truth_used_to_build_constraints": False,
        "native_truth_used_before_collapse_selection": False,
        "coordinate_truth_used_before_collapse_selection": False,
        "native_truth_attached_after_collapse_for_evaluation": True,
        "raw_sequence_exposed": False,
        "strategy_mean_rows": {
            strategy: {
                "row_count": len(rows_for_strategy),
                "mean_collapsed_contact_precision": round(
                    sum(float(row["collapsed_contact_precision"]) for row in rows_for_strategy)
                    / max(1, len(rows_for_strategy)),
                    6,
                ),
                "mean_frontier_native_retention": round(
                    sum(float(row["frontier_native_retention"]) for row in rows_for_strategy)
                    / max(1, len(rows_for_strategy)),
                    6,
                ),
                "mean_collapse_reduction_ratio": round(
                    sum(float(row["collapse_reduction_ratio"]) for row in rows_for_strategy)
                    / max(1, len(rows_for_strategy)),
                    6,
                ),
            }
            for strategy, rows_for_strategy in strategy_reports.items()
        },
        "one_cll_strategy_reports": one_cll_reports,
        "one_cll_self_deciding_report": one_cll_reports.get(SELF_DECIDING_STRATEGY_NAME, {}),
        "one_cll_internal_gap_report": one_cll_reports.get("frontier_internal_gap_balanced", {}),
        "hard_target_rescue_probe": hard_target_rescue_probe,
        "self_deciding_frontier_expansion_rows": expansion_rows_for_report,
        "self_deciding_frontier_expansion_interpretation": (
            "The expansion selector is native-free and accession-agnostic. It is now self-verified by contact collapse: "
            "a candidate low-score frontier region is considered only if the row already contains an accepted broad ridge seed. "
            "Candidate regions are then ranked by native-free self-collapse acceptance score and selected at the largest "
            "internal gap in that row-specific distribution, not by a fixed confidence threshold. Native labels remain audit-only."
        ),
        "one_cll_fixed_budget_probe": one_cll_budget_probe,
        "one_cll_best_long_range_f1_budget": max(
            one_cll_budget_probe,
            key=lambda row: float(row["collapsed_long_range_f1"]),
        ) if one_cll_budget_probe else {},
        "row_reports_by_strategy": strategy_reports,
        "interpretation": (
            "Frontier event success and contact-map success are separate. Collapse reduces each 8x8 event region "
            "to a ranked residue-pair subset, then reports precision/recall instead of treating all region pairs as contacts. "
            "The self-deciding strategy uses only internal score distributions, sequence-inferred phase shape, "
            "long-range candidate geometry, direct external-coupling roots, gap clarity, and degree pressure. "
            "Native labels remain post-selection diagnostics only; no accession-specific or fixed pairs-per-event "
            "budget is used by the main self-deciding collapse path."
        ),
    }
    (OUTPUT_DIR / "contact_collapse_diagnostics_v0.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (OUTPUT_DIR / "1cll_contact_collapse_v0.json").write_text(
        json.dumps(one_cll_reports, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _csv_write(OUTPUT_DIR / "contact_collapse_pairs_v0.csv", pair_rows)
    _csv_write(OUTPUT_DIR / "contact_collapse_event_rows_v0.csv", event_rows)
    _csv_write(
        OUTPUT_DIR / "self_deciding_frontier_expansion_rows_v0.csv",
        expansion_rows_for_report,
    )


if __name__ == "__main__":
    main()
