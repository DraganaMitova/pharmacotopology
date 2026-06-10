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
    select_coupling_trace_loop_self_deciding_frontier_generated_events,
    self_deciding_frontier_expansion_rows,
    self_deciding_frontier_generation_rows,
)
from pharmacotopology.folding_event_region_contact_collapse import (
    DEFAULT_BALANCED_PAIRS_PER_EVENT,
    DEFAULT_PRECISION_MAX_PAIRS_PER_EVENT,
    DEFAULT_RECALL_MAX_PAIRS_PER_EVENT,
    EVENT_REGION_CONTACT_COLLAPSE_BOUNDARY,
    EVENT_REGION_CONTACT_COLLAPSE_KIND,
    SELF_DECIDING_STRATEGY_NAME,
    collapse_event_region_contacts,
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


def _long_native_pairs(row: RealCoordinateVisualRow) -> set[tuple[int, int]]:
    return {pair for pair in row.native_contact_pairs() if abs(pair[1] - pair[0]) >= 24}


def _event_region_long_recall_ceiling(
    *,
    row: RealCoordinateVisualRow,
    events: Sequence[NucleusClosureEvent],
) -> float:
    native_long = _long_native_pairs(row)
    if not native_long:
        return 0.0
    region_pairs: set[tuple[int, int]] = set()
    for event in events:
        region_pairs.update(event.candidate_region_pairs())
    return round(len(region_pairs & native_long) / len(native_long), 6)


def _self_deciding_contact_long_recall_ceiling(
    *,
    row: RealCoordinateVisualRow,
    events: Sequence[NucleusClosureEvent],
    features_by_row: Mapping[str, Sequence[object]],
    constraints_by_row: Mapping[str, Sequence[object]],
) -> dict[str, object]:
    native_pairs = set(row.native_contact_pairs())
    native_long = _long_native_pairs(row)
    predicted: set[tuple[int, int]] = set()
    for event in events:
        pairs, _summary = collapse_event_region_contacts(
            event=event,
            row_features=features_by_row.get(row.row_id, ()),  # type: ignore[arg-type]
            row_constraints=constraints_by_row.get(row.row_id, ()),  # type: ignore[arg-type]
            collapse_strategy=SELF_DECIDING_STRATEGY_NAME,
            min_pairs_per_event=0,
            max_pairs_per_event=0,
        )
        predicted.update(pair.pair() for pair in pairs)
    tp = predicted & native_pairs
    long_tp = predicted & native_long
    return {
        "self_deciding_all_candidate_pair_count": len(predicted),
        "self_deciding_all_candidate_true_positive_contacts": len(tp),
        "self_deciding_all_candidate_precision": round(len(tp) / len(predicted), 6) if predicted else 0.0,
        "self_deciding_all_candidate_long_range_recall": round(len(long_tp) / len(native_long), 6) if native_long else 0.0,
        "native_truth_used_before_frontier_ceiling_audit": False,
        "native_truth_attached_after_frontier_ceiling_for_evaluation": True,
    }


def _frontier_ceiling_audit(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    context,
    event_ids_by_row: Mapping[str, Sequence[str]],
    merged_expansion_event_ids_by_row: Mapping[str, Sequence[str]],
    generated_event_ids_by_row: Mapping[str, Sequence[str]],
    features_by_row: Mapping[str, Sequence[object]],
    constraints_by_row: Mapping[str, Sequence[object]],
) -> dict[str, dict[str, object]]:
    output: dict[str, dict[str, object]] = {}
    competitive_by_row: dict[str, list[NucleusClosureEvent]] = {}
    candidate_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in context.competitive_events:
        competitive_by_row.setdefault(event.row_id, []).append(event)
    for event in context.physical_context.candidate_events:
        candidate_by_row.setdefault(event.row_id, []).append(event)
    for row in rows:
        if row.source_accession not in {"1CLL:A", "1MBN:A", "4AKE:A"}:
            continue
        seed_events = _selected_frontier_events(
            row=row,
            event_ids_by_row=event_ids_by_row,
            event_by_id=context.event_by_id,
        )
        expanded_events = _selected_frontier_events(
            row=row,
            event_ids_by_row=merged_expansion_event_ids_by_row,
            event_by_id=context.event_by_id,
        )
        generated_events = _selected_frontier_events(
            row=row,
            event_ids_by_row=generated_event_ids_by_row,
            event_by_id=context.event_by_id,
        )
        competitive_events = tuple(competitive_by_row.get(row.row_id, ()))
        candidate_events = tuple(candidate_by_row.get(row.row_id, ()))
        contact_ceiling = _self_deciding_contact_long_recall_ceiling(
            row=row,
            events=competitive_events,
            features_by_row=features_by_row,
            constraints_by_row=constraints_by_row,
        )
        output[row.source_accession] = {
            "native_long_range_contact_count": len(_long_native_pairs(row)),
            "seed_frontier_event_count": len(seed_events),
            "expanded_frontier_event_count": len(expanded_events),
            "competitive_event_count": len(competitive_events),
            "candidate_event_count": len(candidate_events),
            "generated_frontier_event_count": len(generated_events),
            "seed_region_long_range_recall_ceiling": _event_region_long_recall_ceiling(row=row, events=seed_events),
            "expanded_region_long_range_recall_ceiling": _event_region_long_recall_ceiling(row=row, events=expanded_events),
            "generated_region_long_range_recall_ceiling": _event_region_long_recall_ceiling(row=row, events=generated_events),
            "competitive_region_long_range_recall_ceiling": _event_region_long_recall_ceiling(row=row, events=competitive_events),
            "candidate_generator_long_range_recall_ceiling": _event_region_long_recall_ceiling(row=row, events=candidate_events),
            "frontier_generation_bottleneck_for_0_40_recall": (
                _event_region_long_recall_ceiling(row=row, events=competitive_events) < 0.40
            ),
            **contact_ceiling,
        }
    return output


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

    generated_frontier_events = (
        select_coupling_trace_loop_self_deciding_frontier_generated_events(
            context,
            seed_events=seed_frontier_events,
        )
    )
    generated_event_ids_by_row: dict[str, list[str]] = {}
    for event in generated_frontier_events:
        generated_event_ids_by_row.setdefault(event.row_id, []).append(event.event_id)

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
            ("self_deciding_frontier_generation_merged", generated_event_ids_by_row),
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

    frontier_ceiling_audit = _frontier_ceiling_audit(
        rows=rows,
        context=context,
        event_ids_by_row=event_ids_by_row,
        merged_expansion_event_ids_by_row=merged_expansion_event_ids_by_row,
        generated_event_ids_by_row=generated_event_ids_by_row,
        features_by_row=features_by_row,
        constraints_by_row=constraints_by_row,
    )
    expansion_rows_for_report = self_deciding_frontier_expansion_rows(
        context,
        seed_events=seed_frontier_events,
    )
    generation_rows_for_report = self_deciding_frontier_generation_rows(
        context,
        seed_events=seed_frontier_events,
    )

    frontier_generation_decisions: dict[str, dict[str, object]] = {}
    for accession, probes in hard_target_rescue_probe.items():
        expanded = probes.get("self_deciding_frontier_expansion_merged", {})
        generated = probes.get("self_deciding_frontier_generation_merged", {})
        expanded_precision = float(expanded.get("collapsed_contact_precision", 0.0)) if isinstance(expanded, Mapping) else 0.0
        generated_precision = float(generated.get("collapsed_contact_precision", 0.0)) if isinstance(generated, Mapping) else 0.0
        generated_recall = float(generated.get("collapsed_long_range_recall", 0.0)) if isinstance(generated, Mapping) else 0.0
        expanded_recall = float(expanded.get("collapsed_long_range_recall", 0.0)) if isinstance(expanded, Mapping) else 0.0
        improves = (
            generated_precision >= expanded_precision
            and generated_recall >= expanded_recall
            and (generated_precision > expanded_precision or generated_recall > expanded_recall)
        )
        equivalent = generated_precision == expanded_precision and generated_recall == expanded_recall
        frontier_generation_decisions[accession] = {
            "frontier_generation_probe_accepted_as_main": improves,
            "frontier_generation_decision": (
                "accepted_generation_improves_precision_or_recall_without_regression"
                if improves
                else "equivalent_no_main_change"
                if equivalent
                else "rejected_precision_collapse_or_no_recall_gain"
            ),
            "expanded_precision_reference": expanded_precision,
            "generated_precision": generated_precision,
            "expanded_long_range_recall_reference": expanded_recall,
            "generated_long_range_recall": generated_recall,
            "native_truth_used_before_frontier_generation_decision": False,
            "native_truth_attached_after_frontier_generation_for_evaluation": True,
        }

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
        "frontier_ceiling_audit": frontier_ceiling_audit,
        "self_deciding_frontier_expansion_rows": expansion_rows_for_report,
        "self_deciding_frontier_generation_rows": generation_rows_for_report,
        "frontier_generation_decisions": frontier_generation_decisions,
        "self_deciding_frontier_generation_interpretation": (
            "A broader candidate-generator frontier was implemented as an honest probe. It ranks candidate events from the full "
            "candidate pool by row-local identity-normalized native-free evidence, applies an internal-gap prefilter, collapses "
            "surviving regions, then applies another internal-gap boundary. If the generated path increases recall but destroys "
            "precision, the report rejects it as the main path instead of calling the target solved."
        ),
        "self_deciding_frontier_expansion_interpretation": (
            "The expansion selector is native-free and accession-agnostic. It is now self-verified by contact collapse: "
            "a candidate low-score frontier region is considered only if the row already contains an accepted broad ridge seed. "
            "Candidate regions are normalized against the seed frontier's own identity baseline, then selected at "
            "the largest internal gap plus any row-local seed-envelope tier. This is not a fixed confidence threshold. "
            "Native labels remain audit-only; the report also records the frontier-ceiling audit so recall limits are not hidden."
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
    _csv_write(
        OUTPUT_DIR / "self_deciding_frontier_generation_rows_v0.csv",
        generation_rows_for_report,
    )


if __name__ == "__main__":
    main()
