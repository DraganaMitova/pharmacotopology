from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.folding_closure_state_builder import build_physical_states
from pharmacotopology.folding_contact_law_features import contact_law_feature_rows
from pharmacotopology.folding_future_frustration import (
    FUTURE_FRUSTRATION_GATE_KIND,
    FUTURE_FRUSTRATION_LIMIT,
    assess_future_frustration,
)
from pharmacotopology.folding_nucleus_closure_search import (
    FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
    NucleusClosureEvent,
    accepted_events,
    nucleus_closure_events,
)
from pharmacotopology.folding_nucleus_competition import select_competitive_events
from pharmacotopology.folding_nucleus_decoy_falsification import (
    decoy_distance,
    matched_decoys_for_selected_events,
)
from pharmacotopology.folding_nucleus_graph_selectivity import (
    NUCLEUS_GRAPH_SELECTIVITY_REPORT_KIND,
    select_graph_events,
)
from pharmacotopology.folding_physical_state import (
    PHYSICAL_CLOSURE_STATE_KIND,
    PhysicalClosureState,
)
from pharmacotopology.folding_physical_term_ablation import (
    ABLATION_NAMES,
    PHYSICAL_TERM_ABLATION_KIND,
    PhysicalTermAblationRow,
    ablation_interpretation,
    classify_ablation_rows,
    composite_selection_score,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_state_decoy_compare import (
    PHYSICAL_STATE_DECOY_COMPARE_KIND,
    physical_state_decoy_comparisons,
    real_vs_decoy_physical_enrichment_ratio,
)


ACTIVE_PHYSICAL_SELECTION_REPORT_KIND = "active_physical_selection_benchmark_v1"
ACTIVE_PHYSICAL_SELECTION_CERTIFICATE_KIND = "active_physical_selection_certificate"
ACTIVE_PHYSICAL_SELECTION_SCORE_KIND = "physical_score_plus_decoy_margin_v1"

GRAPH_EVENTS_PER_ROW = 40
PHYSICAL_GATE_SCORE_MIN = 0.50
PHYSICAL_GATE_UNSATISFIED_POLAR_MAX = 0.24
PHYSICAL_GATE_STERIC_MAX = 0.08

SURVIVAL_FALSE_RATE_MAX = 0.50
SURVIVAL_CLUSTER_PRECISION_MIN = 0.06
SURVIVAL_LONG_RANGE_RECALL_MIN = 0.30
SURVIVAL_DECOY_ENRICHMENT_MIN = 1.30

ROOT_OUTPUT_NAMES = (
    "active_physical_selection_report.json",
    "active_physical_selection_selectors.csv",
    "active_physical_selection_selected_events.csv",
    "active_physical_selection_ablation.csv",
    "active_physical_selection_dashboard.html",
    "active_physical_selection_certificate.json",
)


@dataclass(frozen=True)
class ActivePhysicalContext:
    rows: tuple[RealCoordinateVisualRow, ...]
    candidate_events: tuple[NucleusClosureEvent, ...]
    competitive_events: tuple[NucleusClosureEvent, ...]
    graph_events: tuple[NucleusClosureEvent, ...]
    states: tuple[PhysicalClosureState, ...]
    state_by_event_id: Mapping[str, PhysicalClosureState]
    event_by_id: Mapping[str, NucleusClosureEvent]
    decoy_margin_by_event_id: Mapping[str, float]


@dataclass(frozen=True)
class ActiveSelectionMetric:
    selector_name: str
    selected_event_count: int
    false_nucleus_rate: float
    contact_cluster_precision: float
    long_range_contact_recall: float
    real_vs_decoy_physical_enrichment_ratio: float
    mean_active_physical_score: float
    survives_targets: bool
    native_truth_used_before_active_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(value, 6)


def _region_union(events: Sequence[NucleusClosureEvent]) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    for event in events:
        pairs.update(event.candidate_region_pairs())
    return pairs


def build_active_physical_context(
    rows: Sequence[RealCoordinateVisualRow],
) -> ActivePhysicalContext:
    feature_rows = contact_law_feature_rows(rows)
    candidate_events = nucleus_closure_events(rows, feature_rows)
    accepted = accepted_events(candidate_events, threshold=0.30)
    competitive_events, _ = select_competitive_events(rows, accepted)
    graph_events, _ = select_graph_events(rows, competitive_events)
    states, _ = build_physical_states(rows=rows, events=competitive_events)
    state_by_event = {state.event_id: state for state in states}
    event_by_id = {event.event_id: event for event in competitive_events}
    by_row: dict[str, list[NucleusClosureEvent]] = defaultdict(list)
    for event in competitive_events:
        by_row[event.row_id].append(event)
    margins: dict[str, float] = {}
    for row_events in by_row.values():
        for event in row_events:
            decoy_candidates = [
                candidate
                for candidate in row_events
                if candidate.event_id != event.event_id
            ]
            if not decoy_candidates:
                margins[event.event_id] = 0.0
                continue
            decoy = min(decoy_candidates, key=lambda item: decoy_distance(event, item))
            margins[event.event_id] = _rounded(
                state_by_event[event.event_id].physical_state_score
                - state_by_event[decoy.event_id].physical_state_score
            )
    return ActivePhysicalContext(
        rows=tuple(rows),
        candidate_events=candidate_events,
        competitive_events=competitive_events,
        graph_events=graph_events,
        states=states,
        state_by_event_id=state_by_event,
        event_by_id=event_by_id,
        decoy_margin_by_event_id=margins,
    )


def active_physical_score(
    event: NucleusClosureEvent,
    context: ActivePhysicalContext,
    *,
    remove_term: str | None = None,
) -> float:
    state = context.state_by_event_id[event.event_id]
    score = 0.0
    if remove_term != "without_burial_gain":
        score += state.burial_gain
    if remove_term != "without_loop_strain":
        score -= 0.75 * state.loop_strain
    if remove_term != "without_steric_clash":
        score -= 0.85 * state.steric_clash_score
    if remove_term != "without_unsatisfied_polar_penalty":
        score -= 0.75 * state.unsatisfied_polar_penalty
    if remove_term != "without_future_frustration":
        score -= 0.65 * state.future_frustration_score
    if remove_term != "without_decoy_margin":
        score += 0.20 * context.decoy_margin_by_event_id[event.event_id]
    return _rounded(score)


def _events_by_row(events: Sequence[NucleusClosureEvent]) -> dict[str, list[NucleusClosureEvent]]:
    output: dict[str, list[NucleusClosureEvent]] = defaultdict(list)
    for event in events:
        output[event.row_id].append(event)
    return output


def select_events(
    context: ActivePhysicalContext,
    *,
    selector_name: str,
    remove_term: str | None = None,
) -> tuple[NucleusClosureEvent, ...]:
    competitive_by_row = _events_by_row(context.competitive_events)
    graph_by_row = _events_by_row(context.graph_events)
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        if selector_name == "graph_only":
            selected.extend(graph_by_row.get(row.row_id, ()))
            continue
        candidates = competitive_by_row.get(row.row_id, ())
        if selector_name in {"physical_gate", "future_frustration_gate"}:
            candidates = tuple(
                event
                for event in candidates
                if active_physical_score(event, context, remove_term=remove_term)
                >= PHYSICAL_GATE_SCORE_MIN
                and context.state_by_event_id[
                    event.event_id
                ].unsatisfied_polar_penalty
                <= PHYSICAL_GATE_UNSATISFIED_POLAR_MAX
                and context.state_by_event_id[event.event_id].steric_clash_score
                <= PHYSICAL_GATE_STERIC_MAX
            )
        if selector_name == "future_frustration_gate":
            candidates = tuple(
                event
                for event in candidates
                if assess_future_frustration(
                    context.state_by_event_id[event.event_id],
                    limit=FUTURE_FRUSTRATION_LIMIT,
                ).future_closure_path_preserved
            )
        selected.extend(
            sorted(
                candidates,
                key=lambda event: (
                    -active_physical_score(event, context, remove_term=remove_term),
                    event.segment_a_start,
                    event.segment_b_start,
                    event.event_id,
                ),
            )[:GRAPH_EVENTS_PER_ROW]
        )
    return tuple(selected)


def selector_metrics(
    context: ActivePhysicalContext,
    *,
    selector_name: str,
    selected_events: Sequence[NucleusClosureEvent],
) -> ActiveSelectionMetric:
    by_row = _events_by_row(selected_events)
    false_rates: list[float] = []
    precisions: list[float] = []
    long_recalls: list[float] = []
    for row in context.rows:
        row_events = tuple(by_row.get(row.row_id, ()))
        native_pairs = set(row.native_contact_pairs())
        native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
        region = _region_union(row_events)
        native_hit_count = sum(
            event.native_contact_count_after_scoring for event in row_events
        )
        possible_region_pair_count = len(row_events) * 64
        false_count = sum(
            1 for event in row_events if event.native_contact_count_after_scoring == 0
        )
        false_rates.append(false_count / len(row_events) if row_events else 0.0)
        precisions.append(
            native_hit_count / possible_region_pair_count
            if possible_region_pair_count
            else 0.0
        )
        long_recalls.append(
            len(region & native_long) / len(native_long) if native_long else 1.0
        )
    matches = matched_decoys_for_selected_events(
        selected_events=selected_events,
        candidate_events=context.competitive_events,
    )
    decoy_events = tuple(context.event_by_id[match.decoy_event_id] for match in matches)
    comparisons = physical_state_decoy_comparisons(
        matches=matches,
        real_states=tuple(
            context.state_by_event_id[event.event_id] for event in selected_events
        ),
        decoy_states=tuple(
            context.state_by_event_id[event.event_id] for event in decoy_events
        ),
    )
    false_rate = _rounded(mean(false_rates) if false_rates else 0.0)
    precision = _rounded(mean(precisions) if precisions else 0.0)
    long_recall = _rounded(mean(long_recalls) if long_recalls else 0.0)
    enrichment = real_vs_decoy_physical_enrichment_ratio(comparisons)
    survives = (
        false_rate < SURVIVAL_FALSE_RATE_MAX
        and precision > SURVIVAL_CLUSTER_PRECISION_MIN
        and long_recall > SURVIVAL_LONG_RANGE_RECALL_MIN
        and enrichment > SURVIVAL_DECOY_ENRICHMENT_MIN
    )
    return ActiveSelectionMetric(
        selector_name=selector_name,
        selected_event_count=len(selected_events),
        false_nucleus_rate=false_rate,
        contact_cluster_precision=precision,
        long_range_contact_recall=long_recall,
        real_vs_decoy_physical_enrichment_ratio=enrichment,
        mean_active_physical_score=_rounded(
            mean(
                [
                    active_physical_score(event, context)
                    for event in selected_events
                ]
            )
            if selected_events
            else 0.0
        ),
        survives_targets=survives,
    )


def selected_event_rows(
    context: ActivePhysicalContext,
    selections: Mapping[str, Sequence[NucleusClosureEvent]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for selector_name in sorted(selections):
        for rank, event in enumerate(
            sorted(
                selections[selector_name],
                key=lambda item: (
                    item.row_id,
                    -active_physical_score(item, context),
                    item.segment_a_start,
                    item.segment_b_start,
                ),
            ),
            start=1,
        ):
            state = context.state_by_event_id[event.event_id]
            rows.append(
                {
                    "selector_name": selector_name,
                    "rank": rank,
                    "row_id": event.row_id,
                    "source_accession": event.source_accession,
                    "event_id": event.event_id,
                    "active_physical_score": active_physical_score(event, context),
                    "physical_state_score": state.physical_state_score,
                    "decoy_margin": context.decoy_margin_by_event_id[event.event_id],
                    "loop_strain": state.loop_strain,
                    "steric_clash_score": state.steric_clash_score,
                    "burial_gain": state.burial_gain,
                    "unsatisfied_polar_penalty": state.unsatisfied_polar_penalty,
                    "future_frustration_score": state.future_frustration_score,
                    "native_contact_count_after_scoring": (
                        event.native_contact_count_after_scoring
                    ),
                    "native_truth_used_before_active_selection": False,
                    "raw_sequence_exposed": False,
                }
            )
    return rows


def ablation_rows(context: ActivePhysicalContext) -> tuple[PhysicalTermAblationRow, ...]:
    full_selected = select_events(context, selector_name="physical_rerank")
    full_metrics = selector_metrics(
        context,
        selector_name="physical_rerank",
        selected_events=full_selected,
    )
    full_composite = composite_selection_score(
        false_nucleus_rate=full_metrics.false_nucleus_rate,
        contact_cluster_precision=full_metrics.contact_cluster_precision,
        long_range_contact_recall=full_metrics.long_range_contact_recall,
        real_vs_decoy_physical_enrichment_ratio=(
            full_metrics.real_vs_decoy_physical_enrichment_ratio
        ),
    )
    rows: list[PhysicalTermAblationRow] = []
    for ablation_name in ABLATION_NAMES:
        selected = select_events(
            context,
            selector_name="physical_rerank",
            remove_term=ablation_name,
        )
        metrics = selector_metrics(
            context,
            selector_name=ablation_name,
            selected_events=selected,
        )
        composite = composite_selection_score(
            false_nucleus_rate=metrics.false_nucleus_rate,
            contact_cluster_precision=metrics.contact_cluster_precision,
            long_range_contact_recall=metrics.long_range_contact_recall,
            real_vs_decoy_physical_enrichment_ratio=(
                metrics.real_vs_decoy_physical_enrichment_ratio
            ),
        )
        delta = _rounded(composite - full_composite)
        rows.append(
            PhysicalTermAblationRow(
                ablation_name=ablation_name,
                false_nucleus_rate=metrics.false_nucleus_rate,
                contact_cluster_precision=metrics.contact_cluster_precision,
                long_range_contact_recall=metrics.long_range_contact_recall,
                real_vs_decoy_physical_enrichment_ratio=(
                    metrics.real_vs_decoy_physical_enrichment_ratio
                ),
                composite_selection_score=composite,
                delta_vs_full_selector=delta,
                term_interpretation=ablation_interpretation(delta),
            )
        )
    return tuple(rows)


def build_active_physical_selection_report(
    *,
    context: ActivePhysicalContext,
    selector_rows: Sequence[ActiveSelectionMetric],
    ablations: Sequence[PhysicalTermAblationRow],
    source_benchmark_file: Path,
) -> dict[str, object]:
    selector_lookup = {row.selector_name: row for row in selector_rows}
    term_summary = classify_ablation_rows(ablations)
    active_survives = any(row.survives_targets for row in selector_rows)
    return {
        "report_kind": ACTIVE_PHYSICAL_SELECTION_REPORT_KIND,
        "source_benchmark_kind": REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "source_event_kind": FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
        "source_graph_report_kind": NUCLEUS_GRAPH_SELECTIVITY_REPORT_KIND,
        "physical_state_kind": PHYSICAL_CLOSURE_STATE_KIND,
        "future_frustration_gate_kind": FUTURE_FRUSTRATION_GATE_KIND,
        "physical_state_decoy_compare_kind": PHYSICAL_STATE_DECOY_COMPARE_KIND,
        "physical_term_ablation_kind": PHYSICAL_TERM_ABLATION_KIND,
        "active_physical_selection_score_kind": ACTIVE_PHYSICAL_SELECTION_SCORE_KIND,
        "benchmark_size": len(context.rows),
        "candidate_event_count": len(context.competitive_events),
        "graph_only_false_nucleus_rate": selector_lookup[
            "graph_only"
        ].false_nucleus_rate,
        "physical_rerank_false_nucleus_rate": selector_lookup[
            "physical_rerank"
        ].false_nucleus_rate,
        "physical_gate_false_nucleus_rate": selector_lookup[
            "physical_gate"
        ].false_nucleus_rate,
        "future_frustration_false_nucleus_rate": selector_lookup[
            "future_frustration_gate"
        ].false_nucleus_rate,
        "graph_only_cluster_precision": selector_lookup[
            "graph_only"
        ].contact_cluster_precision,
        "physical_rerank_cluster_precision": selector_lookup[
            "physical_rerank"
        ].contact_cluster_precision,
        "physical_gate_cluster_precision": selector_lookup[
            "physical_gate"
        ].contact_cluster_precision,
        "future_frustration_cluster_precision": selector_lookup[
            "future_frustration_gate"
        ].contact_cluster_precision,
        "graph_only_long_range_recall": selector_lookup[
            "graph_only"
        ].long_range_contact_recall,
        "physical_rerank_long_range_recall": selector_lookup[
            "physical_rerank"
        ].long_range_contact_recall,
        "physical_gate_long_range_recall": selector_lookup[
            "physical_gate"
        ].long_range_contact_recall,
        "future_frustration_long_range_recall": selector_lookup[
            "future_frustration_gate"
        ].long_range_contact_recall,
        "physical_rerank_real_vs_decoy_enrichment_ratio": selector_lookup[
            "physical_rerank"
        ].real_vs_decoy_physical_enrichment_ratio,
        "active_physical_selection_survives": active_survives,
        **term_summary,
        "native_truth_used_before_active_selection": False,
        "native_truth_used_before_physical_scoring": False,
        "native_label_attached_after_active_selection": True,
        "mechanism_discovery_claim_allowed": False,
        "mechanism_discovery_claim_created": False,
        "global_folding_claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "claim_allowed": False,
        "boundary_statement": (
            "Active physical selection can improve decoy enrichment and "
            "long-range recall under reranking, while viability gates reduce "
            "false nuclei at the cost of coverage. No selector passes all "
            "survival gates."
        ),
        "selectors": [row.to_dict() for row in selector_rows],
        "ablations": [row.to_dict() for row in ablations],
    }


def build_active_physical_selection_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": ACTIVE_PHYSICAL_SELECTION_CERTIFICATE_KIND,
        "report_kind": report["report_kind"],
        "active_physical_selection_survives": report[
            "active_physical_selection_survives"
        ],
        "physical_rerank_false_nucleus_rate": report[
            "physical_rerank_false_nucleus_rate"
        ],
        "physical_rerank_cluster_precision": report[
            "physical_rerank_cluster_precision"
        ],
        "physical_rerank_long_range_recall": report[
            "physical_rerank_long_range_recall"
        ],
        "physical_rerank_real_vs_decoy_enrichment_ratio": report[
            "physical_rerank_real_vs_decoy_enrichment_ratio"
        ],
        "best_physical_term": report["best_physical_term"],
        "worst_physical_term": report["worst_physical_term"],
        "native_truth_used_before_active_selection": report[
            "native_truth_used_before_active_selection"
        ],
        "mechanism_discovery_claim_allowed": report[
            "mechanism_discovery_claim_allowed"
        ],
        "global_folding_claim_allowed": report["global_folding_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "output_artifacts": tuple(output_names),
    }


def write_active_physical_selection_outputs(
    *,
    report: Mapping[str, object],
    selector_rows: Sequence[ActiveSelectionMetric],
    selected_rows: Sequence[Mapping[str, object]],
    ablations: Sequence[PhysicalTermAblationRow],
    report_path: Path,
    selectors_path: Path,
    selected_events_path: Path,
    ablation_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows([row.to_dict() for row in selector_rows], selectors_path)
    _write_csv_rows(selected_rows, selected_events_path)
    _write_csv_rows([row.to_dict() for row in ablations], ablation_path)
    dashboard_path.write_text(
        render_active_physical_selection_dashboard(report),
        encoding="utf-8",
    )
    certificate = build_active_physical_selection_certificate(
        report,
        output_names=ROOT_OUTPUT_NAMES,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        selectors_path,
        selected_events_path,
        ablation_path,
        dashboard_path,
        certificate_path,
    )


def run_active_physical_selection_benchmark(
    *,
    benchmark_file: Path,
    report_path: Path,
    selectors_path: Path,
    selected_events_path: Path,
    ablation_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    context = build_active_physical_context(rows)
    selections = {
        "graph_only": select_events(context, selector_name="graph_only"),
        "physical_rerank": select_events(context, selector_name="physical_rerank"),
        "physical_gate": select_events(context, selector_name="physical_gate"),
        "future_frustration_gate": select_events(
            context,
            selector_name="future_frustration_gate",
        ),
    }
    selector_rows = tuple(
        selector_metrics(
            context,
            selector_name=name,
            selected_events=events,
        )
        for name, events in selections.items()
    )
    ablations = ablation_rows(context)
    selected_rows = selected_event_rows(context, selections)
    report = build_active_physical_selection_report(
        context=context,
        selector_rows=selector_rows,
        ablations=ablations,
        source_benchmark_file=benchmark_file,
    )
    return write_active_physical_selection_outputs(
        report=report,
        selector_rows=selector_rows,
        selected_rows=selected_rows,
        ablations=ablations,
        report_path=report_path,
        selectors_path=selectors_path,
        selected_events_path=selected_events_path,
        ablation_path=ablation_path,
        dashboard_path=dashboard_path,
        certificate_path=certificate_path,
    )


def _write_csv_rows(rows: Sequence[Mapping[str, object]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        if not rows:
            return path
        fieldnames = list(rows[0])
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    return path


def _escape(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _metric_cards(report: Mapping[str, object]) -> str:
    labels = (
        "physical_rerank_false_nucleus_rate",
        "physical_rerank_cluster_precision",
        "physical_rerank_long_range_recall",
        "physical_rerank_real_vs_decoy_enrichment_ratio",
        "physical_gate_false_nucleus_rate",
        "physical_gate_cluster_precision",
        "physical_gate_long_range_recall",
        "future_frustration_false_nucleus_rate",
        "best_physical_term",
        "worst_physical_term",
        "active_physical_selection_survives",
        "folding_problem_solved",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _selector_table(report: Mapping[str, object]) -> str:
    rows = report.get("selectors", [])
    if not isinstance(rows, list):
        rows = []
    body = "".join(
        "<tr>"
        f"<td>{_escape(row.get('selector_name'))}</td>"
        f"<td>{_escape(row.get('selected_event_count'))}</td>"
        f"<td>{_escape(row.get('false_nucleus_rate'))}</td>"
        f"<td>{_escape(row.get('contact_cluster_precision'))}</td>"
        f"<td>{_escape(row.get('long_range_contact_recall'))}</td>"
        f"<td>{_escape(row.get('real_vs_decoy_physical_enrichment_ratio'))}</td>"
        "</tr>"
        for row in rows
        if isinstance(row, Mapping)
    )
    return (
        "<section><h2>Selector Comparison</h2>"
        "<table><thead><tr><th>selector</th><th>events</th><th>false</th>"
        "<th>precision</th><th>long-range</th><th>decoy ratio</th></tr></thead>"
        "<tbody>"
        + body
        + "</tbody></table></section>"
    )


def render_active_physical_selection_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Active Physical Selection</title>
  <style>
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f4f6f3;
      color: #202522;
    }}
    header {{
      padding: 34px;
      background: #26332f;
      color: #f7f8f4;
    }}
    main {{
      max-width: 1240px;
      margin: 0 auto;
      padding: 24px;
    }}
    h1, h2, h3 {{
      margin: 0 0 12px;
      letter-spacing: 0;
    }}
    section {{
      margin: 24px 0;
    }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
      gap: 12px;
    }}
    .metric {{
      background: #ffffff;
      border: 1px solid #d3ddd5;
      border-radius: 6px;
      padding: 14px;
    }}
    .metric span {{
      display: block;
      color: #5b6660;
      font-size: 12px;
      overflow-wrap: anywhere;
    }}
    .metric strong {{
      display: block;
      margin-top: 8px;
      font-size: 20px;
      overflow-wrap: anywhere;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: #ffffff;
      border: 1px solid #d3ddd5;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid #e2e8e3;
      text-align: left;
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Active Physical Selection</h1>
    <p>Physical state terms now control closure selection, then the selected nuclei are scored after selection.</p>
  </header>
  <main>
    <section class="metrics">{_metric_cards(report)}</section>
    {_selector_table(report)}
    <section>
      <h2>Claim Boundary</h2>
      <p>{_escape(report.get("boundary_statement", ""))}</p>
    </section>
  </main>
</body>
</html>
"""

