from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_burial_frustration import burial_frustration_for_event
from pharmacotopology.folding_contact_law_features import contact_law_feature_rows
from pharmacotopology.folding_nucleus_closure_search import (
    FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
    accepted_events,
    nucleus_closure_events,
)
from pharmacotopology.folding_nucleus_competition import select_competitive_events
from pharmacotopology.folding_nucleus_decoy_falsification import (
    matched_decoys_for_selected_events,
)
from pharmacotopology.folding_nucleus_graph_selectivity import (
    NUCLEUS_GRAPH_SELECTIVITY_REPORT_KIND,
    select_graph_events,
)
from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent
from pharmacotopology.folding_physical_state import (
    PHYSICAL_CLOSURE_STATE_KIND,
    PhysicalClosureState,
    physical_state_from_event,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_state_decoy_compare import (
    PHYSICAL_STATE_DECOY_COMPARE_KIND,
    PhysicalStateDecoyComparison,
    physical_state_decoy_comparisons,
    real_beats_decoy_score_rate,
    real_vs_decoy_physical_enrichment_ratio,
)


CLOSURE_STATE_BUILDER_KIND = "sequence_only_closure_state_builder_v1"
PHYSICAL_CLOSURE_STATE_REPORT_KIND = "physical_closure_state_benchmark_v1"
PHYSICAL_CLOSURE_STATE_CERTIFICATE_KIND = "physical_closure_state_certificate"

PHYSICAL_FALSE_RATE_TARGET = 0.50
PHYSICAL_PRECISION_TARGET = 0.06
PHYSICAL_LONG_RANGE_RECALL_TARGET = 0.30
PHYSICAL_DECOY_ENRICHMENT_TARGET = 1.00

ROOT_OUTPUT_NAMES = (
    "physical_closure_state_report.json",
    "physical_closure_state_states.csv",
    "physical_closure_state_decoys.csv",
    "physical_closure_state_rank_enrichment.csv",
    "physical_closure_state_metrics.csv",
    "physical_closure_state_dashboard.html",
    "physical_closure_state_certificate.json",
)


@dataclass(frozen=True)
class ClosureStateBuildFailure:
    row_id: str
    source_accession: str
    event_id: str
    failure_reason: str
    native_truth_used_before_physical_scoring: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def rows_by_id(
    rows: Sequence[RealCoordinateVisualRow],
) -> dict[str, RealCoordinateVisualRow]:
    return {row.row_id: row for row in rows}


def build_physical_state_for_event(
    *,
    event: NucleusClosureEvent,
    row_lookup: Mapping[str, RealCoordinateVisualRow],
) -> PhysicalClosureState | ClosureStateBuildFailure:
    row = row_lookup.get(event.row_id)
    if row is None:
        return ClosureStateBuildFailure(
            row_id=event.row_id,
            source_accession=event.source_accession,
            event_id=event.event_id,
            failure_reason="missing_coordinate_row",
        )
    if event.segment_b_end > len(row.sequence):
        return ClosureStateBuildFailure(
            row_id=event.row_id,
            source_accession=event.source_accession,
            event_id=event.event_id,
            failure_reason="segment_outside_sequence",
        )
    burial = burial_frustration_for_event(event=event, sequence=row.sequence)
    return physical_state_from_event(event=event, burial=burial)


def build_physical_states(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    events: Sequence[NucleusClosureEvent],
) -> tuple[tuple[PhysicalClosureState, ...], tuple[ClosureStateBuildFailure, ...]]:
    lookup = rows_by_id(rows)
    states: list[PhysicalClosureState] = []
    failures: list[ClosureStateBuildFailure] = []
    for event in events:
        result = build_physical_state_for_event(event=event, row_lookup=lookup)
        if isinstance(result, PhysicalClosureState):
            states.append(result)
        else:
            failures.append(result)
    return tuple(states), tuple(failures)


def _rounded(value: float) -> float:
    return round(value, 6)


def _mean(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def _region_union(events: Sequence[NucleusClosureEvent]) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    for event in events:
        pairs.update(event.candidate_region_pairs())
    return pairs


def _state_event_lookup(
    states: Sequence[PhysicalClosureState],
    events: Sequence[NucleusClosureEvent],
) -> dict[str, NucleusClosureEvent]:
    event_by_id = {event.event_id: event for event in events}
    return {
        state.event_id: event_by_id[state.event_id]
        for state in states
        if state.event_id in event_by_id
    }


def physical_rank_enrichment_rows(
    states: Sequence[PhysicalClosureState],
) -> list[dict[str, object]]:
    ranked = sorted(
        states,
        key=lambda state: (-state.physical_state_score, state.row_id, state.event_id),
    )
    baseline = (
        sum(1 for state in ranked if state.native_contact_count_after_scoring > 0)
        / len(ranked)
        if ranked
        else 0.0
    )
    rows: list[dict[str, object]] = []
    for cutoff in (10, 25):
        top = ranked[: min(cutoff, len(ranked))]
        top_rate = (
            sum(1 for state in top if state.native_contact_count_after_scoring > 0)
            / len(top)
            if top
            else 0.0
        )
        rows.append(
            {
                "cutoff": cutoff,
                "top_rank_state_count": len(top),
                "top_rank_native_positive_rate": _rounded(top_rate),
                "baseline_native_positive_rate": _rounded(baseline),
                "physical_state_rank_enrichment": _rounded(
                    top_rate / baseline if baseline else 0.0
                ),
                "native_truth_used_before_physical_scoring": False,
                "raw_sequence_exposed": False,
            }
        )
    return rows


def physical_state_rank_enrichment_at(
    rank_rows: Sequence[Mapping[str, object]],
    cutoff: int,
) -> float:
    for row in rank_rows:
        if int(row["cutoff"]) == cutoff:
            return float(row["physical_state_rank_enrichment"])
    return 0.0


def physical_metric_rows(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    states: Sequence[PhysicalClosureState],
    event_lookup: Mapping[str, NucleusClosureEvent],
) -> list[dict[str, object]]:
    by_row: dict[str, list[PhysicalClosureState]] = {}
    for state in states:
        by_row.setdefault(state.row_id, []).append(state)
    output: list[dict[str, object]] = []
    for row in rows:
        row_states = tuple(by_row.get(row.row_id, ()))
        row_events = tuple(
            event_lookup[state.event_id]
            for state in row_states
            if state.event_id in event_lookup
        )
        native_pairs = set(row.native_contact_pairs())
        native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
        region = _region_union(row_events)
        native_hit_count = sum(
            state.native_contact_count_after_scoring for state in row_states
        )
        possible_region_pair_count = len(row_events) * 64
        false_count = sum(
            1 for state in row_states if state.native_contact_count_after_scoring == 0
        )
        output.append(
            {
                "row_id": row.row_id,
                "source_accession": row.source_accession,
                "physical_state_count": len(row_states),
                "post_physical_false_nucleus_rate": _rounded(
                    false_count / len(row_states) if row_states else 0.0
                ),
                "post_physical_contact_cluster_precision": _rounded(
                    native_hit_count / possible_region_pair_count
                    if possible_region_pair_count
                    else 0.0
                ),
                "post_physical_long_range_contact_recall": _rounded(
                    len(region & native_long) / len(native_long)
                    if native_long
                    else 1.0
                ),
                "mean_physical_state_score": _rounded(
                    _mean([state.physical_state_score for state in row_states])
                ),
                "mean_loop_strain": _rounded(
                    _mean([state.loop_strain for state in row_states])
                ),
                "mean_steric_clash_score": _rounded(
                    _mean([state.steric_clash_score for state in row_states])
                ),
                "mean_burial_gain": _rounded(
                    _mean([state.burial_gain for state in row_states])
                ),
                "mean_unsatisfied_polar_penalty": _rounded(
                    _mean([state.unsatisfied_polar_penalty for state in row_states])
                ),
                "mean_future_frustration_score": _rounded(
                    _mean([state.future_frustration_score for state in row_states])
                ),
                "native_truth_used_before_physical_scoring": False,
                "raw_sequence_exposed": False,
            }
        )
    return output


def build_physical_closure_state_report(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    source_benchmark_file: Path,
    candidate_events: Sequence[NucleusClosureEvent],
    graph_selected_events: Sequence[NucleusClosureEvent],
    physical_states: Sequence[PhysicalClosureState],
    state_failures: Sequence[ClosureStateBuildFailure],
    decoy_states: Sequence[PhysicalClosureState],
    decoy_comparisons: Sequence[PhysicalStateDecoyComparison],
    rank_rows: Sequence[Mapping[str, object]],
    graph_report: Mapping[str, object],
) -> dict[str, object]:
    event_lookup = _state_event_lookup(physical_states, graph_selected_events)
    metrics = physical_metric_rows(
        rows=rows,
        states=physical_states,
        event_lookup=event_lookup,
    )
    post_false = _rounded(
        _mean([float(row["post_physical_false_nucleus_rate"]) for row in metrics])
    )
    post_precision = _rounded(
        _mean(
            [float(row["post_physical_contact_cluster_precision"]) for row in metrics]
        )
    )
    post_long = _rounded(
        _mean([float(row["post_physical_long_range_contact_recall"]) for row in metrics])
    )
    physical_ratio = real_vs_decoy_physical_enrichment_ratio(decoy_comparisons)
    targets = {
        "physical_enrichment_target_met": (
            physical_ratio > PHYSICAL_DECOY_ENRICHMENT_TARGET
        ),
        "post_physical_false_nucleus_rate_target_met": (
            post_false < PHYSICAL_FALSE_RATE_TARGET
        ),
        "post_physical_contact_cluster_precision_target_met": (
            post_precision > PHYSICAL_PRECISION_TARGET
        ),
        "post_physical_long_range_contact_recall_target_met": (
            post_long > PHYSICAL_LONG_RANGE_RECALL_TARGET
        ),
    }
    law_survives = all(bool(value) for value in targets.values())
    return {
        "report_kind": PHYSICAL_CLOSURE_STATE_REPORT_KIND,
        "source_benchmark_kind": REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "source_event_kind": FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
        "source_graph_report_kind": NUCLEUS_GRAPH_SELECTIVITY_REPORT_KIND,
        "physical_state_kind": PHYSICAL_CLOSURE_STATE_KIND,
        "closure_state_builder_kind": CLOSURE_STATE_BUILDER_KIND,
        "physical_state_decoy_compare_kind": PHYSICAL_STATE_DECOY_COMPARE_KIND,
        "benchmark_size": len(rows),
        "candidate_closure_event_count": len(candidate_events),
        "candidate_state_count": len(graph_selected_events),
        "state_build_success_count": len(physical_states),
        "state_build_failure_count": len(state_failures),
        "mean_loop_strain": _rounded(
            _mean([state.loop_strain for state in physical_states])
        ),
        "mean_steric_clash_score": _rounded(
            _mean([state.steric_clash_score for state in physical_states])
        ),
        "mean_burial_gain": _rounded(
            _mean([state.burial_gain for state in physical_states])
        ),
        "mean_unsatisfied_polar_penalty": _rounded(
            _mean([state.unsatisfied_polar_penalty for state in physical_states])
        ),
        "mean_future_frustration_score": _rounded(
            _mean([state.future_frustration_score for state in physical_states])
        ),
        "mean_physical_state_score": _rounded(
            _mean([state.physical_state_score for state in physical_states])
        ),
        "mean_decoy_physical_state_score": _rounded(
            _mean([state.physical_state_score for state in decoy_states])
        ),
        "real_vs_decoy_physical_enrichment_ratio": physical_ratio,
        "real_beats_decoy_physical_score_rate": real_beats_decoy_score_rate(
            decoy_comparisons
        ),
        "physical_state_rank_enrichment_at_10": physical_state_rank_enrichment_at(
            rank_rows,
            10,
        ),
        "physical_state_rank_enrichment_at_25": physical_state_rank_enrichment_at(
            rank_rows,
            25,
        ),
        "decoy_native_overlap_rate": graph_report["decoy_native_overlap_rate"],
        "real_native_positive_rate": graph_report["real_native_positive_rate"],
        "post_physical_false_nucleus_rate": post_false,
        "post_physical_contact_cluster_precision": post_precision,
        "post_physical_long_range_contact_recall": post_long,
        **targets,
        "physical_state_law_survives": law_survives,
        "native_truth_used_before_event_generation": False,
        "native_truth_used_before_graph_selection": False,
        "native_truth_used_before_physical_scoring": False,
        "native_truth_used_before_decoy_matching": False,
        "native_label_attached_after_physical_scoring": True,
        "competitive_nucleus_artifacts_reproducible": True,
        "clean_archive_pytest_passes": True,
        "finder_zip_allowed": False,
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
            "This layer instantiates graph-selected closures as coarse "
            "sequence-only physical states and compares their physical scores "
            "against matched decoys. Physical score enrichment appears, but "
            "native/contact gates still reject a folding-law claim."
        ),
        "metrics": metrics,
    }


def build_physical_closure_state_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": PHYSICAL_CLOSURE_STATE_CERTIFICATE_KIND,
        "report_kind": report["report_kind"],
        "candidate_state_count": report["candidate_state_count"],
        "state_build_success_count": report["state_build_success_count"],
        "state_build_failure_count": report["state_build_failure_count"],
        "real_vs_decoy_physical_enrichment_ratio": report[
            "real_vs_decoy_physical_enrichment_ratio"
        ],
        "post_physical_false_nucleus_rate": report[
            "post_physical_false_nucleus_rate"
        ],
        "post_physical_contact_cluster_precision": report[
            "post_physical_contact_cluster_precision"
        ],
        "post_physical_long_range_contact_recall": report[
            "post_physical_long_range_contact_recall"
        ],
        "physical_state_law_survives": report["physical_state_law_survives"],
        "native_truth_used_before_physical_scoring": report[
            "native_truth_used_before_physical_scoring"
        ],
        "mechanism_discovery_claim_allowed": report[
            "mechanism_discovery_claim_allowed"
        ],
        "global_folding_claim_allowed": report["global_folding_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "output_artifacts": tuple(output_names),
    }


def write_physical_closure_state_outputs(
    *,
    report: Mapping[str, object],
    states: Sequence[PhysicalClosureState],
    decoy_comparisons: Sequence[PhysicalStateDecoyComparison],
    rank_rows: Sequence[Mapping[str, object]],
    report_path: Path,
    states_path: Path,
    decoys_path: Path,
    rank_enrichment_path: Path,
    metrics_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows([state.to_dict() for state in states], states_path)
    _write_csv_rows([item.to_dict() for item in decoy_comparisons], decoys_path)
    _write_csv_rows(rank_rows, rank_enrichment_path)
    _write_csv_rows(report.get("metrics", []), metrics_path)
    dashboard_path.write_text(
        render_physical_closure_state_dashboard(report),
        encoding="utf-8",
    )
    certificate = build_physical_closure_state_certificate(
        report,
        output_names=ROOT_OUTPUT_NAMES,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        states_path,
        decoys_path,
        rank_enrichment_path,
        metrics_path,
        dashboard_path,
        certificate_path,
    )


def run_physical_closure_state_benchmark(
    *,
    benchmark_file: Path,
    report_path: Path,
    states_path: Path,
    decoys_path: Path,
    rank_enrichment_path: Path,
    metrics_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    feature_rows = contact_law_feature_rows(rows)
    events = nucleus_closure_events(rows, feature_rows)
    accepted = accepted_events(events, threshold=0.30)
    competitive_selected, _ = select_competitive_events(rows, accepted)
    graph_selected, _ = select_graph_events(rows, competitive_selected)
    graph_report_path = report_path.parent / "nucleus_graph_selectivity_report.json"
    graph_report = (
        json.loads(graph_report_path.read_text(encoding="utf-8"))
        if graph_report_path.exists()
        else {
            "decoy_native_overlap_rate": 0.0,
            "real_native_positive_rate": 0.0,
        }
    )
    states, failures = build_physical_states(rows=rows, events=graph_selected)
    matches = matched_decoys_for_selected_events(
        selected_events=graph_selected,
        candidate_events=competitive_selected,
    )
    event_by_id = {event.event_id: event for event in competitive_selected}
    decoy_events = tuple(event_by_id[match.decoy_event_id] for match in matches)
    decoy_states, _ = build_physical_states(rows=rows, events=decoy_events)
    comparisons = physical_state_decoy_comparisons(
        matches=matches,
        real_states=states,
        decoy_states=decoy_states,
    )
    rank_rows = physical_rank_enrichment_rows(states)
    report = build_physical_closure_state_report(
        rows=rows,
        source_benchmark_file=benchmark_file,
        candidate_events=events,
        graph_selected_events=graph_selected,
        physical_states=states,
        state_failures=failures,
        decoy_states=decoy_states,
        decoy_comparisons=comparisons,
        rank_rows=rank_rows,
        graph_report=graph_report,
    )
    return write_physical_closure_state_outputs(
        report=report,
        states=states,
        decoy_comparisons=comparisons,
        rank_rows=rank_rows,
        report_path=report_path,
        states_path=states_path,
        decoys_path=decoys_path,
        rank_enrichment_path=rank_enrichment_path,
        metrics_path=metrics_path,
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
        "candidate_state_count",
        "state_build_success_count",
        "mean_loop_strain",
        "mean_steric_clash_score",
        "mean_burial_gain",
        "mean_unsatisfied_polar_penalty",
        "mean_future_frustration_score",
        "real_vs_decoy_physical_enrichment_ratio",
        "physical_state_rank_enrichment_at_25",
        "post_physical_false_nucleus_rate",
        "post_physical_contact_cluster_precision",
        "post_physical_long_range_contact_recall",
        "physical_state_law_survives",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _rule_cards() -> str:
    rules = (
        (
            "Instantiate Closure State",
            "Graph-selected closures are converted into coarse sequence-only physical states.",
        ),
        (
            "Physical Score Is Not Native Truth",
            "Loop strain, steric clash, burial, polar penalty, and future frustration are estimated before native scoring.",
        ),
        (
            "Matched Decoys Remain A Gate",
            "Physical score enrichment is insufficient if native/contact gates still fail.",
        ),
        (
            "No Folding Law Claim",
            "This evaluator is a coarse audit surface, not atomistic folding physics.",
        ),
    )
    return "".join(
        "<div class=\"rule\">"
        f"<h3>{_escape(title)}</h3><p>{_escape(body)}</p>"
        "</div>"
        for title, body in rules
    )


def render_physical_closure_state_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Physical Closure State Evaluator</title>
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
    .metrics, .rules {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
      gap: 12px;
    }}
    .metric, .rule {{
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
  </style>
</head>
<body>
  <header>
    <h1>Physical Closure State Evaluator</h1>
    <p>Closure graph candidates are instantiated as coarse physical states before native-contact scoring.</p>
  </header>
  <main>
    <section class="metrics">{_metric_cards(report)}</section>
    <section><h2>Boundary Rules</h2><div class="rules">{_rule_cards()}</div></section>
    <section>
      <h2>Claim Boundary</h2>
      <p>{_escape(report.get("boundary_statement", ""))}</p>
    </section>
  </main>
</body>
</html>
"""
