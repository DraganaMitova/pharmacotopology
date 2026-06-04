import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_closure_state_builder import (  # noqa: E402
    PHYSICAL_CLOSURE_STATE_CERTIFICATE_KIND,
    PHYSICAL_CLOSURE_STATE_REPORT_KIND,
    ROOT_OUTPUT_NAMES,
    build_physical_closure_state_report,
    build_physical_states,
    physical_rank_enrichment_rows,
    run_physical_closure_state_benchmark,
    write_physical_closure_state_outputs,
)
from pharmacotopology.folding_contact_law_features import (  # noqa: E402
    contact_law_feature_rows,
)
from pharmacotopology.folding_nucleus_closure_search import (  # noqa: E402
    accepted_events,
    nucleus_closure_events,
)
from pharmacotopology.folding_nucleus_competition import (  # noqa: E402
    select_competitive_events,
)
from pharmacotopology.folding_nucleus_decoy_falsification import (  # noqa: E402
    matched_decoys_for_selected_events,
)
from pharmacotopology.folding_nucleus_graph_selectivity import (  # noqa: E402
    select_graph_events,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_state_decoy_compare import (  # noqa: E402
    physical_state_decoy_comparisons,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
REL_BENCHMARK_8 = Path("data/folding_real_coordinate_visual_8.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
GRAPH_REPORT = RUN_DIR / "nucleus_graph_selectivity_report.json"
REPORT = RUN_DIR / "physical_closure_state_report.json"
STATES = RUN_DIR / "physical_closure_state_states.csv"
DECOYS = RUN_DIR / "physical_closure_state_decoys.csv"
RANK_ENRICHMENT = RUN_DIR / "physical_closure_state_rank_enrichment.csv"
METRICS = RUN_DIR / "physical_closure_state_metrics.csv"
DASHBOARD = RUN_DIR / "physical_closure_state_dashboard.html"
CERTIFICATE = RUN_DIR / "physical_closure_state_certificate.json"


_CACHE = None


def _generated():
    global _CACHE
    if _CACHE is None:
        rows = load_real_coordinate_visual_rows(REL_BENCHMARK_8)
        features = contact_law_feature_rows(rows)
        events = nucleus_closure_events(rows, features)
        accepted = accepted_events(events, threshold=0.3)
        competitive, _ = select_competitive_events(rows, accepted)
        graph_selected, _ = select_graph_events(rows, competitive)
        states, failures = build_physical_states(rows=rows, events=graph_selected)
        matches = matched_decoys_for_selected_events(
            selected_events=graph_selected,
            candidate_events=competitive,
        )
        event_by_id = {event.event_id: event for event in competitive}
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
            source_benchmark_file=REL_BENCHMARK_8,
            candidate_events=events,
            graph_selected_events=graph_selected,
            physical_states=states,
            state_failures=failures,
            decoy_states=decoy_states,
            decoy_comparisons=comparisons,
            rank_rows=rank_rows,
            graph_report=json.loads(GRAPH_REPORT.read_text(encoding="utf-8")),
        )
        _CACHE = (rows, events, graph_selected, states, failures, decoy_states, comparisons, rank_rows, report)
    return _CACHE


def test_physical_closure_state_report_tracks_physical_failure() -> None:
    _, events, graph_selected, states, failures, decoy_states, comparisons, rank_rows, report = _generated()

    assert report["report_kind"] == PHYSICAL_CLOSURE_STATE_REPORT_KIND
    assert len(events) == 5041
    assert len(graph_selected) == 320
    assert len(states) == 320
    assert len(failures) == 0
    assert len(decoy_states) == 320
    assert len(comparisons) == 320
    assert len(rank_rows) == 2
    assert report["candidate_state_count"] == 320
    assert report["state_build_success_count"] == 320
    assert report["state_build_failure_count"] == 0
    assert report["mean_loop_strain"] == 0.012782
    assert report["mean_steric_clash_score"] == 0.019839
    assert report["mean_burial_gain"] == 0.583596
    assert report["mean_unsatisfied_polar_penalty"] == 0.126116
    assert report["mean_future_frustration_score"] == 0.282609
    assert report["mean_physical_state_score"] == 0.41558
    assert report["mean_decoy_physical_state_score"] == 0.354731
    assert report["real_vs_decoy_physical_enrichment_ratio"] == 1.171535
    assert report["real_beats_decoy_physical_score_rate"] == 0.65625
    assert report["physical_state_rank_enrichment_at_10"] == 1.024
    assert report["physical_state_rank_enrichment_at_25"] == 1.1264
    assert report["post_physical_false_nucleus_rate"] == 0.609375
    assert report["post_physical_contact_cluster_precision"] == 0.043311
    assert report["post_physical_long_range_contact_recall"] == 0.372496
    assert report["physical_enrichment_target_met"] is True
    assert report["post_physical_false_nucleus_rate_target_met"] is False
    assert report["post_physical_contact_cluster_precision_target_met"] is False
    assert report["post_physical_long_range_contact_recall_target_met"] is True
    assert report["physical_state_law_survives"] is False
    assert report["native_truth_used_before_physical_scoring"] is False
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["folding_problem_solved"] is False


def test_checked_in_physical_closure_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    states = list(csv.DictReader(STATES.read_text(encoding="utf-8").splitlines()))
    decoys = list(csv.DictReader(DECOYS.read_text(encoding="utf-8").splitlines()))
    rank_rows = list(
        csv.DictReader(RANK_ENRICHMENT.read_text(encoding="utf-8").splitlines())
    )
    metrics = list(csv.DictReader(METRICS.read_text(encoding="utf-8").splitlines()))
    dashboard = DASHBOARD.read_text(encoding="utf-8")

    assert report["report_kind"] == PHYSICAL_CLOSURE_STATE_REPORT_KIND
    assert certificate["certificate_kind"] == PHYSICAL_CLOSURE_STATE_CERTIFICATE_KIND
    assert len(states) == 320
    assert len(decoys) == 320
    assert len(rank_rows) == 2
    assert len(metrics) == 8
    assert certificate["physical_state_law_survives"] is False
    assert certificate["native_truth_used_before_physical_scoring"] is False
    assert certificate["mechanism_discovery_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "Instantiate Closure State" in dashboard
    assert "Physical Score Is Not Native Truth" in dashboard
    assert "Matched Decoys Remain A Gate" in dashboard


def test_physical_closure_outputs_do_not_export_raw_sequences() -> None:
    coordinate_rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    generated_paths = list(RUN_DIR.glob("physical_closure_state_*"))

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in coordinate_rows:
            assert row.sequence not in text


def test_physical_closure_artifacts_are_reproducible(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    _, _, _, states, _, _, comparisons, rank_rows, report = _generated()
    outputs = {
        "report": tmp_path / REPORT.name,
        "states": tmp_path / STATES.name,
        "decoys": tmp_path / DECOYS.name,
        "rank": tmp_path / RANK_ENRICHMENT.name,
        "metrics": tmp_path / METRICS.name,
        "dashboard": tmp_path / DASHBOARD.name,
        "certificate": tmp_path / CERTIFICATE.name,
    }
    write_physical_closure_state_outputs(
        report=report,
        states=states,
        decoy_comparisons=comparisons,
        rank_rows=rank_rows,
        report_path=outputs["report"],
        states_path=outputs["states"],
        decoys_path=outputs["decoys"],
        rank_enrichment_path=outputs["rank"],
        metrics_path=outputs["metrics"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT,
        "states": STATES,
        "decoys": DECOYS,
        "rank": RANK_ENRICHMENT,
        "metrics": METRICS,
        "dashboard": DASHBOARD,
        "certificate": CERTIFICATE,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")


def test_physical_closure_runner_and_helpers_keep_safe_shapes(tmp_path) -> None:
    outputs = run_physical_closure_state_benchmark(
        benchmark_file=REL_BENCHMARK_8,
        report_path=tmp_path / REPORT.name,
        states_path=tmp_path / STATES.name,
        decoys_path=tmp_path / DECOYS.name,
        rank_enrichment_path=tmp_path / RANK_ENRICHMENT.name,
        metrics_path=tmp_path / METRICS.name,
        dashboard_path=tmp_path / DASHBOARD.name,
        certificate_path=tmp_path / CERTIFICATE.name,
    )
    report = json.loads(outputs[0].read_text(encoding="utf-8"))
    _, _, _, states, _, _, comparisons, rank_rows, _ = _generated()

    assert set(ROOT_OUTPUT_NAMES) == {
        "physical_closure_state_report.json",
        "physical_closure_state_states.csv",
        "physical_closure_state_decoys.csv",
        "physical_closure_state_rank_enrichment.csv",
        "physical_closure_state_metrics.csv",
        "physical_closure_state_dashboard.html",
        "physical_closure_state_certificate.json",
    }
    assert report["claim_allowed"] is False
    assert all(state.native_truth_used_before_physical_scoring is False for state in states)
    assert all(state.raw_sequence_exposed is False for state in states)
    assert all(item.native_truth_used_before_physical_scoring is False for item in comparisons)
    assert all(item.raw_sequence_exposed is False for item in comparisons)
    assert all(row["raw_sequence_exposed"] is False for row in rank_rows)

