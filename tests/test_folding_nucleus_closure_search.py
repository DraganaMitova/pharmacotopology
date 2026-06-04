import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_contact_law_features import (  # noqa: E402
    contact_law_feature_rows,
)
from pharmacotopology.folding_nucleus_closure_search import (  # noqa: E402
    FOLDING_NUCLEUS_CLOSURE_CERTIFICATE_KIND,
    FOLDING_NUCLEUS_CLOSURE_REPORT_KIND,
    ROOT_OUTPUT_NAMES,
    build_folding_nucleus_closure_report,
    event_rows,
    failure_rows_from_metrics,
    metric_rows_from_report,
    nucleus_closure_events,
    trajectory_rows,
    write_folding_nucleus_closure_outputs,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
REL_BENCHMARK_8 = Path("data/folding_real_coordinate_visual_8.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
REPORT = RUN_DIR / "folding_nucleus_closure_report.json"
EVENTS = RUN_DIR / "folding_nucleus_closure_events.csv"
TRAJECTORY = RUN_DIR / "folding_nucleus_closure_trajectory.csv"
METRICS = RUN_DIR / "folding_nucleus_closure_metrics.csv"
FAILURES = RUN_DIR / "folding_nucleus_closure_failures.csv"
DASHBOARD = RUN_DIR / "folding_nucleus_closure_dashboard.html"
CERTIFICATE = RUN_DIR / "folding_nucleus_closure_certificate.json"


_CACHE = None


def _generated():
    global _CACHE
    if _CACHE is None:
        rows = load_real_coordinate_visual_rows(REL_BENCHMARK_8)
        features = contact_law_feature_rows(rows)
        events = nucleus_closure_events(rows, features)
        report = build_folding_nucleus_closure_report(
            rows=rows,
            feature_rows=features,
            events=events,
            source_benchmark_file=REL_BENCHMARK_8,
        )
        trajectories = trajectory_rows(
            rows=rows,
            events=events,
            threshold=float(report["selected_threshold"]),
        )
        _CACHE = (rows, features, events, trajectories, report)
    return _CACHE


def test_folding_nucleus_closure_report_tracks_signal_and_failure() -> None:
    _, _, events, trajectories, report = _generated()

    assert report["report_kind"] == FOLDING_NUCLEUS_CLOSURE_REPORT_KIND
    assert report["candidate_closure_event_count"] == 5041
    assert len(events) == 5041
    assert report["nucleus_threshold_min"] == 0.3
    assert report["selected_threshold"] == 0.3
    assert report["accepted_event_count"] == 3846
    assert len(trajectories) == 3846
    assert report["native_nucleus_recall"] == 0.50826
    assert report["false_nucleus_rate"] == 0.692721
    assert report["long_range_contact_recall_after_nucleus"] == 0.913242
    assert report["pair_level_mean_long_range_contact_recall"] == 0.0
    assert report["long_range_recall_delta_vs_pair_level"] == 0.913242
    assert report["contact_cluster_precision"] == 0.032261
    assert report["closure_event_stability"] == 0.388034
    assert report["trap_event_count"] == 3056
    assert report["trajectory_native_gain"] == 0.48909
    assert report["nucleus_level_long_range_beats_pair_level"] is True
    assert report["cooperative_closure_supported"] is True
    assert report["nucleus_law_survives"] is False
    assert report["candidate_law_failure_count"] == 6
    assert report["native_truth_used_before_event_generation"] is False
    assert report["native_label_attached_after_event_generation"] is True
    assert report["row_specific_nucleus_thresholds_forbidden"] is True
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["global_folding_claim_allowed"] is False
    assert report["folding_problem_solved"] is False


def test_nucleus_event_rows_do_not_export_contact_pairs_or_sequences() -> None:
    _, _, events, _, _ = _generated()
    rows = event_rows(events)
    first = rows[0]

    assert len(rows) == 5041
    assert "sequence" not in first
    assert "raw_sequence" not in first
    assert "candidate_region_pairs" not in first
    assert first["native_truth_used_before_event_generation"] is False
    assert first["raw_sequence_exposed"] is False


def test_checked_in_folding_nucleus_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    events = list(csv.DictReader(EVENTS.read_text(encoding="utf-8").splitlines()))
    trajectories = list(
        csv.DictReader(TRAJECTORY.read_text(encoding="utf-8").splitlines())
    )
    metrics = list(csv.DictReader(METRICS.read_text(encoding="utf-8").splitlines()))
    failures = list(csv.DictReader(FAILURES.read_text(encoding="utf-8").splitlines()))
    dashboard = DASHBOARD.read_text(encoding="utf-8")

    assert report["report_kind"] == FOLDING_NUCLEUS_CLOSURE_REPORT_KIND
    assert certificate["certificate_kind"] == FOLDING_NUCLEUS_CLOSURE_CERTIFICATE_KIND
    assert len(events) == 5041
    assert len(trajectories) == 3846
    assert len(metrics) == 8
    assert len(failures) == 6
    assert certificate["nucleus_level_long_range_beats_pair_level"] is True
    assert certificate["nucleus_law_survives"] is False
    assert certificate["native_truth_used_before_event_generation"] is False
    assert certificate["row_specific_nucleus_thresholds_forbidden"] is True
    assert certificate["mechanism_discovery_claim_allowed"] is False
    assert certificate["global_folding_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "Closure Events, Not Pair Events" in dashboard
    assert "Native Labels After Event Generation" in dashboard
    assert "Long-Range Recovery Is The First Target" in dashboard
    assert "No Folding Law Claim" in dashboard


def test_folding_nucleus_outputs_do_not_export_raw_sequences() -> None:
    coordinate_rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    generated_paths = list(RUN_DIR.glob("folding_nucleus_closure_*"))

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in coordinate_rows:
            assert row.sequence not in text


def test_folding_nucleus_closure_artifacts_are_reproducible(
    tmp_path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(ROOT)
    _, _, events, trajectories, report = _generated()
    outputs = {
        "report": tmp_path / REPORT.name,
        "events": tmp_path / EVENTS.name,
        "trajectory": tmp_path / TRAJECTORY.name,
        "metrics": tmp_path / METRICS.name,
        "failures": tmp_path / FAILURES.name,
        "dashboard": tmp_path / DASHBOARD.name,
        "certificate": tmp_path / CERTIFICATE.name,
    }
    write_folding_nucleus_closure_outputs(
        report=report,
        events=events,
        trajectories=trajectories,
        report_path=outputs["report"],
        events_path=outputs["events"],
        trajectory_path=outputs["trajectory"],
        metrics_path=outputs["metrics"],
        failures_path=outputs["failures"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT,
        "events": EVENTS,
        "trajectory": TRAJECTORY,
        "metrics": METRICS,
        "failures": FAILURES,
        "dashboard": DASHBOARD,
        "certificate": CERTIFICATE,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")


def test_folding_nucleus_helpers_keep_safe_shapes() -> None:
    _, _, events, _, report = _generated()
    metric_rows = metric_rows_from_report(report)
    failures = failure_rows_from_metrics(metric_rows)

    assert len(metric_rows) == 8
    assert len(failures) == 6
    assert set(ROOT_OUTPUT_NAMES) == {
        "folding_nucleus_closure_report.json",
        "folding_nucleus_closure_events.csv",
        "folding_nucleus_closure_trajectory.csv",
        "folding_nucleus_closure_metrics.csv",
        "folding_nucleus_closure_failures.csv",
        "folding_nucleus_closure_dashboard.html",
        "folding_nucleus_closure_certificate.json",
    }
    assert all(event.raw_sequence_exposed is False for event in events)
    assert all(
        event.native_truth_used_before_event_generation is False
        for event in events
    )
