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
    build_folding_nucleus_closure_report,
    nucleus_closure_events,
)
from pharmacotopology.folding_nucleus_competition import (  # noqa: E402
    COMPETITIVE_NUCLEUS_SELECTION_CERTIFICATE_KIND,
    COMPETITIVE_NUCLEUS_SELECTION_REPORT_KIND,
    ROOT_OUTPUT_NAMES,
    build_competitive_nucleus_selection_report,
    select_competitive_events,
    trajectory_rows,
    write_competitive_nucleus_selection_outputs,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
REL_BENCHMARK_8 = Path("data/folding_real_coordinate_visual_8.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
REPORT = RUN_DIR / "competitive_nucleus_selection_report.json"
SELECTED = RUN_DIR / "competitive_nucleus_selection_selected_events.csv"
REJECTIONS = RUN_DIR / "competitive_nucleus_selection_rejections.csv"
COMPATIBILITY = RUN_DIR / "competitive_nucleus_selection_compatibility.csv"
TRAJECTORY = RUN_DIR / "competitive_nucleus_selection_trajectory.csv"
METRICS = RUN_DIR / "competitive_nucleus_selection_metrics.csv"
DASHBOARD = RUN_DIR / "competitive_nucleus_selection_dashboard.html"
CERTIFICATE = RUN_DIR / "competitive_nucleus_selection_certificate.json"


_CACHE = None


def _generated():
    global _CACHE
    if _CACHE is None:
        rows = load_real_coordinate_visual_rows(REL_BENCHMARK_8)
        features = contact_law_feature_rows(rows)
        events = nucleus_closure_events(rows, features)
        pre_report = build_folding_nucleus_closure_report(
            rows=rows,
            feature_rows=features,
            events=events,
            source_benchmark_file=REL_BENCHMARK_8,
        )
        selected, decisions = select_competitive_events(
            rows,
            tuple(
                event
                for event in events
                if event.nucleus_score >= 0.3 and event.frustration_cost < 0.75
            ),
        )
        report = build_competitive_nucleus_selection_report(
            rows=rows,
            events=events,
            selected_events=selected,
            decisions=decisions,
            source_benchmark_file=REL_BENCHMARK_8,
            pre_competition_report=pre_report,
        )
        trajectories = trajectory_rows(rows, selected, decisions)
        _CACHE = (rows, events, selected, decisions, trajectories, report)
    return _CACHE


def test_competitive_nucleus_selection_reports_useful_failure() -> None:
    _, events, selected, decisions, trajectories, report = _generated()

    assert report["report_kind"] == COMPETITIVE_NUCLEUS_SELECTION_REPORT_KIND
    assert len(events) == 5041
    assert report["pre_competition_event_count"] == 3846
    assert report["post_competition_selected_event_count"] == 686
    assert len(selected) == 686
    assert len(trajectories) == 686
    assert len(decisions) == 3846
    assert report["event_reduction_ratio"] == 0.821633
    assert report["pre_false_nucleus_rate"] == 0.692721
    assert report["post_false_nucleus_rate"] == 0.594592
    assert report["pre_contact_cluster_precision"] == 0.032261
    assert report["post_contact_cluster_precision"] == 0.044608
    assert report["pre_long_range_contact_recall"] == 0.913242
    assert report["post_long_range_contact_recall"] == 0.588907
    assert report["frustration_rejection_count"] == 4
    assert report["geometry_rejection_count"] == 62
    assert report["competition_rejection_count"] == 309
    assert report["overlap_rejection_count"] == 2743
    assert report["trap_rejection_count"] == 42
    assert report["selected_event_count_target_met"] is True
    assert report["long_range_contact_recall_target_met"] is True
    assert report["false_nucleus_rate_target_met"] is False
    assert report["contact_cluster_precision_target_met"] is False
    assert report["nucleus_competition_law_survives"] is False
    assert report["native_truth_used_before_selection"] is False
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["folding_problem_solved"] is False


def test_checked_in_competitive_nucleus_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    selected = list(csv.DictReader(SELECTED.read_text(encoding="utf-8").splitlines()))
    rejections = list(
        csv.DictReader(REJECTIONS.read_text(encoding="utf-8").splitlines())
    )
    compatibility = list(
        csv.DictReader(COMPATIBILITY.read_text(encoding="utf-8").splitlines())
    )
    trajectories = list(
        csv.DictReader(TRAJECTORY.read_text(encoding="utf-8").splitlines())
    )
    metrics = list(csv.DictReader(METRICS.read_text(encoding="utf-8").splitlines()))
    dashboard = DASHBOARD.read_text(encoding="utf-8")

    assert report["report_kind"] == COMPETITIVE_NUCLEUS_SELECTION_REPORT_KIND
    assert certificate["certificate_kind"] == (
        COMPETITIVE_NUCLEUS_SELECTION_CERTIFICATE_KIND
    )
    assert len(selected) == 686
    assert len(rejections) == 3160
    assert len(compatibility) == 30904
    assert len(trajectories) == 686
    assert len(metrics) == 8
    assert certificate["nucleus_competition_law_survives"] is False
    assert certificate["native_truth_used_before_selection"] is False
    assert certificate["mechanism_discovery_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "Competition Filters Event Floods" in dashboard
    assert "False Nucleus Rate Remains A Survival Gate" in dashboard
    assert "No Mechanism Discovery Claim" in dashboard
    assert {
        "compatible",
        "overlapping",
        "competing",
        "topologically_conflicting",
    } <= {row["compatibility_label"] for row in compatibility}


def test_competitive_nucleus_outputs_do_not_export_raw_sequences() -> None:
    coordinate_rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    generated_paths = list(RUN_DIR.glob("competitive_nucleus_selection_*"))

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in coordinate_rows:
            assert row.sequence not in text


def test_competitive_nucleus_artifacts_are_reproducible(
    tmp_path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(ROOT)
    _, _, selected, decisions, trajectories, report = _generated()
    outputs = {
        "report": tmp_path / REPORT.name,
        "selected": tmp_path / SELECTED.name,
        "rejections": tmp_path / REJECTIONS.name,
        "compatibility": tmp_path / COMPATIBILITY.name,
        "trajectory": tmp_path / TRAJECTORY.name,
        "metrics": tmp_path / METRICS.name,
        "dashboard": tmp_path / DASHBOARD.name,
        "certificate": tmp_path / CERTIFICATE.name,
    }
    write_competitive_nucleus_selection_outputs(
        report=report,
        selected_events=selected,
        decisions=decisions,
        trajectories=trajectories,
        report_path=outputs["report"],
        selected_events_path=outputs["selected"],
        rejections_path=outputs["rejections"],
        compatibility_path=outputs["compatibility"],
        trajectory_path=outputs["trajectory"],
        metrics_path=outputs["metrics"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT,
        "selected": SELECTED,
        "rejections": REJECTIONS,
        "compatibility": COMPATIBILITY,
        "trajectory": TRAJECTORY,
        "metrics": METRICS,
        "dashboard": DASHBOARD,
        "certificate": CERTIFICATE,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")


def test_competitive_nucleus_helpers_keep_safe_shapes() -> None:
    _, _, selected, decisions, _, report = _generated()

    assert set(ROOT_OUTPUT_NAMES) == {
        "competitive_nucleus_selection_report.json",
        "competitive_nucleus_selection_selected_events.csv",
        "competitive_nucleus_selection_rejections.csv",
        "competitive_nucleus_selection_compatibility.csv",
        "competitive_nucleus_selection_trajectory.csv",
        "competitive_nucleus_selection_metrics.csv",
        "competitive_nucleus_selection_dashboard.html",
        "competitive_nucleus_selection_certificate.json",
    }
    assert all(decision.native_truth_used_before_selection is False for decision in decisions)
    assert all(decision.raw_sequence_exposed is False for decision in decisions)
    assert all(decision.native_label_attached_after_selection is True for decision in decisions)
    assert all(event.raw_sequence_exposed is False for event in selected)
    assert report["claim_allowed"] is False

