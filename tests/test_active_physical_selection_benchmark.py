import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_physical_selection import (  # noqa: E402
    ACTIVE_PHYSICAL_SELECTION_CERTIFICATE_KIND,
    ACTIVE_PHYSICAL_SELECTION_REPORT_KIND,
    ROOT_OUTPUT_NAMES,
    build_active_physical_context,
    build_active_physical_selection_report,
    run_active_physical_selection_benchmark,
    select_events,
    selected_event_rows,
    selector_metrics,
    ablation_rows,
    write_active_physical_selection_outputs,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
REL_BENCHMARK_8 = Path("data/folding_real_coordinate_visual_8.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
REPORT = RUN_DIR / "active_physical_selection_report.json"
SELECTORS = RUN_DIR / "active_physical_selection_selectors.csv"
SELECTED_EVENTS = RUN_DIR / "active_physical_selection_selected_events.csv"
ABLATION = RUN_DIR / "active_physical_selection_ablation.csv"
DASHBOARD = RUN_DIR / "active_physical_selection_dashboard.html"
CERTIFICATE = RUN_DIR / "active_physical_selection_certificate.json"


_CACHE = None


def _generated():
    global _CACHE
    if _CACHE is None:
        rows = load_real_coordinate_visual_rows(REL_BENCHMARK_8)
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
            source_benchmark_file=REL_BENCHMARK_8,
        )
        _CACHE = (rows, context, selections, selector_rows, ablations, selected_rows, report)
    return _CACHE


def test_active_physical_selection_report_rejects_selector_law() -> None:
    _, context, selections, selector_rows, ablations, selected_rows, report = _generated()

    assert report["report_kind"] == ACTIVE_PHYSICAL_SELECTION_REPORT_KIND
    assert report["benchmark_size"] == 8
    assert report["candidate_event_count"] == 686
    assert len(context.candidate_events) == 5041
    assert len(context.competitive_events) == 686
    assert len(context.graph_events) == 320
    assert len(context.states) == 686
    assert {name: len(events) for name, events in selections.items()} == {
        "graph_only": 320,
        "physical_rerank": 320,
        "physical_gate": 74,
        "future_frustration_gate": 62,
    }
    assert len(selector_rows) == 4
    assert len(ablations) == 6
    assert len(selected_rows) == 776

    assert report["graph_only_false_nucleus_rate"] == 0.609375
    assert report["physical_rerank_false_nucleus_rate"] == 0.559375
    assert report["physical_gate_false_nucleus_rate"] == 0.316338
    assert report["future_frustration_false_nucleus_rate"] == 0.309455
    assert report["graph_only_cluster_precision"] == 0.043311
    assert report["physical_rerank_cluster_precision"] == 0.050488
    assert report["physical_gate_cluster_precision"] == 0.069038
    assert report["future_frustration_cluster_precision"] == 0.071567
    assert report["graph_only_long_range_recall"] == 0.372496
    assert report["physical_rerank_long_range_recall"] == 0.405008
    assert report["physical_gate_long_range_recall"] == 0.088436
    assert report["future_frustration_long_range_recall"] == 0.070742
    assert report["physical_rerank_real_vs_decoy_enrichment_ratio"] == 1.585882
    assert report["best_physical_term"] == "burial_gain"
    assert report["worst_physical_term"] == "future_frustration"
    assert report["physical_terms_with_positive_ablation_effect"] == ("burial_gain",)
    assert report["physical_terms_rejected_as_noise"] == (
        "loop_strain",
        "steric_clash",
        "unsatisfied_polar_penalty",
        "future_frustration",
        "decoy_margin",
    )
    assert report["active_physical_selection_survives"] is False
    assert report["native_truth_used_before_active_selection"] is False
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert report["drug_design_created"] is False
    assert report["protein_sequence_design_created"] is False


def test_checked_in_active_physical_selection_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    selectors = list(csv.DictReader(SELECTORS.read_text(encoding="utf-8").splitlines()))
    selected_events = list(
        csv.DictReader(SELECTED_EVENTS.read_text(encoding="utf-8").splitlines())
    )
    ablations = list(csv.DictReader(ABLATION.read_text(encoding="utf-8").splitlines()))
    dashboard = DASHBOARD.read_text(encoding="utf-8")

    assert report["report_kind"] == ACTIVE_PHYSICAL_SELECTION_REPORT_KIND
    assert certificate["certificate_kind"] == ACTIVE_PHYSICAL_SELECTION_CERTIFICATE_KIND
    assert len(selectors) == 4
    assert len(selected_events) == 776
    assert len(ablations) == 6
    assert certificate["active_physical_selection_survives"] is False
    assert certificate["mechanism_discovery_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "Active Physical Selection" in dashboard
    assert "Selector Comparison" in dashboard
    assert "Claim Boundary" in dashboard


def test_active_physical_selection_outputs_do_not_export_raw_sequences() -> None:
    coordinate_rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    generated_paths = list(RUN_DIR.glob("active_physical_selection_*"))

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in coordinate_rows:
            assert row.sequence not in text


def test_active_physical_selection_artifacts_are_reproducible(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    _, _, _, selector_rows, ablations, selected_rows, report = _generated()
    outputs = {
        "report": tmp_path / REPORT.name,
        "selectors": tmp_path / SELECTORS.name,
        "selected": tmp_path / SELECTED_EVENTS.name,
        "ablation": tmp_path / ABLATION.name,
        "dashboard": tmp_path / DASHBOARD.name,
        "certificate": tmp_path / CERTIFICATE.name,
    }
    write_active_physical_selection_outputs(
        report=report,
        selector_rows=selector_rows,
        selected_rows=selected_rows,
        ablations=ablations,
        report_path=outputs["report"],
        selectors_path=outputs["selectors"],
        selected_events_path=outputs["selected"],
        ablation_path=outputs["ablation"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT,
        "selectors": SELECTORS,
        "selected": SELECTED_EVENTS,
        "ablation": ABLATION,
        "dashboard": DASHBOARD,
        "certificate": CERTIFICATE,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")


def test_active_physical_selection_runner_and_helpers_keep_safe_shapes(tmp_path) -> None:
    outputs = run_active_physical_selection_benchmark(
        benchmark_file=REL_BENCHMARK_8,
        report_path=tmp_path / REPORT.name,
        selectors_path=tmp_path / SELECTORS.name,
        selected_events_path=tmp_path / SELECTED_EVENTS.name,
        ablation_path=tmp_path / ABLATION.name,
        dashboard_path=tmp_path / DASHBOARD.name,
        certificate_path=tmp_path / CERTIFICATE.name,
    )
    report = json.loads(outputs[0].read_text(encoding="utf-8"))
    _, _, _, selector_rows, ablations, selected_rows, _ = _generated()

    assert set(ROOT_OUTPUT_NAMES) == {
        "active_physical_selection_report.json",
        "active_physical_selection_selectors.csv",
        "active_physical_selection_selected_events.csv",
        "active_physical_selection_ablation.csv",
        "active_physical_selection_dashboard.html",
        "active_physical_selection_certificate.json",
    }
    assert report["claim_allowed"] is False
    assert all(row.native_truth_used_before_active_selection is False for row in selector_rows)
    assert all(row.raw_sequence_exposed is False for row in selector_rows)
    assert all(row.term_interpretation for row in ablations)
    assert all(row["native_truth_used_before_active_selection"] is False for row in selected_rows)
    assert all(row["raw_sequence_exposed"] is False for row in selected_rows)
