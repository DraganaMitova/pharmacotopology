import csv
from dataclasses import replace
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.artifact_io import stable_float_text  # noqa: E402
from pharmacotopology.folding_coupling_nucleus_selector import (  # noqa: E402
    COUPLING_NUCLEUS_SELECTOR_CERTIFICATE_KIND,
    COUPLING_NUCLEUS_SELECTOR_REPORT_KIND,
    EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_KIND,
    ROOT_OUTPUT_NAMES,
    build_coupling_nucleus_context,
    build_coupling_nucleus_selector_report,
    decoy_rows,
    run_coupling_nucleus_selector_benchmark,
    select_coupling_events,
    selected_event_rows,
    selector_metrics,
    write_coupling_nucleus_selector_outputs,
)
from pharmacotopology.folding_evolutionary_constraints import (  # noqa: E402
    EVOLUTIONARY_COUPLING_LAYER_KIND,
    load_coupling_dataset,
    validate_coupling_dataset,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
REL_BENCHMARK_8 = Path("data/folding_real_coordinate_visual_8.locked.json")
COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_couplings.locked.json"
REL_COUPLINGS = Path("data/folding_real_coordinate_visual_8_couplings.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
REPORT = RUN_DIR / "coupling_nucleus_selector_report.json"
SELECTORS = RUN_DIR / "coupling_nucleus_selector_selectors.csv"
SELECTED_EVENTS = RUN_DIR / "coupling_nucleus_selector_selected_events.csv"
ASSESSMENTS = RUN_DIR / "coupling_nucleus_selector_assessments.csv"
DECOYS = RUN_DIR / "coupling_nucleus_selector_decoys.csv"
DASHBOARD = RUN_DIR / "coupling_nucleus_selector_dashboard.html"
CERTIFICATE = RUN_DIR / "coupling_nucleus_selector_certificate.json"


_CACHE = None


def _generated():
    global _CACHE
    if _CACHE is None:
        rows = load_real_coordinate_visual_rows(REL_BENCHMARK_8)
        couplings = load_coupling_dataset(REL_COUPLINGS)
        context = build_coupling_nucleus_context(
            rows=rows,
            coupling_dataset=couplings,
        )
        selections = {
            "graph_only": select_coupling_events(
                context,
                selector_name="graph_only",
            ),
            "physical_rerank": select_coupling_events(
                context,
                selector_name="physical_rerank",
            ),
            "coupling_rerank": select_coupling_events(
                context,
                selector_name="coupling_rerank",
            ),
            "coupling_trace_loop": select_coupling_events(
                context,
                selector_name="coupling_trace_loop",
            ),
            "coupling_future_gate": select_coupling_events(
                context,
                selector_name="coupling_future_gate",
            ),
            "coupling_decoy_falsifier": select_coupling_events(
                context,
                selector_name="coupling_decoy_falsifier",
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
        selected_rows = selected_event_rows(context, selections)
        decoys = decoy_rows(context, selections["coupling_rerank"])
        report = build_coupling_nucleus_selector_report(
            context=context,
            selector_rows=selector_rows,
            source_benchmark_file=REL_BENCHMARK_8,
            coupling_file=REL_COUPLINGS,
        )
        _CACHE = (rows, couplings, context, selections, selector_rows, selected_rows, decoys, report)
    return _CACHE


def test_coupling_nucleus_selector_finds_oracle_control_boundary() -> None:
    _, couplings, context, selections, selector_rows, selected_rows, decoys, report = _generated()

    assert couplings.layer_kind == EVOLUTIONARY_COUPLING_LAYER_KIND
    assert len(couplings.constraints) == 745
    assert len(context.competitive_events) == 686
    assert len(context.assessments) == 686
    assert {name: len(events) for name, events in selections.items()} == {
        "graph_only": 320,
        "physical_rerank": 320,
        "coupling_rerank": 320,
        "coupling_trace_loop": 54,
        "coupling_future_gate": 72,
        "coupling_decoy_falsifier": 48,
    }
    assert len(selector_rows) == 6
    assert len(selected_rows) == 1134
    assert len(decoys) == 320

    assert report["report_kind"] == COUPLING_NUCLEUS_SELECTOR_REPORT_KIND
    assert report["benchmark_size"] == 8
    assert report["candidate_event_count"] == 686
    assert report["coupling_constraint_count"] == 745
    assert report["graph_only_false_nucleus_rate"] == 0.609375
    assert report["physical_rerank_false_nucleus_rate"] == 0.559375
    assert report["coupling_rerank_false_nucleus_rate"] == 0.34375
    assert report["coupling_trace_loop_false_nucleus_rate"] == 0.0
    assert report["coupling_future_gate_false_nucleus_rate"] == 0.0
    assert report["coupling_decoy_falsifier_false_nucleus_rate"] == 0.0
    assert report["coupling_rerank_cluster_precision"] == 0.072461
    assert report["coupling_trace_loop_cluster_precision"] == 0.164128
    assert report["coupling_rerank_long_range_recall"] == 0.565164
    assert report["coupling_trace_loop_long_range_recall"] == 0.397896
    assert report["coupling_rerank_constraint_recall"] == 0.589824
    assert report["coupling_trace_loop_constraint_recall"] == 0.466536
    assert report["coupling_rerank_real_vs_decoy_enrichment_ratio"] == 1.931956
    assert report["coupling_trace_loop_real_vs_decoy_enrichment_ratio"] == 1.693887
    assert report["coupling_selector_targets_met"] is True
    assert (
        report["external_batch_kind"]
        == EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_KIND
    )
    assert report["external_evolutionary_couplings_used"] is False
    assert report["per_constraint_coordinate_truth_used"] is True
    assert report["coordinate_truth_used_to_build_constraints"] is True
    assert report["native_truth_used_before_coupling_selection"] is True
    assert report["oracle_constraint_control"] is True
    assert report["claim_mode_validation_passed"] is False
    assert set(report["claim_mode_validation_failures"]) == {
        "external_evolutionary_couplings_used=false",
        "coupling_source_kind=coordinate_oracle_surrogate_for_missing_evolutionary_channel_v1",
        "coordinate_truth_used_to_build_constraints=true",
        "native_truth_used_before_coupling_selection=true",
        "oracle_constraint_control=true",
    }
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert report["claim_allowed"] is False


def test_checked_in_coupling_nucleus_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    selectors = list(csv.DictReader(SELECTORS.read_text(encoding="utf-8").splitlines()))
    selected_events = list(
        csv.DictReader(SELECTED_EVENTS.read_text(encoding="utf-8").splitlines())
    )
    assessments = list(
        csv.DictReader(ASSESSMENTS.read_text(encoding="utf-8").splitlines())
    )
    decoys = list(csv.DictReader(DECOYS.read_text(encoding="utf-8").splitlines()))
    dashboard = DASHBOARD.read_text(encoding="utf-8")

    assert report["report_kind"] == COUPLING_NUCLEUS_SELECTOR_REPORT_KIND
    assert certificate["certificate_kind"] == (
        COUPLING_NUCLEUS_SELECTOR_CERTIFICATE_KIND
    )
    assert len(selectors) == 6
    assert len(selected_events) == 1134
    assert len(assessments) == 686
    assert len(decoys) == 320
    assert certificate["coupling_selector_targets_met"] is True
    assert certificate["oracle_constraint_control"] is True
    assert certificate["per_constraint_coordinate_truth_used"] is True
    assert certificate["claim_mode_validation_passed"] is False
    assert certificate["external_evolutionary_couplings_used"] is False
    assert certificate["mechanism_discovery_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "Coupling-Preserved Nucleus Selector" in dashboard
    assert "Selector Comparison" in dashboard
    assert "Claim Boundary" in dashboard


def test_coupling_nucleus_outputs_do_not_export_raw_sequences() -> None:
    coordinate_rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    generated_paths = list(RUN_DIR.glob("coupling_nucleus_selector_*")) + [COUPLINGS]

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in coordinate_rows:
            assert row.sequence not in text


def test_coupling_nucleus_artifacts_are_reproducible(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    _, _, context, _, selector_rows, selected_rows, decoys, report = _generated()
    outputs = {
        "report": tmp_path / REPORT.name,
        "selectors": tmp_path / SELECTORS.name,
        "selected": tmp_path / SELECTED_EVENTS.name,
        "assessments": tmp_path / ASSESSMENTS.name,
        "decoys": tmp_path / DECOYS.name,
        "dashboard": tmp_path / DASHBOARD.name,
        "certificate": tmp_path / CERTIFICATE.name,
    }
    write_coupling_nucleus_selector_outputs(
        report=report,
        selector_rows=selector_rows,
        selected_rows=selected_rows,
        assessments=context.assessments,
        decoys=decoys,
        report_path=outputs["report"],
        selectors_path=outputs["selectors"],
        selected_events_path=outputs["selected"],
        assessments_path=outputs["assessments"],
        decoys_path=outputs["decoys"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT,
        "selectors": SELECTORS,
        "selected": SELECTED_EVENTS,
        "assessments": ASSESSMENTS,
        "decoys": DECOYS,
        "dashboard": DASHBOARD,
        "certificate": CERTIFICATE,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")


def test_coupling_nucleus_runner_and_helpers_keep_safe_shapes(tmp_path) -> None:
    outputs = run_coupling_nucleus_selector_benchmark(
        benchmark_file=REL_BENCHMARK_8,
        coupling_file=REL_COUPLINGS,
        report_path=tmp_path / REPORT.name,
        selectors_path=tmp_path / SELECTORS.name,
        selected_events_path=tmp_path / SELECTED_EVENTS.name,
        assessments_path=tmp_path / ASSESSMENTS.name,
        decoys_path=tmp_path / DECOYS.name,
        dashboard_path=tmp_path / DASHBOARD.name,
        certificate_path=tmp_path / CERTIFICATE.name,
    )
    report = json.loads(outputs[0].read_text(encoding="utf-8"))
    _, _, _, _, selector_rows, selected_rows, decoys, _ = _generated()

    assert set(ROOT_OUTPUT_NAMES) == {
        "coupling_nucleus_selector_report.json",
        "coupling_nucleus_selector_selectors.csv",
        "coupling_nucleus_selector_selected_events.csv",
        "coupling_nucleus_selector_assessments.csv",
        "coupling_nucleus_selector_decoys.csv",
        "coupling_nucleus_selector_dashboard.html",
        "coupling_nucleus_selector_certificate.json",
    }
    assert report["coupling_selector_targets_met"] is True
    assert report["claim_allowed"] is False
    assert all(row.raw_sequence_exposed is False for row in selector_rows)
    assert all(
        row.native_truth_used_before_coupling_selection is True
        for row in selector_rows
    )
    assert all(row["raw_sequence_exposed"] is False for row in selected_rows)
    assert all(row["raw_sequence_exposed"] is False for row in decoys)


def test_stable_csv_float_text_keeps_locked_artifact_style() -> None:
    assert stable_float_text(1.0) == "1.0"
    assert stable_float_text(0.6200000000000001) == "0.62"
    assert stable_float_text(0.32633599999999996) == "0.326336"
    assert stable_float_text(-0.0000001) == "0.0"


def test_per_constraint_coordinate_truth_taints_claim_mode() -> None:
    _, couplings, context, selections, _, _, _, _ = _generated()
    external_labeled_oracle = replace(
        couplings,
        coordinate_truth_used_to_build_constraints=False,
        native_truth_used_before_coupling_selection=False,
        external_evolutionary_couplings_used=True,
        coupling_source_kind="external_msa_dca_couplings_v1",
    )
    tainted_context = replace(context, coupling_dataset=external_labeled_oracle)
    selector_rows = tuple(
        selector_metrics(
            tainted_context,
            selector_name=name,
            selected_events=events,
        )
        for name, events in selections.items()
    )
    report = build_coupling_nucleus_selector_report(
        context=tainted_context,
        selector_rows=selector_rows,
        source_benchmark_file=REL_BENCHMARK_8,
        coupling_file=REL_COUPLINGS,
    )

    assert external_labeled_oracle.per_constraint_coordinate_truth_used is True
    assert external_labeled_oracle.coordinate_truth_tainted is True
    assert report["external_evolutionary_couplings_used"] is True
    assert report["per_constraint_coordinate_truth_used"] is True
    assert report["coordinate_truth_used_to_build_constraints"] is True
    assert report["oracle_constraint_control"] is True
    assert report["claim_mode_validation_passed"] is False
    assert "coordinate_truth_used_to_build_constraints=true" in report[
        "claim_mode_validation_failures"
    ]
    assert report["claim_allowed"] is False
    assert all(
        row.coordinate_truth_used_to_build_constraints is True
        for row in selector_rows
    )


def test_constraint_raw_sequence_exposure_is_rejected() -> None:
    rows, couplings, *_ = _generated()
    unsafe_constraint = replace(couplings.constraints[0], raw_sequence_exposed=True)
    unsafe_dataset = replace(couplings, constraints=(unsafe_constraint,))

    with pytest.raises(ValueError, match="must not expose raw sequence text"):
        validate_coupling_dataset(rows=rows, dataset=unsafe_dataset)
