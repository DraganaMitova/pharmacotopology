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
    nucleus_closure_events,
)
from pharmacotopology.folding_nucleus_competition import (  # noqa: E402
    select_competitive_events,
)
from pharmacotopology.folding_nucleus_graph_selectivity import (  # noqa: E402
    NUCLEUS_GRAPH_SELECTIVITY_CERTIFICATE_KIND,
    NUCLEUS_GRAPH_SELECTIVITY_REPORT_KIND,
    ROOT_OUTPUT_NAMES,
    build_nucleus_graph_selectivity_report,
    graph_rows,
    metric_rows_from_report,
    select_graph_events,
    write_nucleus_graph_selectivity_outputs,
)
from pharmacotopology.folding_nucleus_decoy_falsification import (  # noqa: E402
    matched_decoys_for_selected_events,
)
from pharmacotopology.folding_nucleus_rank_enrichment import (  # noqa: E402
    RANK_ENRICHMENT_CUTOFFS,
    rank_enrichment_rows_for_row,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
REL_BENCHMARK_8 = Path("data/folding_real_coordinate_visual_8.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
COMPETITIVE_REPORT = RUN_DIR / "competitive_nucleus_selection_report.json"
REPORT = RUN_DIR / "nucleus_graph_selectivity_report.json"
GRAPHS = RUN_DIR / "nucleus_graph_selectivity_graphs.csv"
SELECTED = RUN_DIR / "nucleus_graph_selectivity_selected_events.csv"
REJECTIONS = RUN_DIR / "nucleus_graph_selectivity_rejections.csv"
DECOYS = RUN_DIR / "nucleus_graph_selectivity_decoys.csv"
RANK_ENRICHMENT = RUN_DIR / "nucleus_graph_selectivity_rank_enrichment.csv"
METRICS = RUN_DIR / "nucleus_graph_selectivity_metrics.csv"
DASHBOARD = RUN_DIR / "nucleus_graph_selectivity_dashboard.html"
CERTIFICATE = RUN_DIR / "nucleus_graph_selectivity_certificate.json"


_CACHE = None


def _generated():
    global _CACHE
    if _CACHE is None:
        rows = load_real_coordinate_visual_rows(REL_BENCHMARK_8)
        features = contact_law_feature_rows(rows)
        events = nucleus_closure_events(rows, features)
        accepted = tuple(
            event
            for event in events
            if event.nucleus_score >= 0.3 and event.frustration_cost < 0.75
        )
        competitive_selected, _ = select_competitive_events(rows, accepted)
        selected, decisions = select_graph_events(rows, competitive_selected)
        score_by_event = {
            decision.event_id: decision.graph_core_score for decision in decisions
        }
        competitive_by_row = {}
        for event in competitive_selected:
            competitive_by_row.setdefault(event.row_id, []).append(event)
        rank_rows = []
        for row in rows:
            rank_rows.extend(
                rank_enrichment_rows_for_row(
                    row_id=row.row_id,
                    source_accession=row.source_accession,
                    candidate_events=tuple(
                        competitive_by_row.get(row.row_id, ())
                    ),
                    score_function=lambda event, lookup=score_by_event: lookup[
                        event.event_id
                    ],
                    cutoffs=RANK_ENRICHMENT_CUTOFFS,
                )
            )
        decoys = matched_decoys_for_selected_events(
            selected_events=selected,
            candidate_events=competitive_selected,
        )
        graphs = graph_rows(rows, selected, decisions)
        report = build_nucleus_graph_selectivity_report(
            rows=rows,
            events=events,
            pre_graph_events=competitive_selected,
            selected_events=selected,
            decisions=decisions,
            decoys=decoys,
            rank_rows=tuple(rank_rows),
            source_benchmark_file=REL_BENCHMARK_8,
            pre_graph_report=json.loads(COMPETITIVE_REPORT.read_text(encoding="utf-8")),
        )
        _CACHE = (rows, events, competitive_selected, selected, decisions, decoys, rank_rows, graphs, report)
    return _CACHE


def test_nucleus_graph_selectivity_reports_decoy_falsification() -> None:
    _, events, competitive, selected, decisions, decoys, rank_rows, graphs, report = _generated()

    assert report["report_kind"] == NUCLEUS_GRAPH_SELECTIVITY_REPORT_KIND
    assert len(events) == 5041
    assert len(competitive) == 686
    assert len(selected) == 320
    assert len(decisions) == 686
    assert len(decoys) == 320
    assert len(rank_rows) == 24
    assert len(graphs) == 8
    assert report["pre_graph_selected_event_count"] == 686
    assert report["post_graph_selected_event_count"] == 320
    assert report["pre_false_nucleus_rate"] == 0.594592
    assert report["post_false_nucleus_rate"] == 0.609375
    assert report["pre_contact_cluster_precision"] == 0.044608
    assert report["post_contact_cluster_precision"] == 0.043311
    assert report["pre_long_range_contact_recall"] == 0.588907
    assert report["post_long_range_contact_recall"] == 0.372496
    assert report["rank_enrichment_at_10"] == 0.97687
    assert report["rank_enrichment_at_25"] == 0.959248
    assert report["rank_enrichment_at_50"] == 0.946458
    assert report["native_positive_top_rank_rate"] == 0.395
    assert report["decoy_native_overlap_rate"] == 0.484375
    assert report["real_native_positive_rate"] == 0.390625
    assert report["real_vs_decoy_enrichment_ratio"] == 0.806452
    assert report["isolated_event_rejection_count"] == 1
    assert report["hydrophobic_only_rejection_count"] == 19
    assert report["unsupported_long_span_rejection_count"] == 0
    assert report["topology_conflict_rejection_count"] == 10
    assert report["trap_graph_rejection_count"] == 0
    assert report["post_graph_selected_event_target_met"] is True
    assert report["post_long_range_contact_recall_target_met"] is True
    assert report["post_false_nucleus_rate_target_met"] is False
    assert report["post_contact_cluster_precision_target_met"] is False
    assert report["decoy_enrichment_target_met"] is False
    assert report["nucleus_graph_law_survives"] is False
    assert report["competitive_nucleus_artifacts_reproducible"] is True
    assert report["clean_archive_pytest_passes"] is True
    assert report["finder_zip_allowed"] is False
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["folding_problem_solved"] is False


def test_checked_in_nucleus_graph_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    graphs = list(csv.DictReader(GRAPHS.read_text(encoding="utf-8").splitlines()))
    selected = list(csv.DictReader(SELECTED.read_text(encoding="utf-8").splitlines()))
    rejections = list(
        csv.DictReader(REJECTIONS.read_text(encoding="utf-8").splitlines())
    )
    decoys = list(csv.DictReader(DECOYS.read_text(encoding="utf-8").splitlines()))
    rank_rows = list(
        csv.DictReader(RANK_ENRICHMENT.read_text(encoding="utf-8").splitlines())
    )
    metrics = list(csv.DictReader(METRICS.read_text(encoding="utf-8").splitlines()))
    dashboard = DASHBOARD.read_text(encoding="utf-8")

    assert report["report_kind"] == NUCLEUS_GRAPH_SELECTIVITY_REPORT_KIND
    assert certificate["certificate_kind"] == (
        NUCLEUS_GRAPH_SELECTIVITY_CERTIFICATE_KIND
    )
    assert len(graphs) == 8
    assert len(selected) == 320
    assert len(rejections) == 366
    assert len(decoys) == 320
    assert len(rank_rows) == 24
    assert len(metrics) == 8
    assert certificate["nucleus_graph_law_survives"] is False
    assert certificate["finder_zip_allowed"] is False
    assert certificate["native_truth_used_before_graph_selection"] is False
    assert certificate["native_truth_used_before_decoy_matching"] is False
    assert certificate["mechanism_discovery_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "Score Graph Cores, Not Events Alone" in dashboard
    assert "Matched Decoys Must Be Beaten" in dashboard
    assert "False Nuclei Remain A Gate" in dashboard


def test_nucleus_graph_outputs_do_not_export_raw_sequences() -> None:
    coordinate_rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    generated_paths = list(RUN_DIR.glob("nucleus_graph_selectivity_*"))

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in coordinate_rows:
            assert row.sequence not in text


def test_nucleus_graph_selectivity_artifacts_are_reproducible(
    tmp_path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(ROOT)
    _, _, _, selected, decisions, decoys, rank_rows, graphs, report = _generated()
    outputs = {
        "report": tmp_path / REPORT.name,
        "graphs": tmp_path / GRAPHS.name,
        "selected": tmp_path / SELECTED.name,
        "rejections": tmp_path / REJECTIONS.name,
        "decoys": tmp_path / DECOYS.name,
        "rank_enrichment": tmp_path / RANK_ENRICHMENT.name,
        "metrics": tmp_path / METRICS.name,
        "dashboard": tmp_path / DASHBOARD.name,
        "certificate": tmp_path / CERTIFICATE.name,
    }
    write_nucleus_graph_selectivity_outputs(
        report=report,
        graphs=graphs,
        selected_events=selected,
        decisions=decisions,
        decoys=decoys,
        rank_rows=tuple(rank_rows),
        report_path=outputs["report"],
        graphs_path=outputs["graphs"],
        selected_events_path=outputs["selected"],
        rejections_path=outputs["rejections"],
        decoys_path=outputs["decoys"],
        rank_enrichment_path=outputs["rank_enrichment"],
        metrics_path=outputs["metrics"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT,
        "graphs": GRAPHS,
        "selected": SELECTED,
        "rejections": REJECTIONS,
        "decoys": DECOYS,
        "rank_enrichment": RANK_ENRICHMENT,
        "metrics": METRICS,
        "dashboard": DASHBOARD,
        "certificate": CERTIFICATE,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")


def test_nucleus_graph_helpers_keep_safe_shapes() -> None:
    _, _, _, selected, decisions, decoys, rank_rows, _, report = _generated()

    assert set(ROOT_OUTPUT_NAMES) == {
        "nucleus_graph_selectivity_report.json",
        "nucleus_graph_selectivity_graphs.csv",
        "nucleus_graph_selectivity_selected_events.csv",
        "nucleus_graph_selectivity_rejections.csv",
        "nucleus_graph_selectivity_decoys.csv",
        "nucleus_graph_selectivity_rank_enrichment.csv",
        "nucleus_graph_selectivity_metrics.csv",
        "nucleus_graph_selectivity_dashboard.html",
        "nucleus_graph_selectivity_certificate.json",
    }
    assert len(metric_rows_from_report(report)) == 8
    assert all(decision.native_truth_used_before_graph_selection is False for decision in decisions)
    assert all(decision.raw_sequence_exposed is False for decision in decisions)
    assert all(decoy.native_truth_used_before_decoy_matching is False for decoy in decoys)
    assert all(decoy.raw_sequence_exposed is False for decoy in decoys)
    assert all(row.native_truth_used_before_ranking is False for row in rank_rows)
    assert all(row.raw_sequence_exposed is False for row in rank_rows)
    assert all(event.raw_sequence_exposed is False for event in selected)
    assert report["claim_allowed"] is False

