import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_contact_law_features import (  # noqa: E402
    CONTACT_LAW_FEATURE_KIND,
    contact_law_feature_rows,
    feature_rows_by_row_id,
)
from pharmacotopology.folding_contact_threshold_search import (  # noqa: E402
    CONTACT_LAW_THRESHOLD_CERTIFICATE_KIND,
    CONTACT_LAW_THRESHOLD_REPORT_KIND,
    ROOT_OUTPUT_NAMES,
    build_contact_law_threshold_report,
    failure_rows,
    leave_one_out_rows,
    threshold_grid_rows,
    write_contact_law_threshold_outputs,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
REL_BENCHMARK_8 = Path("data/folding_real_coordinate_visual_8.locked.json")
RUN_DIR = ROOT / "first_contact_clean_pharmacotopology_layer_run"
REPORT = RUN_DIR / "contact_law_threshold_report.json"
ROWS = RUN_DIR / "contact_law_threshold_rows.csv"
GRID = RUN_DIR / "contact_law_threshold_grid.csv"
HOLDOUT = RUN_DIR / "contact_law_threshold_holdout.csv"
FAILURES = RUN_DIR / "contact_law_threshold_failures.csv"
DASHBOARD = RUN_DIR / "contact_law_threshold_dashboard.html"
CERTIFICATE = RUN_DIR / "contact_law_threshold_certificate.json"


_CACHE = None


def _generated():
    global _CACHE
    if _CACHE is None:
        coordinate_rows = load_real_coordinate_visual_rows(REL_BENCHMARK_8)
        feature_rows = contact_law_feature_rows(coordinate_rows)
        row_groups = feature_rows_by_row_id(feature_rows)
        grid = threshold_grid_rows(row_groups)
        holdout = leave_one_out_rows(row_groups)
        report = build_contact_law_threshold_report(
            feature_rows=feature_rows,
            source_benchmark_file=REL_BENCHMARK_8,
            grid_rows=grid,
            holdout_rows=holdout,
        )
        failures = failure_rows(
            holdout_rows=holdout,
            best_model_id=str(report["best_law_candidate_model"]),
        )
        _CACHE = (coordinate_rows, feature_rows, grid, holdout, failures, report)
    return _CACHE


def test_contact_law_threshold_report_rejects_current_scalar_law() -> None:
    _, feature_rows, grid, holdout, failures, report = _generated()

    assert report["report_kind"] == CONTACT_LAW_THRESHOLD_REPORT_KIND
    assert report["feature_kind"] == CONTACT_LAW_FEATURE_KIND
    assert report["pair_feature_row_count"] == 94546
    assert len(feature_rows) == 94546
    assert report["native_pair_label_count"] == 2988
    assert report["threshold_grid_row_count"] == 505
    assert len(grid) == 505
    assert report["holdout_row_count"] == 32
    assert len(holdout) == 32
    assert report["law_search_completed"] is True
    assert report["artifact_reproducible"] is True
    assert report["uploaded_zip_pytest_passes"] is True
    assert report["native_truth_used_before_feature_generation"] is False
    assert report["row_specific_thresholds_forbidden"] is True

    assert report["current_scalar_score_best_global_threshold"] == 0.52
    assert report["current_scalar_score_best_global_f1"] == 0.150906
    assert report["current_scalar_score_best_global_micro_f1"] == 0.147772
    assert report["current_scalar_score_threshold_std"] == 0.217715
    assert report["current_scalar_score_threshold_stable"] is False
    assert report["current_scalar_score_law_rejected"] is True

    assert report["pair_only_best_f1"] == 0.178009
    assert report["best_law_candidate_model"] == "pair_plus_entropy_score"
    assert report["best_law_candidate_loo_mean_test_f1"] == 0.233569
    assert report["best_law_candidate_loo_threshold_std"] == 0.007071
    assert report["best_law_candidate_survives"] is False
    assert report["law_generalizes"] is False
    assert len(failures) == 4
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["global_folding_claim_allowed"] is False
    assert report["folding_problem_solved"] is False


def test_contact_law_feature_rows_have_required_sequence_only_fields() -> None:
    _, feature_rows, _, _, _, _ = _generated()
    first = feature_rows[0].to_dict()

    for field in (
        "sequence_separation",
        "normalized_separation",
        "local_i_to_i4_support",
        "helix_window_support",
        "beta_window_support",
        "hydrophobic_pair_support",
        "aromatic_anchor_support",
        "opposite_charge_support",
        "same_charge_penalty",
        "breaker_penalty",
        "loop_entropy_cost",
        "cluster_neighbor_support",
        "parallel_contact_support",
        "isolation_penalty",
        "native_contact",
        "pair_plus_cluster_plus_entropy_score",
    ):
        assert field in first
    assert "sequence" not in first
    assert "raw_sequence" not in first
    assert sum(1 for row in feature_rows if row.native_contact) == 2988


def test_checked_in_contact_law_threshold_outputs_have_expected_surfaces() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(ROWS.read_text(encoding="utf-8").splitlines()))
    grid = list(csv.DictReader(GRID.read_text(encoding="utf-8").splitlines()))
    holdout = list(csv.DictReader(HOLDOUT.read_text(encoding="utf-8").splitlines()))
    failures = list(csv.DictReader(FAILURES.read_text(encoding="utf-8").splitlines()))
    dashboard = DASHBOARD.read_text(encoding="utf-8")

    assert report["report_kind"] == CONTACT_LAW_THRESHOLD_REPORT_KIND
    assert certificate["certificate_kind"] == CONTACT_LAW_THRESHOLD_CERTIFICATE_KIND
    assert len(rows) == 94546
    assert len(grid) == 505
    assert len(holdout) == 32
    assert len(failures) == 4
    assert "sequence" not in rows[0]
    assert "raw_sequence" not in rows[0]
    assert certificate["current_scalar_score_law_rejected"] is True
    assert certificate["law_generalizes"] is False
    assert certificate["native_truth_used_before_feature_generation"] is False
    assert certificate["row_specific_thresholds_forbidden"] is True
    assert certificate["mechanism_discovery_claim_allowed"] is False
    assert certificate["global_folding_claim_allowed"] is False
    assert certificate["folding_problem_solved"] is False
    assert "Native Labels After Features" in dashboard
    assert "No Row-Specific Thresholds" in dashboard
    assert "Scalar Law Must Survive Falsification" in dashboard
    assert "No Discovery Claim" in dashboard


def test_contact_law_threshold_outputs_do_not_export_raw_sequences() -> None:
    coordinate_rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    generated_paths = list(RUN_DIR.glob("contact_law_threshold_*"))

    for path in generated_paths:
        text = path.read_text(encoding="utf-8")
        for row in coordinate_rows:
            assert row.sequence not in text


def test_contact_law_threshold_artifacts_are_reproducible(
    tmp_path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(ROOT)
    _, feature_rows, grid, holdout, failures, report = _generated()
    outputs = {
        "report": tmp_path / REPORT.name,
        "rows": tmp_path / ROWS.name,
        "grid": tmp_path / GRID.name,
        "holdout": tmp_path / HOLDOUT.name,
        "failures": tmp_path / FAILURES.name,
        "dashboard": tmp_path / DASHBOARD.name,
        "certificate": tmp_path / CERTIFICATE.name,
    }
    write_contact_law_threshold_outputs(
        report=report,
        feature_rows=feature_rows,
        grid_rows=grid,
        holdout_rows=holdout,
        failures=failures,
        report_path=outputs["report"],
        rows_path=outputs["rows"],
        grid_path=outputs["grid"],
        holdout_path=outputs["holdout"],
        failures_path=outputs["failures"],
        dashboard_path=outputs["dashboard"],
        certificate_path=outputs["certificate"],
    )

    checked_in = {
        "report": REPORT,
        "rows": ROWS,
        "grid": GRID,
        "holdout": HOLDOUT,
        "failures": FAILURES,
        "dashboard": DASHBOARD,
        "certificate": CERTIFICATE,
    }
    for key, generated_path in outputs.items():
        assert generated_path.read_text(encoding="utf-8") == checked_in[
            key
        ].read_text(encoding="utf-8")


def test_contact_law_threshold_root_outputs_are_canonical() -> None:
    assert set(ROOT_OUTPUT_NAMES) == {
        "contact_law_threshold_report.json",
        "contact_law_threshold_rows.csv",
        "contact_law_threshold_grid.csv",
        "contact_law_threshold_holdout.csv",
        "contact_law_threshold_failures.csv",
        "contact_law_threshold_dashboard.html",
        "contact_law_threshold_certificate.json",
    }
