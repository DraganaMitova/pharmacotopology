import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from pharmacotopology.folding_reference_loader import (  # noqa: E402
    load_folding_reference_dataset,
)
from pharmacotopology.folding_topology import (  # noqa: E402
    BROAD_FOLD_CLASSES,
    normalize_sequence,
    run_folding_topology_benchmark,
)
from render_folding_benchmark_dashboard import (  # noqa: E402
    render_folding_benchmark_dashboard,
)
from run_folding_topology_benchmark import (  # noqa: E402
    write_certificate_output,
    write_confusion_matrix_output,
    write_failures_output,
    write_folding_benchmark_outputs,
)


REAL_10_PATH = ROOT / "data" / "folding_benchmarks_real_10.locked.json"
EXPECTED_SOURCES = ["RCSB_PDB", "CATH", "DisProt"]


def test_tracked_real_10_dataset_is_locked_external_and_balanced() -> None:
    payload = json.loads(REAL_10_PATH.read_text(encoding="utf-8"))
    rows = payload["references"]

    assert payload["benchmark_kind"] == "real_external_folding_topology_benchmark"
    assert payload["benchmark_size"] == 10
    assert payload["target_benchmark_size"] == 10
    assert payload["benchmark_sources"] == EXPECTED_SOURCES
    assert payload["locked_after_generation"] is True
    assert payload["no_retuning_flag"] is True
    assert payload["folding_problem_solved"] is False
    assert payload["lock_certificate"]["lock_blockers"] == []
    assert set(row["reference_fold_class"] for row in rows) == set(
        BROAD_FOLD_CLASSES
    )

    per_class = {
        fold_class: sum(
            1 for row in rows if row["reference_fold_class"] == fold_class
        )
        for fold_class in BROAD_FOLD_CLASSES
    }
    assert set(per_class.values()) == {2}

    for row in rows:
        assert row["is_external_reference"] is True
        assert row["source_database"]
        assert row["source_accession"]
        assert row["reference_label_source"]
        assert row["curation_notes"]
        assert row["sequence_length"] == len(normalize_sequence(row["sequence"]))


def test_real_10_benchmark_writes_exact_report_and_sidecars(
    tmp_path: Path,
) -> None:
    dataset = load_folding_reference_dataset(REAL_10_PATH, require_external=True)
    comparisons = run_folding_topology_benchmark(dataset.references)

    report_path, csv_path = write_folding_benchmark_outputs(
        comparisons,
        tmp_path / "real_folding_10_report.json",
        tmp_path / "real_folding_10_rows.csv",
        reference_dataset_validation=dataset.validation,
        reference_dataset_metadata=dataset.metadata,
    )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    certificate_path = write_certificate_output(
        report,
        tmp_path / "real_folding_10_certificate.json",
    )
    failures_path = write_failures_output(
        comparisons,
        tmp_path / "real_folding_10_failures.csv",
    )
    confusion_path = write_confusion_matrix_output(
        comparisons,
        tmp_path / "real_folding_10_confusion_matrix.csv",
    )
    dashboard_path = render_folding_benchmark_dashboard(
        report_path,
        csv_path,
        tmp_path / "real_folding_10_dashboard.html",
    )

    assert report["benchmark_kind"] == "real_external_folding_topology_benchmark"
    assert report["benchmark_size"] == 10
    assert report["external_rows"] == 10
    assert report["benchmark_sources"] == EXPECTED_SOURCES
    assert report["evidence_readiness_summary"] == {
        "external_benchmark_attached": 10
    }
    assert report["match_count"] == 4
    assert report["mismatch_count"] == 6
    assert report["accuracy"] == 0.4
    assert report["locked_after_generation"] is True
    assert report["no_retuning_flag"] is True
    assert report["folding_problem_solved"] is False
    assert report["reference_dataset_validation"]["warnings"] == []

    certificate = json.loads(certificate_path.read_text(encoding="utf-8"))
    assert certificate["external_rows"] == 10
    assert certificate["folding_problem_solved"] is False
    assert certificate["lock_certificate"]["benchmark_sources"] == EXPECTED_SOURCES

    failure_rows = list(
        csv.DictReader(failures_path.read_text(encoding="utf-8").splitlines())
    )
    assert len(failure_rows) == 7
    assert {row["failure_reason"] for row in failure_rows} == {
        "fold_class_proxy_mismatch",
        "low_topology_similarity",
    }

    confusion_text = confusion_path.read_text(encoding="utf-8")
    assert "reference_fold_class" in confusion_text
    assert "alpha_rich" in confusion_text
    assert "disordered_flexible" in confusion_text

    dashboard = dashboard_path.read_text(encoding="utf-8")
    assert "External rows" in dashboard
    assert "Accuracy" in dashboard
    assert "Sources: RCSB_PDB, CATH, DisProt" in dashboard
