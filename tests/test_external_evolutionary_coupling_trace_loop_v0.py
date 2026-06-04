import csv
import json
from pathlib import Path
import sys
from typing import Optional

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from build_real_external_coupling_file_v0 import (  # noqa: E402
    build_real_external_coupling_file_v0,
)
from pharmacotopology.folding_coupling_negative_controls import (  # noqa: E402
    EXTERNAL_COUPLING_CONTROL_NAMES,
    generate_external_coupling_negative_controls,
)
from pharmacotopology.folding_external_coupling_importer import (  # noqa: E402
    import_external_coupling_dataset,
)
from pharmacotopology.folding_external_coupling_sources import (  # noqa: E402
    EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
)
from pharmacotopology.folding_external_coupling_trace_loop import (  # noqa: E402
    EXTERNAL_COUPLING_TRACE_LOOP_REPORT_KIND,
    ROOT_OUTPUT_NAMES,
    classify_external_probe_result,
    run_external_evolutionary_coupling_trace_loop_benchmark,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
REL_BENCHMARK_8 = Path("data/folding_real_coordinate_visual_8.locked.json")
REL_ORACLE_COUPLINGS = Path(
    "data/folding_real_coordinate_visual_8_couplings.locked.json"
)


def _external_constraint(row, *, rank: int, source_kind: str) -> dict[str, object]:
    separation_span = max(1, row.sequence_length - 24)
    separation = min(24 + ((rank - 1) % separation_span), row.sequence_length - 1)
    start_count = max(1, row.sequence_length - separation)
    i = ((rank - 1) // separation_span) % start_count + 1
    j = i + separation
    return {
        "row_id": row.row_id,
        "source_accession": row.source_accession,
        "constraint_id": f"external_{row.source_accession.replace(':', '_')}_{i}_{j}_{rank}",
        "i": i,
        "j": j,
        "alignment_i": i + 2,
        "alignment_j": j + 2,
        "sequence_separation": j - i,
        "normalized_separation": round((j - i) / row.sequence_length, 6),
        "confidence": round(0.98 - (rank % 20) * 0.001, 6),
        "raw_score": round(4.0 - (rank % 30) * 0.01, 6),
        "apc_corrected_score": round(3.7 - (rank % 30) * 0.01, 6),
        "rank": rank,
        "rank_fraction": round(rank / row.sequence_length, 6),
        "constraint_class": "external_dca_coupling",
        "source_kind": source_kind,
        "msa_source_kind": "pfam_or_uniref",
        "msa_sha256": f"sha256:test-{row.row_id}",
        "msa_depth": 4200,
        "effective_sequence_count": round(row.sequence_length * 6.2, 6),
        "effective_sequence_count_over_length": 6.2,
        "target_coverage": 0.94,
        "focus_sequence_mapping_confidence": 1.0,
        "top_L_couplings_available": True,
        "coordinate_truth_used_to_build_constraint": False,
        "native_truth_used_before_coupling_selection": False,
        "structure_model_used": False,
        "raw_sequence_exposed": False,
    }


def _write_external_fixture(
    path: Path,
    *,
    source_kind: str = "external_msa_dca_plmc_v1",
    tainted_row_index: Optional[int] = None,
) -> Path:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    constraints: list[dict[str, object]] = []
    for row_index, row in enumerate(rows):
        if row_index >= 4:
            continue
        for rank in range(1, row.sequence_length + 1):
            raw = _external_constraint(row, rank=rank, source_kind=source_kind)
            if tainted_row_index == row_index and rank == 1:
                raw["coordinate_truth_used_to_build_constraint"] = True
            constraints.append(raw)
    payload = {
        "layer_kind": "locked_safe_coupling_constraint_layer_v1",
        "batch_id": EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
        "constraint_kind": "safe_residue_pair_coupling_constraint_v1",
        "coupling_source_kind": source_kind,
        "external_evolutionary_couplings_used": True,
        "coordinate_truth_used_to_build_constraints": False,
        "native_truth_used_before_coupling_selection": False,
        "oracle_constraint_control": False,
        "raw_sequence_exposed": False,
        "source_benchmark_file": str(REL_BENCHMARK_8),
        "benchmark_row_ids_preregistered": [row.row_id for row in rows],
        "constraints": constraints,
    }
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def test_external_importer_preregisters_rows_and_quality_gates(tmp_path) -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    external_file = _write_external_fixture(tmp_path / "external.json")
    result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=external_file,
    )

    statuses = {status.row_id: status for status in result.row_statuses}
    assert len(statuses) == 8
    assert sum(
        status.row_external_status == "external_couplings_available"
        for status in statuses.values()
    ) == 4
    assert sum(
        status.row_external_status == "external_couplings_rejected_low_depth"
        for status in statuses.values()
    ) == 4
    assert result.dataset.external_evolutionary_couplings_used is True
    assert result.dataset.coordinate_truth_tainted is False
    assert result.dataset.native_truth_tainted is False
    assert result.dataset.oracle_constraint_control is False
    assert len(result.dataset.constraints) == sum(row.sequence_length for row in rows[:4])
    assert result.constraint_audits[0].native_contact_supported in {True, False}
    assert result.constraint_audits[0].benchmark_counts_as_false_positive in {
        True,
        False,
    }


def test_external_importer_rejects_disallowed_source_kind(tmp_path) -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    external_file = _write_external_fixture(
        tmp_path / "bad_source.json",
        source_kind="pdb_contact_map",
    )

    with pytest.raises(ValueError, match="rejected external coupling source kind"):
        import_external_coupling_dataset(
            rows=rows,
            external_coupling_file=external_file,
        )


def test_external_importer_requires_provenance_fields(tmp_path) -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    external_file = _write_external_fixture(tmp_path / "missing_field.json")
    payload = json.loads(external_file.read_text(encoding="utf-8"))
    del payload["constraints"][0]["apc_corrected_score"]
    external_file.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="missing required field"):
        import_external_coupling_dataset(
            rows=rows,
            external_coupling_file=external_file,
        )


def test_external_importer_rejects_duplicate_pairs(tmp_path) -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    external_file = _write_external_fixture(tmp_path / "duplicate.json")
    payload = json.loads(external_file.read_text(encoding="utf-8"))
    duplicate = dict(payload["constraints"][0])
    duplicate["constraint_id"] = "duplicate_pair"
    duplicate["rank"] = 999
    payload["constraints"].append(duplicate)
    external_file.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate external coupling pair"):
        import_external_coupling_dataset(
            rows=rows,
            external_coupling_file=external_file,
        )


def test_external_importer_marks_coordinate_taint_without_laundering(tmp_path) -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    external_file = _write_external_fixture(
        tmp_path / "tainted.json",
        tainted_row_index=0,
    )
    result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=external_file,
    )

    tainted = result.row_statuses[0]
    assert tainted.row_external_status == "external_couplings_rejected_coordinate_taint"
    assert result.any_coordinate_taint is True
    assert result.dataset.coordinate_truth_tainted is True
    assert result.dataset.oracle_constraint_control is True
    assert all(
        constraint.row_id != tainted.row_id for constraint in result.dataset.constraints
    )


def test_negative_controls_preserve_counts_and_clear_taint(tmp_path) -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    external_file = _write_external_fixture(tmp_path / "external.json")
    result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=external_file,
    )
    controls = generate_external_coupling_negative_controls(
        rows=rows,
        dataset=result.dataset,
    )

    assert set(controls) == set(EXTERNAL_COUPLING_CONTROL_NAMES)
    for control in controls.values():
        assert control.constraint_count == len(result.dataset.constraints)
        assert control.dataset.coordinate_truth_tainted is False
        assert control.dataset.native_truth_tainted is False
        assert control.dataset.raw_sequence_exposed is False


def test_external_trace_loop_runner_writes_claim_locked_outputs(tmp_path) -> None:
    external_file = _write_external_fixture(tmp_path / "external.json")
    outputs = {
        "report": tmp_path / "external_coupling_trace_loop_report.json",
        "certificate": tmp_path / "external_coupling_trace_loop_certificate.json",
        "selectors": tmp_path / "external_coupling_trace_loop_selectors.csv",
        "selected": tmp_path / "external_coupling_trace_loop_selected_events.csv",
        "controls": tmp_path / "external_coupling_trace_loop_controls.csv",
        "row_status": tmp_path / "external_coupling_trace_loop_row_status.csv",
        "dashboard": tmp_path / "external_coupling_trace_loop_dashboard.html",
    }
    paths = run_external_evolutionary_coupling_trace_loop_benchmark(
        benchmark_file=REL_BENCHMARK_8,
        external_coupling_file=external_file,
        oracle_coupling_file=REL_ORACLE_COUPLINGS,
        report_path=outputs["report"],
        certificate_path=outputs["certificate"],
        selectors_path=outputs["selectors"],
        selected_events_path=outputs["selected"],
        controls_path=outputs["controls"],
        row_status_path=outputs["row_status"],
        dashboard_path=outputs["dashboard"],
    )
    report = json.loads(outputs["report"].read_text(encoding="utf-8"))
    certificate = json.loads(outputs["certificate"].read_text(encoding="utf-8"))
    selectors = list(
        csv.DictReader(outputs["selectors"].read_text(encoding="utf-8").splitlines())
    )
    controls = list(
        csv.DictReader(outputs["controls"].read_text(encoding="utf-8").splitlines())
    )
    row_status = list(
        csv.DictReader(outputs["row_status"].read_text(encoding="utf-8").splitlines())
    )
    dashboard = outputs["dashboard"].read_text(encoding="utf-8")

    assert tuple(path.name for path in paths) == ROOT_OUTPUT_NAMES
    assert report["report_kind"] == EXTERNAL_COUPLING_TRACE_LOOP_REPORT_KIND
    assert report["batch_id"] == EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID
    assert report["result"] in {
        "external_channel_not_yet_supported",
        "external_channel_supported_in_v0",
    }
    assert report["external_couplings_available_rows"] == 4
    assert report["external_rows_rejected_low_depth"] == 4
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert report["claim_allowed"] is False
    assert certificate["claim_allowed"] is False
    assert len(selectors) == 8
    assert len(controls) == 8
    assert len(row_status) == 8
    assert "External Evolutionary Coupling Trace Loop V0" in dashboard


def test_real_external_builder_can_emit_zero_usable_rows_without_fabrication(tmp_path) -> None:
    output = tmp_path / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
    run_manifest = tmp_path / "external_coupling_target_manifest_v0.json"
    build_log = tmp_path / "external_coupling_build_log_v0.csv"
    paths = build_real_external_coupling_file_v0(
        benchmark_file=REL_BENCHMARK_8,
        raw_external_coupling_file=None,
        output=output,
        run_manifest_output=run_manifest,
        build_log_output=build_log,
        target_manifest_path=Path("data/external_coupling_target_manifest_v0.locked.json"),
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    manifest = json.loads(run_manifest.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(build_log.read_text(encoding="utf-8").splitlines()))

    assert paths == (output, run_manifest, build_log)
    assert payload["batch_id"] == EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID
    assert payload["reject_duplicate_coupling_pairs"] is True
    assert payload["constraints"] == []
    assert len(manifest) == 8
    assert len(rows) == 8
    assert {row["external_coupling_status"] for row in rows} == {
        "external_couplings_rejected_no_sequence_mapping"
    }
    assert all(row["duplicate_count_dropped"] == "0" for row in rows)


def test_external_probe_result_classification() -> None:
    assert (
        classify_external_probe_result(
            available_rows=0,
            external_real_beats_physical=True,
            external_real_beats_matched_controls=True,
        )
        == "insufficient_external_signal"
    )
    assert (
        classify_external_probe_result(
            available_rows=4,
            external_real_beats_physical=True,
            external_real_beats_matched_controls=False,
        )
        == "external_channel_not_yet_supported"
    )
    assert (
        classify_external_probe_result(
            available_rows=4,
            external_real_beats_physical=True,
            external_real_beats_matched_controls=True,
        )
        == "external_channel_supported_in_v0"
    )
