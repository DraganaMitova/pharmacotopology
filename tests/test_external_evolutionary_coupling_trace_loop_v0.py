import csv
import json
from pathlib import Path
import sys
from types import SimpleNamespace
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
from build_query_centered_hmmer_plmc_couplings_v0 import (  # noqa: E402
    DEFAULT_MINIMUM_SEQUENCE_SEPARATION,
    build_query_centered_hmmer_plmc_couplings_v0,
)
from build_query_centered_hmmer_apc_couplings_v0 import (  # noqa: E402
    build_query_centered_hmmer_apc_couplings_v0,
)
from build_query_centered_hmmer_plmc_apc_consensus_couplings_v0 import (  # noqa: E402
    build_query_centered_hmmer_plmc_apc_consensus_couplings_v0,
)
from pharmacotopology.folding_coupling_negative_controls import (  # noqa: E402
    EXTERNAL_ADVERSARIAL_CALIBRATED_CONTROL_NAMES,
    EXTERNAL_COUPLING_CONTROL_NAMES,
    generate_adversarial_calibrated_external_coupling_controls,
    generate_external_coupling_negative_controls,
)
from pharmacotopology import folding_coupling_nucleus_selector as selector_module  # noqa: E402
from pharmacotopology.folding_evolutionary_constraints import (  # noqa: E402
    COUPLING_CONSTRAINT_KIND,
    CouplingClosureAssessment,
    EVOLUTIONARY_COUPLING_LAYER_KIND,
    load_coupling_dataset,
)
from pharmacotopology.folding_nucleus_closure_search import (  # noqa: E402
    NucleusClosureEvent,
)
from pharmacotopology.folding_physical_selection import (  # noqa: E402
    build_active_physical_context,
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
    _build_multiscale_physical_contexts,
    _run_multiscale_critical_boundary_selector,
    _run_multiscale_future_preserved_selector,
    classify_external_probe_result,
    run_external_evolutionary_coupling_trace_loop_benchmark,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
LOCKED_EXTERNAL_COUPLINGS = (
    ROOT
    / "data"
    / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
)
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


def test_locked_external_hmmer_plmc_artifact_covers_all_rows_without_taint() -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=LOCKED_EXTERNAL_COUPLINGS,
    )
    statuses = {status.source_accession: status for status in result.row_statuses}
    protein_g = statuses["1PGA:A"]

    assert len(result.dataset.constraints) == 1139
    assert sum(
        status.row_external_status == "external_couplings_available"
        for status in result.row_statuses
    ) == 8
    assert protein_g.row_external_status == "external_couplings_available"
    assert protein_g.raw_constraint_count == 56
    assert protein_g.accepted_constraint_count == 56
    assert protein_g.target_coverage == 0.910714
    assert protein_g.focus_sequence_mapping_confidence == 1.0
    assert protein_g.effective_sequence_count_over_length == 8.785714
    assert result.dataset.external_evolutionary_couplings_used is True
    assert result.dataset.coordinate_truth_tainted is False
    assert result.dataset.native_truth_tainted is False
    assert result.dataset.structure_model_tainted is False
    assert result.dataset.oracle_constraint_control is False


def test_locked_boundary_field_replacement_probe_adds_zero_false_frontier() -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    dataset = load_coupling_dataset(LOCKED_EXTERNAL_COUPLINGS)
    context = selector_module.build_coupling_nucleus_context(
        rows=rows,
        coupling_dataset=dataset,
    )
    terminal = (
        selector_module.select_coupling_trace_loop_terminal_bridge_expanded_events(
            context
        )
    )
    replacement = (
        selector_module.select_coupling_trace_loop_boundary_field_replacement_probe_events(
            context
        )
    )
    terminal_metric = selector_module.selector_metrics(
        context,
        selector_name="terminal",
        selected_events=terminal,
    )
    replacement_metric = selector_module.selector_metrics(
        context,
        selector_name="replacement",
        selected_events=replacement,
    )
    terminal_ids = {event.event_id for event in terminal}
    added = tuple(
        event for event in replacement if event.event_id not in terminal_ids
    )

    assert (
        replacement_metric.selected_event_count
        == terminal_metric.selected_event_count
    )
    assert replacement_metric.false_nucleus_rate == 0.0
    assert replacement_metric.contact_cluster_precision == 0.167806
    assert replacement_metric.long_range_contact_recall == 0.350282
    assert (
        replacement_metric.long_range_contact_recall
        > terminal_metric.long_range_contact_recall
    )
    assert {event.event_id for event in added} == {
        "819f28bce14e526a",
        "de7a4160201080a9",
        "5c4251d4a0d1c221",
        "6dc052d830b6d8f3",
        "687f5226dff29418",
        "716dac20bd4eedb3",
    }
    assert all(event.native_contact_count_after_scoring > 0 for event in added)


def test_locked_macro_scale_future_preserved_selector_lifts_zero_false_recall() -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    dataset = load_coupling_dataset(LOCKED_EXTERNAL_COUPLINGS)
    physical_context = build_active_physical_context(
        rows,
        segment_length=20,
        segment_stride=4,
    )
    context = selector_module.build_coupling_nucleus_context(
        rows=rows,
        coupling_dataset=dataset,
        physical_context=physical_context,
    )
    selected = (
        selector_module.select_coupling_trace_loop_macro_scale_future_preserved_events(
            context
        )
    )
    metric = selector_module.selector_metrics(
        context,
        selector_name="macro_scale_future_preserved",
        selected_events=selected,
    )
    controls = generate_external_coupling_negative_controls(
        rows=rows,
        dataset=dataset,
    )
    adversarial_controls = generate_adversarial_calibrated_external_coupling_controls(
        rows=rows,
        dataset=dataset,
    )
    control_metrics = []
    for name, control in controls.items():
        control_context = selector_module.build_coupling_nucleus_context(
            rows=rows,
            coupling_dataset=control.dataset,
            physical_context=physical_context,
        )
        control_selected = (
            selector_module.select_coupling_trace_loop_macro_scale_future_preserved_events(
                control_context
            )
        )
        control_metrics.append(
            selector_module.selector_metrics(
                control_context,
                selector_name=f"macro_control_{name}",
                selected_events=control_selected,
            )
        )
    adversarial_metrics = []
    for name, control in adversarial_controls.items():
        control_context = selector_module.build_coupling_nucleus_context(
            rows=rows,
            coupling_dataset=control.dataset,
            physical_context=physical_context,
        )
        control_selected = (
            selector_module.select_coupling_trace_loop_macro_scale_future_preserved_events(
                control_context
            )
        )
        adversarial_metrics.append(
            selector_module.selector_metrics(
                control_context,
                selector_name=f"macro_adversarial_{name}",
                selected_events=control_selected,
            )
        )
    selected_control_metrics = [
        item for item in [*control_metrics, *adversarial_metrics]
        if item.selected_event_count
    ]

    assert metric.selected_event_count == 28
    assert metric.false_nucleus_rate == 0.0
    assert metric.contact_cluster_precision == 0.043505
    assert metric.long_range_contact_recall == 0.564435
    assert all(event.native_contact_count_after_scoring > 0 for event in selected)
    assert metric.contact_cluster_precision > max(
        item.contact_cluster_precision for item in selected_control_metrics
    )
    assert metric.long_range_contact_recall > max(
        item.long_range_contact_recall for item in selected_control_metrics
    )
    assert all(
        metric.false_nucleus_rate <= item.false_nucleus_rate
        for item in selected_control_metrics
    )


def test_locked_multiscale_future_preserved_selector_extends_recall_frontier() -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    dataset = load_coupling_dataset(LOCKED_EXTERNAL_COUPLINGS)
    run = _run_multiscale_future_preserved_selector(
        rows=rows,
        dataset=dataset,
        selector_name="external_multiscale_future_preserved",
        control_kind="external_real_multiscale_future_preserved",
        physical_contexts=_build_multiscale_physical_contexts(rows),
    )

    assert run.metric.selected_event_count == 114
    assert run.metric.false_nucleus_rate == 0.0
    assert run.metric.contact_cluster_precision == 0.03814
    assert run.metric.long_range_contact_recall == 0.745269
    assert all(
        event.native_contact_count_after_scoring > 0
        for event in run.selected_events
    )


def test_locked_multiscale_critical_boundary_self_selects_without_static_cap() -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    dataset = load_coupling_dataset(LOCKED_EXTERNAL_COUPLINGS)
    run = _run_multiscale_critical_boundary_selector(
        rows=rows,
        dataset=dataset,
        selector_name="external_multiscale_critical_boundary",
        control_kind="external_real_multiscale_critical_boundary",
        physical_contexts=_build_multiscale_physical_contexts(rows),
    )

    assert run.metric.selected_event_count == 114
    assert run.metric.false_nucleus_rate == 0.026786
    assert run.metric.contact_cluster_precision == 0.045534
    assert run.metric.long_range_contact_recall == 0.716901
    assert sum(
        event.native_contact_count_after_scoring == 0
        for event in run.selected_events
    ) == 3


def test_coupling_dataset_reuses_constraint_row_grouping() -> None:
    dataset = load_coupling_dataset(LOCKED_EXTERNAL_COUPLINGS)

    first = dataset.constraints_by_row_id()
    second = dataset.constraints_by_row_id()

    assert first is second
    assert set(first) == {
        row.row_id for row in load_real_coordinate_visual_rows(BENCHMARK_8)
    }


def test_query_centered_hmmer_plmc_builder_filters_local_pairs_before_top_l(
    tmp_path,
) -> None:
    row = load_real_coordinate_visual_rows(BENCHMARK_8)[0]
    focus_fasta = tmp_path / "focus.fasta"
    focus_fasta.write_text(f">synthetic_focus\n{row.sequence}\n", encoding="utf-8")
    base_external = tmp_path / "base_external.json"
    base_external.write_text(
        json.dumps(
            {
                "external_evolutionary_couplings_used": True,
                "coordinate_truth_used_to_build_constraints": False,
                "native_truth_used_before_coupling_selection": False,
                "oracle_constraint_control": False,
                "constraints": [],
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    plmc_rows: list[str] = []
    for index in range(1, min(row.sequence_length, 20)):
        plmc_rows.append(f"{index} A {index + 1} C 0 {1000 - index}.0")
    long_range_count = 0
    for i in range(1, row.sequence_length + 1):
        for j in range(i + DEFAULT_MINIMUM_SEQUENCE_SEPARATION, row.sequence_length + 1):
            long_range_count += 1
            plmc_rows.append(f"{i} A {j} C 0 {100 - long_range_count * 0.01:.6f}")
            if long_range_count >= row.sequence_length:
                break
        if long_range_count >= row.sequence_length:
            break
    plmc_couplings = tmp_path / "focus_plmc.couplings"
    plmc_couplings.write_text("\n".join(plmc_rows) + "\n", encoding="utf-8")

    output = build_query_centered_hmmer_plmc_couplings_v0(
        benchmark_file=BENCHMARK_8,
        base_external_coupling_file=base_external,
        row_id=row.row_id,
        focus_fasta_file=focus_fasta,
        plmc_couplings_file=plmc_couplings,
        output=tmp_path / "query_centered_external.json",
        msa_depth=2000,
        effective_sequence_count=2000.0,
        hmmer_job_id="synthetic_hmmer",
        hmmer_iteration_id="iter3",
    )

    payload = json.loads(output.read_text(encoding="utf-8"))
    constraints = payload["constraints"]

    assert len(constraints) == row.sequence_length
    assert all(
        constraint["sequence_separation"] >= DEFAULT_MINIMUM_SEQUENCE_SEPARATION
        for constraint in constraints
    )
    assert constraints[0]["raw_score"] == 99.99
    assert constraints[0]["minimum_sequence_separation"] == (
        DEFAULT_MINIMUM_SEQUENCE_SEPARATION
    )
    assert (
        payload[
            f"hmmer_query_centered_{row.source_accession.replace(':', '_')}_minimum_sequence_separation"
        ]
        == DEFAULT_MINIMUM_SEQUENCE_SEPARATION
    )


def test_query_centered_hmmer_apc_builder_emits_external_long_range_channel(
    tmp_path,
) -> None:
    rows = tuple(load_real_coordinate_visual_rows(BENCHMARK_8))
    row = next(row for row in rows if row.source_accession == "1PGA:A")
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    records = [f">synthetic_focus\n{row.sequence}"]
    for record_index in range(1, 96):
        sequence = "".join(
            alphabet[
                (
                    record_index * ((position % 5) + 1)
                    + position
                    + (position // 9)
                )
                % len(alphabet)
            ]
            for position in range(1, row.sequence_length + 1)
        )
        records.append(f">synthetic_{record_index}\n{sequence}")
    focus_fasta = tmp_path / "focus_apc.fasta"
    focus_fasta.write_text("\n".join(records) + "\n", encoding="utf-8")

    base_external = tmp_path / "base_external.json"
    base_external.write_text(
        json.dumps(
            {
                "layer_kind": EVOLUTIONARY_COUPLING_LAYER_KIND,
                "batch_id": EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
                "constraint_kind": COUPLING_CONSTRAINT_KIND,
                "coupling_source_kind": "external_uniref_msa_dca_v1",
                "external_evolutionary_couplings_used": True,
                "coordinate_truth_used_to_build_constraints": False,
                "native_truth_used_before_coupling_selection": False,
                "oracle_constraint_control": False,
                "raw_sequence_exposed": False,
                "source_benchmark_file": str(REL_BENCHMARK_8),
                "benchmark_row_ids_preregistered": [row.row_id for row in rows],
                "constraints": [],
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    output = build_query_centered_hmmer_apc_couplings_v0(
        benchmark_file=BENCHMARK_8,
        base_external_coupling_file=base_external,
        row_id=row.row_id,
        focus_fasta_file=focus_fasta,
        output=tmp_path / "query_centered_apc_external.json",
        msa_depth=400,
        effective_sequence_count=400.0,
        max_records=96,
    )

    payload = json.loads(output.read_text(encoding="utf-8"))
    constraints = payload["constraints"]

    assert len(constraints) == row.sequence_length
    assert all(
        constraint["sequence_separation"] >= DEFAULT_MINIMUM_SEQUENCE_SEPARATION
        for constraint in constraints
    )
    assert all(
        constraint["constraint_class"]
        == "external_query_centered_hmmer_mi_apc_coupling"
        for constraint in constraints
    )
    assert all(
        constraint["source_kind"] == "external_uniref_msa_dca_v1"
        for constraint in constraints
    )
    assert all(
        constraint["coordinate_truth_used_to_build_constraint"] is False
        and constraint["native_truth_used_before_coupling_selection"] is False
        and constraint["structure_model_used"] is False
        for constraint in constraints
    )
    assert constraints[0]["minimum_sequence_separation"] == (
        DEFAULT_MINIMUM_SEQUENCE_SEPARATION
    )
    assert constraints[0]["apc_focus_id"] == "synthetic_focus"
    result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=output,
    )
    statuses = {status.source_accession: status for status in result.row_statuses}
    assert statuses["1PGA:A"].row_external_status == "external_couplings_available"
    assert result.dataset.coordinate_truth_tainted is False
    assert result.dataset.native_truth_tainted is False
    assert result.dataset.structure_model_tainted is False


def test_query_centered_plmc_apc_consensus_reweights_plmc_spine(
    tmp_path,
) -> None:
    rows = tuple(load_real_coordinate_visual_rows(BENCHMARK_8))
    row = next(row for row in rows if row.source_accession == "1PGA:A")
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    records = [f">synthetic_focus\n{row.sequence}"]
    for record_index in range(1, 96):
        sequence = "".join(
            alphabet[
                (
                    record_index * ((position % 5) + 1)
                    + position
                    + (position // 9)
                )
                % len(alphabet)
            ]
            for position in range(1, row.sequence_length + 1)
        )
        records.append(f">synthetic_{record_index}\n{sequence}")
    focus_fasta = tmp_path / "focus_apc.fasta"
    focus_fasta.write_text("\n".join(records) + "\n", encoding="utf-8")
    base_constraints = [
        _external_constraint(
            row,
            rank=rank,
            source_kind="external_msa_dca_plmc_v1",
        )
        for rank in range(1, row.sequence_length + 1)
    ]
    base_external = tmp_path / "base_external.json"
    base_external.write_text(
        json.dumps(
            {
                "layer_kind": EVOLUTIONARY_COUPLING_LAYER_KIND,
                "batch_id": EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
                "constraint_kind": COUPLING_CONSTRAINT_KIND,
                "coupling_source_kind": "external_msa_dca_plmc_v1",
                "external_evolutionary_couplings_used": True,
                "coordinate_truth_used_to_build_constraints": False,
                "native_truth_used_before_coupling_selection": False,
                "oracle_constraint_control": False,
                "raw_sequence_exposed": False,
                "source_benchmark_file": str(REL_BENCHMARK_8),
                "benchmark_row_ids_preregistered": [row.row_id for row in rows],
                "constraints": base_constraints,
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    output = build_query_centered_hmmer_plmc_apc_consensus_couplings_v0(
        benchmark_file=BENCHMARK_8,
        base_external_coupling_file=base_external,
        row_id=row.row_id,
        focus_fasta_file=focus_fasta,
        output=tmp_path / "query_centered_consensus_external.json",
        consensus_weight=0.75,
        max_records=96,
    )

    payload = json.loads(output.read_text(encoding="utf-8"))
    constraints = payload["constraints"]
    base_pairs = {(constraint["i"], constraint["j"]) for constraint in base_constraints}
    consensus_pairs = {
        (constraint["i"], constraint["j"]) for constraint in constraints
    }

    assert len(constraints) == row.sequence_length
    assert consensus_pairs <= base_pairs
    assert all(
        constraint["constraint_class"]
        == "external_query_centered_hmmer_plmc_apc_consensus_coupling"
        for constraint in constraints
    )
    assert all(
        constraint["plmc_apc_consensus_method"]
        == "plmc_plus_weighted_apc_agreement"
        for constraint in constraints
    )
    assert constraints[0]["plmc_apc_consensus_weight"] == 0.75
    assert all(
        constraint["coordinate_truth_used_to_build_constraint"] is False
        and constraint["native_truth_used_before_coupling_selection"] is False
        and constraint["structure_model_used"] is False
        for constraint in constraints
    )
    result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=output,
    )
    statuses = {status.source_accession: status for status in result.row_statuses}
    assert statuses["1PGA:A"].row_external_status == "external_couplings_available"
    assert result.dataset.oracle_constraint_control is False


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


def test_adversarial_calibrated_controls_repair_provenance_metadata(tmp_path) -> None:
    rows = load_real_coordinate_visual_rows(BENCHMARK_8)
    external_file = _write_external_fixture(tmp_path / "external.json")
    result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=external_file,
    )
    controls = generate_adversarial_calibrated_external_coupling_controls(
        rows=rows,
        dataset=result.dataset,
    )

    assert set(controls) == set(EXTERNAL_ADVERSARIAL_CALIBRATED_CONTROL_NAMES)
    for control in controls.values():
        assert control.constraint_count == len(result.dataset.constraints)
        assert control.dataset.coordinate_truth_tainted is False
        assert control.dataset.native_truth_tainted is False
        assert all(constraint.rank > 0 for constraint in control.dataset.constraints)
        assert all(
            constraint.apc_corrected_score == constraint.confidence
            for constraint in control.dataset.constraints
        )


def _event(
    event_id: str,
    *,
    contact_cluster_gain: float,
    segment_a_start: int,
    segment_b_start: int,
    secondary_structure_compatibility: float = 0.7,
) -> NucleusClosureEvent:
    return NucleusClosureEvent(
        row_id="row_1",
        source_accession="TEST:A",
        sequence_hash="test",
        sequence_length=96,
        event_id=event_id,
        segment_a_start=segment_a_start,
        segment_a_end=segment_a_start + 7,
        segment_b_start=segment_b_start,
        segment_b_end=segment_b_start + 7,
        sequence_span=segment_b_start - segment_a_start + 8,
        normalized_span=0.5,
        candidate_contact_count=64,
        contact_cluster_gain=contact_cluster_gain,
        secondary_structure_compatibility=secondary_structure_compatibility,
        hydrophobic_burial_gain=0.6,
        registry_support=0.6,
        loop_entropy_cost=0.1,
        geometry_violation_cost=0.0,
        frustration_cost=0.0,
        isolation_penalty=0.0,
        nucleus_score=0.6,
        closure_event_stability=0.6,
        native_contact_count_after_scoring=0,
        native_long_range_contact_count_after_scoring=0,
        native_label_attached_after_event_generation=True,
    )


def _assessment(
    event: NucleusClosureEvent,
    *,
    direct_support_score: float,
    future_preservation_score: float,
    blocked_future_pressure: float,
) -> CouplingClosureAssessment:
    return CouplingClosureAssessment(
        row_id=event.row_id,
        source_accession=event.source_accession,
        event_id=event.event_id,
        direct_coupling_count=2,
        direct_coupling_confidence=direct_support_score,
        direct_support_score=direct_support_score,
        future_coupling_count=8,
        future_preserved_count=6,
        future_preservation_score=future_preservation_score,
        blocked_future_count=1,
        blocked_future_confidence=blocked_future_pressure,
        blocked_future_pressure=blocked_future_pressure,
        coupling_selectivity_score=direct_support_score + future_preservation_score,
        constraint_pairs_total=10,
        coordinate_truth_used_to_build_constraints=False,
        native_truth_used_before_coupling_selection=False,
    )


def test_rank_consistent_recovery_gate_adds_supported_lower_cluster_event(
    monkeypatch,
) -> None:
    core = _event(
        "core",
        contact_cluster_gain=0.35,
        segment_a_start=1,
        segment_b_start=41,
    )
    recovered = _event(
        "recovered",
        contact_cluster_gain=0.33,
        segment_a_start=13,
        segment_b_start=53,
    )
    rejected = _event(
        "rejected",
        contact_cluster_gain=0.33,
        segment_a_start=25,
        segment_b_start=65,
    )
    context = SimpleNamespace(
        rows=(SimpleNamespace(row_id="row_1"),),
        assessment_by_event_id={
            core.event_id: _assessment(
                core,
                direct_support_score=0.72,
                future_preservation_score=0.80,
                blocked_future_pressure=0.02,
            ),
            recovered.event_id: _assessment(
                recovered,
                direct_support_score=0.70,
                future_preservation_score=0.70,
                blocked_future_pressure=0.04,
            ),
            rejected.event_id: _assessment(
                rejected,
                direct_support_score=0.63,
                future_preservation_score=0.95,
                blocked_future_pressure=0.04,
            ),
        },
        coupling_decoy_margin_by_event_id={
            core.event_id: 0.10,
            recovered.event_id: 0.31,
            rejected.event_id: 0.40,
        },
    )
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_core_expanded_events",
        lambda *args, **kwargs: (core, recovered, rejected),
    )

    selected = selector_module.select_coupling_trace_loop_rank_consistent_cluster_gated_events(
        context
    )

    assert tuple(event.event_id for event in selected) == ("core", "recovered")


def test_persistent_rank_consistent_gate_recovers_supported_trace_fragment(
    monkeypatch,
) -> None:
    core = _event(
        "core",
        contact_cluster_gain=0.35,
        segment_a_start=1,
        segment_b_start=41,
    )
    persistent = _event(
        "persistent",
        contact_cluster_gain=0.33,
        segment_a_start=13,
        segment_b_start=57,
    )
    neighbor_a = _event(
        "neighbor_a",
        contact_cluster_gain=0.325,
        segment_a_start=21,
        segment_b_start=65,
    )
    neighbor_b = _event(
        "neighbor_b",
        contact_cluster_gain=0.325,
        segment_a_start=29,
        segment_b_start=73,
    )
    isolated = _event(
        "isolated",
        contact_cluster_gain=0.33,
        segment_a_start=60,
        segment_b_start=88,
    )
    context = SimpleNamespace(
        rows=(SimpleNamespace(row_id="row_1"),),
        assessment_by_event_id={
            core.event_id: _assessment(
                core,
                direct_support_score=0.72,
                future_preservation_score=0.82,
                blocked_future_pressure=0.02,
            ),
            persistent.event_id: _assessment(
                persistent,
                direct_support_score=0.68,
                future_preservation_score=0.82,
                blocked_future_pressure=0.04,
            ),
            neighbor_a.event_id: _assessment(
                neighbor_a,
                direct_support_score=0.70,
                future_preservation_score=0.55,
                blocked_future_pressure=0.04,
            ),
            neighbor_b.event_id: _assessment(
                neighbor_b,
                direct_support_score=0.69,
                future_preservation_score=0.56,
                blocked_future_pressure=0.04,
            ),
            isolated.event_id: _assessment(
                isolated,
                direct_support_score=0.68,
                future_preservation_score=0.82,
                blocked_future_pressure=0.04,
            ),
        },
        coupling_decoy_margin_by_event_id={
            core.event_id: 0.10,
            persistent.event_id: -0.02,
            neighbor_a.event_id: 0.0,
            neighbor_b.event_id: 0.0,
            isolated.event_id: -0.02,
        },
    )
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_core_expanded_events",
        lambda *args, **kwargs: (
            core,
            persistent,
            neighbor_a,
            neighbor_b,
            isolated,
        ),
    )

    evidence = selector_module.trace_loop_persistence_evidence(
        persistent,
        context,
        (core, persistent, neighbor_a, neighbor_b, isolated),
    )
    selected = (
        selector_module.select_coupling_trace_loop_persistent_rank_consistent_cluster_gated_events(
            context
        )
    )

    assert evidence.trace_loop_persistence_score >= 0.70
    assert evidence.persistent_neighbor_count >= 2
    assert tuple(event.event_id for event in selected) == ("core", "persistent")


def test_score_margin_expanded_selector_adds_guarded_trace_candidate(
    monkeypatch,
) -> None:
    core = _event(
        "core",
        contact_cluster_gain=0.35,
        segment_a_start=1,
        segment_b_start=41,
    )
    expanded = _event(
        "expanded",
        contact_cluster_gain=0.48,
        segment_a_start=13,
        segment_b_start=57,
    )
    rejected = _event(
        "rejected",
        contact_cluster_gain=0.48,
        segment_a_start=25,
        segment_b_start=73,
    )
    weak_direct = _event(
        "weak_direct",
        contact_cluster_gain=0.48,
        segment_a_start=31,
        segment_b_start=81,
    )
    decoy = _event(
        "decoy",
        contact_cluster_gain=0.48,
        segment_a_start=37,
        segment_b_start=89,
    )
    context = SimpleNamespace(
        rows=(SimpleNamespace(row_id="row_1"),),
        competitive_events=(core, expanded, rejected, decoy),
        assessment_by_event_id={
            core.event_id: _assessment(
                core,
                direct_support_score=0.72,
                future_preservation_score=0.82,
                blocked_future_pressure=0.02,
            ),
            expanded.event_id: _assessment(
                expanded,
                direct_support_score=0.30,
                future_preservation_score=0.30,
                blocked_future_pressure=0.02,
            ),
            rejected.event_id: _assessment(
                rejected,
                direct_support_score=0.30,
                future_preservation_score=0.30,
                blocked_future_pressure=0.12,
            ),
            weak_direct.event_id: _assessment(
                weak_direct,
                direct_support_score=0.30,
                future_preservation_score=0.30,
                blocked_future_pressure=0.02,
            ),
            decoy.event_id: _assessment(
                decoy,
                direct_support_score=0.20,
                future_preservation_score=0.20,
                blocked_future_pressure=0.02,
            ),
        },
    )
    scores = {
        core.event_id: 0.50,
        expanded.event_id: 0.62,
        rejected.event_id: 0.62,
        weak_direct.event_id: 0.62,
        decoy.event_id: 0.40,
    }
    direct_evidence = {
        expanded.event_id: 2,
        rejected.event_id: 2,
        weak_direct.event_id: 1,
    }
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_persistent_rank_consistent_cluster_gated_events",
        lambda *args, **kwargs: (core,),
    )
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_events",
        lambda *args, **kwargs: (core, expanded, rejected, weak_direct),
    )
    monkeypatch.setattr(
        selector_module,
        "coupling_nucleus_score",
        lambda event, _context: scores[event.event_id],
    )
    monkeypatch.setattr(
        selector_module,
        "_direct_constraint_trace_evidence",
        lambda event, _context: {
            "direct_constraint_count": direct_evidence.get(event.event_id, 0),
            "direct_constraint_confidence_sum": 1.0,
            "direct_top_10pct_rank_count": 0,
        },
    )
    monkeypatch.setattr(
        selector_module,
        "decoy_distance",
        lambda _event, candidate: 0.0
        if candidate.event_id == decoy.event_id
        else 10.0,
    )

    selected = selector_module.select_coupling_trace_loop_score_margin_expanded_events(
        context
    )

    assert tuple(event.event_id for event in selected) == ("core", "expanded")


def test_boundary_continuity_expanded_selector_rescues_structured_low_cluster_candidate(
    monkeypatch,
) -> None:
    core = _event(
        "core",
        contact_cluster_gain=0.35,
        segment_a_start=1,
        segment_b_start=41,
    )
    rescued = _event(
        "rescued",
        contact_cluster_gain=0.31,
        segment_a_start=13,
        segment_b_start=57,
        secondary_structure_compatibility=0.58,
    )
    weak_shape = _event(
        "weak_shape",
        contact_cluster_gain=0.31,
        segment_a_start=25,
        segment_b_start=73,
        secondary_structure_compatibility=0.48,
    )
    overconfident = _event(
        "overconfident",
        contact_cluster_gain=0.31,
        segment_a_start=37,
        segment_b_start=89,
        secondary_structure_compatibility=0.70,
    )
    decoy = _event(
        "decoy",
        contact_cluster_gain=0.31,
        segment_a_start=49,
        segment_b_start=93,
    )
    context = SimpleNamespace(
        rows=(SimpleNamespace(row_id="row_1"),),
        competitive_events=(core, rescued, weak_shape, overconfident, decoy),
        assessment_by_event_id={
            event.event_id: _assessment(
                event,
                direct_support_score=0.42,
                future_preservation_score=0.72,
                blocked_future_pressure=0.03,
            )
            for event in (core, rescued, weak_shape, overconfident, decoy)
        },
    )
    scores = {
        core.event_id: 0.50,
        rescued.event_id: 0.54,
        weak_shape.event_id: 0.55,
        overconfident.event_id: 0.66,
        decoy.event_id: 0.32,
    }
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_score_margin_expanded_events",
        lambda *args, **kwargs: (core,),
    )
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_events",
        lambda *args, **kwargs: (core, rescued, weak_shape, overconfident),
    )
    monkeypatch.setattr(
        selector_module,
        "coupling_nucleus_score",
        lambda event, _context: scores[event.event_id],
    )
    monkeypatch.setattr(
        selector_module,
        "decoy_distance",
        lambda _event, candidate: 0.0
        if candidate.event_id == decoy.event_id
        else 10.0,
    )

    selected = (
        selector_module.select_coupling_trace_loop_boundary_continuity_expanded_events(
            context
        )
    )

    assert tuple(event.event_id for event in selected) == ("core", "rescued")


def test_edge_continuity_expanded_selector_rescues_modest_high_cluster_edge(
    monkeypatch,
) -> None:
    core = _event(
        "core",
        contact_cluster_gain=0.35,
        segment_a_start=1,
        segment_b_start=41,
    )
    rescued = _event(
        "rescued",
        contact_cluster_gain=0.44,
        segment_a_start=13,
        segment_b_start=57,
        secondary_structure_compatibility=0.60,
    )
    overconfident = _event(
        "overconfident",
        contact_cluster_gain=0.44,
        segment_a_start=25,
        segment_b_start=73,
        secondary_structure_compatibility=0.60,
    )
    low_margin = _event(
        "low_margin",
        contact_cluster_gain=0.44,
        segment_a_start=37,
        segment_b_start=89,
        secondary_structure_compatibility=0.60,
    )
    decoy = _event(
        "decoy",
        contact_cluster_gain=0.44,
        segment_a_start=49,
        segment_b_start=93,
    )
    context = SimpleNamespace(
        rows=(SimpleNamespace(row_id="row_1"),),
        competitive_events=(core, rescued, overconfident, low_margin, decoy),
        assessment_by_event_id={
            event.event_id: _assessment(
                event,
                direct_support_score=0.32,
                future_preservation_score=0.58,
                blocked_future_pressure=0.03,
            )
            for event in (core, rescued, overconfident, low_margin, decoy)
        },
    )
    scores = {
        core.event_id: 0.50,
        rescued.event_id: 0.45,
        overconfident.event_id: 0.62,
        low_margin.event_id: 0.39,
        decoy.event_id: 0.30,
    }
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_boundary_continuity_expanded_events",
        lambda *args, **kwargs: (core,),
    )
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_events",
        lambda *args, **kwargs: (core, rescued, overconfident, low_margin),
    )
    monkeypatch.setattr(
        selector_module,
        "coupling_nucleus_score",
        lambda event, _context: scores[event.event_id],
    )
    monkeypatch.setattr(
        selector_module,
        "decoy_distance",
        lambda _event, candidate: 0.0
        if candidate.event_id == decoy.event_id
        else 10.0,
    )

    selected = selector_module.select_coupling_trace_loop_edge_continuity_expanded_events(
        context
    )

    assert tuple(event.event_id for event in selected) == ("core", "rescued")


def test_pressure_release_expanded_selector_rescues_dense_blocked_trace(
    monkeypatch,
) -> None:
    core = _event(
        "core",
        contact_cluster_gain=0.35,
        segment_a_start=1,
        segment_b_start=41,
    )
    rescued = _event(
        "rescued",
        contact_cluster_gain=0.39,
        segment_a_start=13,
        segment_b_start=57,
        secondary_structure_compatibility=0.52,
    )
    low_pressure = _event(
        "low_pressure",
        contact_cluster_gain=0.39,
        segment_a_start=25,
        segment_b_start=73,
        secondary_structure_compatibility=0.52,
    )
    weak_density = _event(
        "weak_density",
        contact_cluster_gain=0.39,
        segment_a_start=37,
        segment_b_start=89,
        secondary_structure_compatibility=0.52,
    )
    decoy = _event(
        "decoy",
        contact_cluster_gain=0.39,
        segment_a_start=49,
        segment_b_start=93,
    )
    context = SimpleNamespace(
        rows=(SimpleNamespace(row_id="row_1"),),
        competitive_events=(core, rescued, low_pressure, weak_density, decoy),
        assessment_by_event_id={
            core.event_id: _assessment(
                core,
                direct_support_score=0.42,
                future_preservation_score=0.32,
                blocked_future_pressure=0.16,
            ),
            rescued.event_id: _assessment(
                rescued,
                direct_support_score=0.42,
                future_preservation_score=0.32,
                blocked_future_pressure=0.16,
            ),
            low_pressure.event_id: _assessment(
                low_pressure,
                direct_support_score=0.42,
                future_preservation_score=0.32,
                blocked_future_pressure=0.07,
            ),
            weak_density.event_id: _assessment(
                weak_density,
                direct_support_score=0.42,
                future_preservation_score=0.32,
                blocked_future_pressure=0.16,
            ),
            decoy.event_id: _assessment(
                decoy,
                direct_support_score=0.10,
                future_preservation_score=0.10,
                blocked_future_pressure=0.16,
            ),
        },
    )
    scores = {
        core.event_id: 0.50,
        rescued.event_id: 0.44,
        low_pressure.event_id: 0.44,
        weak_density.event_id: 0.44,
        decoy.event_id: 0.28,
    }
    direct_evidence = {
        rescued.event_id: {
            "direct_constraint_count": 3,
            "direct_constraint_confidence_sum": 1.20,
            "direct_top_10pct_rank_count": 1,
        },
        low_pressure.event_id: {
            "direct_constraint_count": 3,
            "direct_constraint_confidence_sum": 1.20,
            "direct_top_10pct_rank_count": 1,
        },
        weak_density.event_id: {
            "direct_constraint_count": 2,
            "direct_constraint_confidence_sum": 1.19,
            "direct_top_10pct_rank_count": 1,
        },
    }
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_edge_continuity_expanded_events",
        lambda *args, **kwargs: (core,),
    )
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_events",
        lambda *args, **kwargs: (core, rescued, low_pressure, weak_density),
    )
    monkeypatch.setattr(
        selector_module,
        "coupling_nucleus_score",
        lambda event, _context: scores[event.event_id],
    )
    monkeypatch.setattr(
        selector_module,
        "decoy_distance",
        lambda _event, candidate: 0.0
        if candidate.event_id == decoy.event_id
        else 10.0,
    )
    monkeypatch.setattr(
        selector_module,
        "_direct_constraint_trace_evidence",
        lambda event, _context: direct_evidence.get(
            event.event_id,
            {
                "direct_constraint_count": 0,
                "direct_constraint_confidence_sum": 0.0,
                "direct_top_10pct_rank_count": 0,
            },
        ),
    )

    selected = (
        selector_module.select_coupling_trace_loop_pressure_release_expanded_events(
            context
        )
    )

    assert tuple(event.event_id for event in selected) == ("core", "rescued")


def test_registry_extension_expanded_selector_rescues_residual_trace(
    monkeypatch,
) -> None:
    anchor = _event(
        "anchor",
        contact_cluster_gain=0.50,
        segment_a_start=17,
        segment_b_start=81,
        secondary_structure_compatibility=0.60,
    )
    registry = _event(
        "registry",
        contact_cluster_gain=0.49,
        segment_a_start=33,
        segment_b_start=81,
        secondary_structure_compatibility=0.55,
    )
    tail = _event(
        "tail",
        contact_cluster_gain=0.43,
        segment_a_start=81,
        segment_b_start=129,
        secondary_structure_compatibility=0.54,
    )
    weak_tail = _event(
        "weak_tail",
        contact_cluster_gain=0.43,
        segment_a_start=97,
        segment_b_start=145,
        secondary_structure_compatibility=0.54,
    )
    decoy = _event(
        "decoy",
        contact_cluster_gain=0.49,
        segment_a_start=129,
        segment_b_start=153,
    )
    context = SimpleNamespace(
        rows=(SimpleNamespace(row_id="row_1"),),
        competitive_events=(anchor, registry, tail, weak_tail, decoy),
        assessment_by_event_id={
            anchor.event_id: _assessment(
                anchor,
                direct_support_score=0.30,
                future_preservation_score=0.32,
                blocked_future_pressure=0.04,
            ),
            registry.event_id: _assessment(
                registry,
                direct_support_score=0.12,
                future_preservation_score=0.20,
                blocked_future_pressure=0.07,
            ),
            tail.event_id: _assessment(
                tail,
                direct_support_score=0.27,
                future_preservation_score=0.20,
                blocked_future_pressure=0.05,
            ),
            weak_tail.event_id: _assessment(
                weak_tail,
                direct_support_score=0.27,
                future_preservation_score=0.20,
                blocked_future_pressure=0.05,
            ),
            decoy.event_id: _assessment(
                decoy,
                direct_support_score=0.10,
                future_preservation_score=0.10,
                blocked_future_pressure=0.05,
            ),
        },
        coupling_decoy_margin_by_event_id={
            anchor.event_id: 0.0,
            registry.event_id: -0.02,
            tail.event_id: 0.16,
            weak_tail.event_id: 0.16,
            decoy.event_id: 0.0,
        },
    )
    scores = {
        anchor.event_id: 0.52,
        registry.event_id: 0.42,
        tail.event_id: 0.43,
        weak_tail.event_id: 0.43,
        decoy.event_id: 0.30,
    }
    direct_evidence = {
        tail.event_id: {
            "direct_constraint_count": 4,
            "direct_constraint_confidence_sum": 0.82,
            "direct_top_10pct_rank_count": 0,
        },
        weak_tail.event_id: {
            "direct_constraint_count": 3,
            "direct_constraint_confidence_sum": 0.82,
            "direct_top_10pct_rank_count": 0,
        },
    }
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_pressure_release_expanded_events",
        lambda *args, **kwargs: (anchor,),
    )
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_events",
        lambda *args, **kwargs: (anchor, registry, tail, weak_tail),
    )
    monkeypatch.setattr(
        selector_module,
        "coupling_nucleus_score",
        lambda event, _context: scores[event.event_id],
    )
    monkeypatch.setattr(
        selector_module,
        "decoy_distance",
        lambda _event, candidate: 0.0
        if candidate.event_id == decoy.event_id
        else 10.0,
    )
    monkeypatch.setattr(
        selector_module,
        "_direct_constraint_trace_evidence",
        lambda event, _context: direct_evidence.get(
            event.event_id,
            {
                "direct_constraint_count": 0,
                "direct_constraint_confidence_sum": 0.0,
                "direct_top_10pct_rank_count": 0,
            },
        ),
    )

    selected = (
        selector_module.select_coupling_trace_loop_registry_extension_expanded_events(
            context
        )
    )

    assert tuple(event.event_id for event in selected) == (
        "anchor",
        "registry",
        "tail",
    )


def test_terminal_bridge_expanded_selector_rescues_final_bridge_and_tail(
    monkeypatch,
) -> None:
    anchor = _event(
        "anchor",
        contact_cluster_gain=0.45,
        segment_a_start=33,
        segment_b_start=61,
        secondary_structure_compatibility=0.60,
    )
    bridge = _event(
        "bridge",
        contact_cluster_gain=0.37,
        segment_a_start=61,
        segment_b_start=89,
        secondary_structure_compatibility=0.66,
    )
    tail = _event(
        "tail",
        contact_cluster_gain=0.42,
        segment_a_start=25,
        segment_b_start=69,
        secondary_structure_compatibility=0.54,
    )
    weak_bridge = _event(
        "weak_bridge",
        contact_cluster_gain=0.37,
        segment_a_start=89,
        segment_b_start=117,
        secondary_structure_compatibility=0.64,
    )
    weak_tail = _event(
        "weak_tail",
        contact_cluster_gain=0.42,
        segment_a_start=45,
        segment_b_start=89,
        secondary_structure_compatibility=0.54,
    )
    context = SimpleNamespace(
        rows=(SimpleNamespace(row_id="row_1"),),
        competitive_events=(anchor, bridge, tail, weak_bridge, weak_tail),
        assessment_by_event_id={
            anchor.event_id: _assessment(
                anchor,
                direct_support_score=0.30,
                future_preservation_score=0.32,
                blocked_future_pressure=0.04,
            ),
            bridge.event_id: _assessment(
                bridge,
                direct_support_score=0.21,
                future_preservation_score=0.23,
                blocked_future_pressure=0.06,
            ),
            tail.event_id: _assessment(
                tail,
                direct_support_score=0.16,
                future_preservation_score=0.31,
                blocked_future_pressure=0.30,
            ),
            weak_bridge.event_id: _assessment(
                weak_bridge,
                direct_support_score=0.21,
                future_preservation_score=0.23,
                blocked_future_pressure=0.06,
            ),
            weak_tail.event_id: _assessment(
                weak_tail,
                direct_support_score=0.16,
                future_preservation_score=0.31,
                blocked_future_pressure=0.30,
            ),
        },
        coupling_decoy_margin_by_event_id={
            anchor.event_id: 0.0,
            bridge.event_id: 0.06,
            tail.event_id: 0.05,
            weak_bridge.event_id: 0.06,
            weak_tail.event_id: 0.05,
        },
    )
    scores = {
        anchor.event_id: 0.48,
        bridge.event_id: 0.32,
        tail.event_id: 0.36,
        weak_bridge.event_id: 0.32,
        weak_tail.event_id: 0.36,
    }
    selector_margins = {
        bridge.event_id: -0.005,
        tail.event_id: 0.06,
        weak_bridge.event_id: -0.005,
        weak_tail.event_id: 0.06,
    }
    direct_evidence = {
        bridge.event_id: {
            "direct_constraint_count": 3,
            "direct_constraint_confidence_sum": 0.71,
            "direct_top_10pct_rank_count": 0,
        },
        tail.event_id: {
            "direct_constraint_count": 3,
            "direct_constraint_confidence_sum": 0.51,
            "direct_top_10pct_rank_count": 0,
        },
        weak_bridge.event_id: {
            "direct_constraint_count": 3,
            "direct_constraint_confidence_sum": 0.69,
            "direct_top_10pct_rank_count": 0,
        },
        weak_tail.event_id: {
            "direct_constraint_count": 2,
            "direct_constraint_confidence_sum": 0.51,
            "direct_top_10pct_rank_count": 0,
        },
    }
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_registry_extension_expanded_events",
        lambda *args, **kwargs: (anchor,),
    )
    monkeypatch.setattr(
        selector_module,
        "select_coupling_trace_loop_events",
        lambda *args, **kwargs: (anchor, bridge, tail, weak_bridge, weak_tail),
    )
    monkeypatch.setattr(
        selector_module,
        "coupling_nucleus_score",
        lambda event, _context: scores[event.event_id],
    )
    monkeypatch.setattr(
        selector_module,
        "_selector_score_decoy_margin",
        lambda event, _context: selector_margins.get(event.event_id, 0.0),
    )
    monkeypatch.setattr(
        selector_module,
        "_direct_constraint_trace_evidence",
        lambda event, _context: direct_evidence.get(
            event.event_id,
            {
                "direct_constraint_count": 0,
                "direct_constraint_confidence_sum": 0.0,
                "direct_top_10pct_rank_count": 0,
            },
        ),
    )

    selected = (
        selector_module.select_coupling_trace_loop_terminal_bridge_expanded_events(
            context
        )
    )

    assert tuple(event.event_id for event in selected) == (
        "anchor",
        "bridge",
        "tail",
    )


def test_external_trace_loop_runner_writes_claim_locked_outputs(tmp_path) -> None:
    external_file = _write_external_fixture(tmp_path / "external.json")
    outputs = {
        "report": tmp_path / "external_coupling_trace_loop_report.json",
        "certificate": tmp_path / "external_coupling_trace_loop_certificate.json",
        "selectors": tmp_path / "external_coupling_trace_loop_selectors.csv",
        "selected": tmp_path / "external_coupling_trace_loop_selected_events.csv",
        "frontier": tmp_path / "external_coupling_trace_loop_frontier.csv",
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
        frontier_path=outputs["frontier"],
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
    frontier = list(
        csv.DictReader(outputs["frontier"].read_text(encoding="utf-8").splitlines())
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
    assert "external_margin_gated_beats_matched_controls" in report
    assert "external_top_rank_gated_beats_matched_controls" in report
    assert "external_core_expanded_beats_matched_controls" in report
    assert "external_cluster_gated_core_expanded_beats_matched_controls" in report
    assert "external_rank_consistent_cluster_gated_beats_matched_controls" in report
    assert (
        "external_rank_consistent_cluster_gated_beats_adversarial_calibrated_controls"
        in report
    )
    assert "external_rank_consistent_cluster_gated_probe_passed" in report
    assert (
        "external_persistent_rank_consistent_cluster_gated_probe_passed"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_selector_score_probe_passed"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_vs_control_nucleus_score_enrichment_ratio"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_recovered_event_count"
        in report
    )
    assert "external_score_margin_expanded_selected_event_count" in report
    assert "external_score_margin_expanded_added_event_count" in report
    assert (
        "external_score_margin_expanded_added_native_long_range_contact_count"
        in report
    )
    assert "external_score_margin_expanded_added_false_event_count" in report
    assert "external_score_margin_expanded_long_range_recall" in report
    assert (
        "external_score_margin_expanded_long_range_recall_delta_vs_persistent"
        in report
    )
    assert "external_score_margin_expanded_beats_matched_controls" in report
    assert (
        "external_score_margin_expanded_beats_adversarial_calibrated_controls"
        in report
    )
    assert report["external_score_margin_expanded_claim_allowed"] is False
    assert "external_boundary_continuity_expanded_selected_event_count" in report
    assert "external_boundary_continuity_expanded_added_event_count" in report
    assert (
        "external_boundary_continuity_expanded_added_native_long_range_contact_count"
        in report
    )
    assert "external_boundary_continuity_expanded_added_false_event_count" in report
    assert "external_boundary_continuity_expanded_long_range_recall" in report
    assert (
        "external_boundary_continuity_expanded_long_range_recall_delta_vs_score_margin"
        in report
    )
    assert "external_boundary_continuity_expanded_beats_matched_controls" in report
    assert (
        "external_boundary_continuity_expanded_beats_adversarial_calibrated_controls"
        in report
    )
    assert report["external_boundary_continuity_expanded_claim_allowed"] is False
    assert "external_edge_continuity_expanded_selected_event_count" in report
    assert "external_edge_continuity_expanded_added_event_count" in report
    assert (
        "external_edge_continuity_expanded_added_native_long_range_contact_count"
        in report
    )
    assert "external_edge_continuity_expanded_added_false_event_count" in report
    assert "external_edge_continuity_expanded_long_range_recall" in report
    assert (
        "external_edge_continuity_expanded_long_range_recall_delta_vs_boundary_continuity"
        in report
    )
    assert "external_edge_continuity_expanded_beats_matched_controls" in report
    assert (
        "external_edge_continuity_expanded_beats_adversarial_calibrated_controls"
        in report
    )
    assert report["external_edge_continuity_expanded_claim_allowed"] is False
    assert "external_pressure_release_expanded_selected_event_count" in report
    assert "external_pressure_release_expanded_added_event_count" in report
    assert (
        "external_pressure_release_expanded_added_native_long_range_contact_count"
        in report
    )
    assert "external_pressure_release_expanded_added_false_event_count" in report
    assert "external_pressure_release_expanded_long_range_recall" in report
    assert (
        "external_pressure_release_expanded_long_range_recall_delta_vs_edge_continuity"
        in report
    )
    assert "external_pressure_release_expanded_beats_matched_controls" in report
    assert (
        "external_pressure_release_expanded_beats_adversarial_calibrated_controls"
        in report
    )
    assert report["external_pressure_release_expanded_claim_allowed"] is False
    assert "external_registry_extension_expanded_selected_event_count" in report
    assert "external_registry_extension_expanded_added_event_count" in report
    assert (
        "external_registry_extension_expanded_added_native_long_range_contact_count"
        in report
    )
    assert "external_registry_extension_expanded_added_false_event_count" in report
    assert "external_registry_extension_expanded_long_range_recall" in report
    assert (
        "external_registry_extension_expanded_long_range_recall_delta_vs_pressure_release"
        in report
    )
    assert "external_registry_extension_expanded_beats_matched_controls" in report
    assert (
        "external_registry_extension_expanded_beats_adversarial_calibrated_controls"
        in report
    )
    assert report["external_registry_extension_expanded_claim_allowed"] is False
    assert "external_terminal_bridge_expanded_selected_event_count" in report
    assert "external_terminal_bridge_expanded_added_event_count" in report
    assert (
        "external_terminal_bridge_expanded_added_native_long_range_contact_count"
        in report
    )
    assert "external_terminal_bridge_expanded_added_false_event_count" in report
    assert "external_terminal_bridge_expanded_long_range_recall" in report
    assert (
        "external_terminal_bridge_expanded_long_range_recall_delta_vs_registry_extension"
        in report
    )
    assert "external_terminal_bridge_expanded_beats_matched_controls" in report
    assert (
        "external_terminal_bridge_expanded_beats_adversarial_calibrated_controls"
        in report
    )
    assert report["external_terminal_bridge_expanded_claim_allowed"] is False
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_row_count"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_false_candidate_count"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_candidate_count"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_matched_controls"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_row_count"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_matched_controls"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_candidate_count"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_adversarial_controls"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_row_count"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_adversarial_controls"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_signal_seen"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_matched_controls"
        in report
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_adversarial_controls"
        in report
    )
    assert "hard_adversarial_calibrated_probe_passed" in report
    assert (
        "external_rank_consistent_cluster_gated_native_positive_frontier_count"
        in report
    )
    assert (
        report["external_rank_consistent_cluster_gated_frontier_claim_allowed"]
        is False
    )
    assert "external_boundary_field_replacement_probe_added_event_count" in report
    assert "external_boundary_field_replacement_probe_long_range_recall" in report
    assert (
        "external_boundary_field_replacement_probe_long_range_recall_delta_vs_terminal_bridge"
        in report
    )
    assert report["external_boundary_field_replacement_probe_claim_allowed"] is False
    assert report["external_macro_scale_future_preserved_segment_length"] == 20
    assert report["external_macro_scale_future_preserved_segment_stride"] == 4
    assert "external_macro_scale_future_preserved_selected_event_count" in report
    assert "external_macro_scale_future_preserved_false_nucleus_rate" in report
    assert "external_macro_scale_future_preserved_cluster_precision" in report
    assert "external_macro_scale_future_preserved_long_range_recall" in report
    assert (
        "external_macro_scale_future_preserved_long_range_recall_delta_vs_boundary_replacement"
        in report
    )
    assert (
        "external_macro_scale_future_preserved_cluster_precision_margin_vs_matched_controls"
        in report
    )
    assert (
        "external_macro_scale_future_preserved_long_range_recall_margin_vs_matched_controls"
        in report
    )
    assert (
        "external_macro_scale_future_preserved_cluster_precision_margin_vs_adversarial_controls"
        in report
    )
    assert (
        "external_macro_scale_future_preserved_long_range_recall_margin_vs_adversarial_controls"
        in report
    )
    assert (
        "external_macro_scale_future_preserved_beats_matched_controls" in report
    )
    assert (
        "external_macro_scale_future_preserved_beats_adversarial_calibrated_controls"
        in report
    )
    assert report["external_macro_scale_future_preserved_claim_allowed"] is False
    assert "external_multiscale_future_preserved_segment_lengths" in report
    assert "external_multiscale_future_preserved_max_events_per_row" in report
    assert "external_multiscale_future_preserved_selected_event_count" in report
    assert "external_multiscale_future_preserved_false_nucleus_rate" in report
    assert "external_multiscale_future_preserved_cluster_precision" in report
    assert "external_multiscale_future_preserved_long_range_recall" in report
    assert (
        "external_multiscale_future_preserved_long_range_recall_delta_vs_macro"
        in report
    )
    assert (
        "external_multiscale_future_preserved_long_range_recall_margin_vs_matched_controls"
        in report
    )
    assert (
        "external_multiscale_future_preserved_long_range_recall_margin_vs_adversarial_controls"
        in report
    )
    assert "external_multiscale_future_preserved_beats_matched_controls" in report
    assert (
        "external_multiscale_future_preserved_beats_adversarial_calibrated_controls"
        in report
    )
    assert report["external_multiscale_future_preserved_claim_allowed"] is False
    assert report["external_multiscale_critical_boundary_kind"] == (
        "largest_internal_coherence_gap_v0"
    )
    assert (
        report["external_multiscale_critical_boundary_static_future_threshold_used"]
        is False
    )
    assert (
        report["external_multiscale_critical_boundary_static_event_cap_used"]
        is False
    )
    assert "external_multiscale_critical_boundary_selected_event_count" in report
    assert "external_multiscale_critical_boundary_false_event_count" in report
    assert "external_multiscale_critical_boundary_false_nucleus_rate" in report
    assert "external_multiscale_critical_boundary_cluster_precision" in report
    assert "external_multiscale_critical_boundary_long_range_recall" in report
    assert (
        "external_multiscale_critical_boundary_long_range_recall_delta_vs_tuned_multiscale"
        in report
    )
    assert (
        "external_multiscale_critical_boundary_long_range_recall_margin_vs_matched_controls"
        in report
    )
    assert "external_multiscale_critical_boundary_beats_matched_controls" in report
    assert (
        "external_multiscale_critical_boundary_beats_adversarial_calibrated_controls"
        in report
    )
    assert report["external_multiscale_critical_boundary_claim_allowed"] is False
    assert "external_terminal_bridge_replacement_frontier_count" in report
    assert (
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum"
        in report
    )
    assert (
        "external_terminal_bridge_replacement_frontier_external_count_delta_positive_count"
        in report
    )
    assert (
        "external_terminal_bridge_replacement_frontier_external_confidence_delta_positive_count"
        in report
    )
    assert (
        "external_terminal_bridge_replacement_frontier_external_count_delta_sum"
        in report
    )
    assert (
        "external_terminal_bridge_replacement_frontier_external_confidence_delta_sum"
        in report
    )
    assert (
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum_with_external_count_gain"
        in report
    )
    assert (
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum_with_external_confidence_gain"
        in report
    )
    assert (
        report[
            "external_terminal_bridge_replacement_frontier_external_count_delta_positive_count"
        ]
        <= report["external_terminal_bridge_replacement_frontier_count"]
    )
    assert (
        report[
            "external_terminal_bridge_replacement_frontier_external_confidence_delta_positive_count"
        ]
        <= report["external_terminal_bridge_replacement_frontier_count"]
    )
    assert (
        report["external_terminal_bridge_replacement_frontier_claim_allowed"]
        is False
    )
    assert report["external_margin_gated_claim_allowed"] is False
    assert report["external_top_rank_gated_claim_allowed"] is False
    assert report["external_core_expanded_claim_allowed"] is False
    assert report["external_cluster_gated_core_expanded_claim_allowed"] is False
    assert report["external_rank_consistent_cluster_gated_claim_allowed"] is False
    assert (
        report["external_persistent_rank_consistent_cluster_gated_claim_allowed"]
        is False
    )
    assert (
        report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_claim_allowed"
        ]
        is False
    )
    assert report["mechanism_discovery_claim_allowed"] is False
    assert report["folding_problem_solved"] is False
    assert report["claim_allowed"] is False
    assert certificate["claim_allowed"] is False
    assert (
        certificate[
            "external_rank_consistent_cluster_gated_native_positive_frontier_count"
        ]
        == report[
            "external_rank_consistent_cluster_gated_native_positive_frontier_count"
        ]
    )
    assert (
        certificate[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_adversarial_controls"
        ]
        == report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_adversarial_controls"
        ]
    )
    assert (
        certificate[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_adversarial_controls"
        ]
        == report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_adversarial_controls"
        ]
    )
    assert (
        certificate["external_score_margin_expanded_added_event_count"]
        == report["external_score_margin_expanded_added_event_count"]
    )
    assert (
        certificate["external_boundary_continuity_expanded_added_event_count"]
        == report["external_boundary_continuity_expanded_added_event_count"]
    )
    assert (
        certificate["external_edge_continuity_expanded_added_event_count"]
        == report["external_edge_continuity_expanded_added_event_count"]
    )
    assert (
        certificate["external_pressure_release_expanded_added_event_count"]
        == report["external_pressure_release_expanded_added_event_count"]
    )
    assert (
        certificate["external_registry_extension_expanded_added_event_count"]
        == report["external_registry_extension_expanded_added_event_count"]
    )
    assert (
        certificate["external_terminal_bridge_expanded_added_event_count"]
        == report["external_terminal_bridge_expanded_added_event_count"]
    )
    assert (
        certificate[
            "external_boundary_field_replacement_probe_long_range_recall_delta_vs_terminal_bridge"
        ]
        == report[
            "external_boundary_field_replacement_probe_long_range_recall_delta_vs_terminal_bridge"
        ]
    )
    assert (
        certificate["external_macro_scale_future_preserved_long_range_recall"]
        == report["external_macro_scale_future_preserved_long_range_recall"]
    )
    assert (
        certificate[
            "external_macro_scale_future_preserved_long_range_recall_delta_vs_boundary_replacement"
        ]
        == report[
            "external_macro_scale_future_preserved_long_range_recall_delta_vs_boundary_replacement"
        ]
    )
    assert (
        certificate[
            "external_macro_scale_future_preserved_long_range_recall_margin_vs_matched_controls"
        ]
        == report[
            "external_macro_scale_future_preserved_long_range_recall_margin_vs_matched_controls"
        ]
    )
    assert (
        certificate[
            "external_macro_scale_future_preserved_beats_matched_controls"
        ]
        == report["external_macro_scale_future_preserved_beats_matched_controls"]
    )
    assert (
        certificate["external_multiscale_future_preserved_long_range_recall"]
        == report["external_multiscale_future_preserved_long_range_recall"]
    )
    assert (
        certificate[
            "external_multiscale_future_preserved_long_range_recall_margin_vs_matched_controls"
        ]
        == report[
            "external_multiscale_future_preserved_long_range_recall_margin_vs_matched_controls"
        ]
    )
    assert (
        certificate["external_multiscale_critical_boundary_long_range_recall"]
        == report["external_multiscale_critical_boundary_long_range_recall"]
    )
    assert (
        certificate[
            "external_multiscale_critical_boundary_long_range_recall_delta_vs_tuned_multiscale"
        ]
        == report[
            "external_multiscale_critical_boundary_long_range_recall_delta_vs_tuned_multiscale"
        ]
    )
    assert (
        certificate["external_terminal_bridge_replacement_frontier_count"]
        == report["external_terminal_bridge_replacement_frontier_count"]
    )
    assert (
        certificate[
            "external_terminal_bridge_replacement_frontier_external_count_delta_positive_count"
        ]
        == report[
            "external_terminal_bridge_replacement_frontier_external_count_delta_positive_count"
        ]
    )
    assert (
        certificate[
            "external_terminal_bridge_replacement_frontier_external_confidence_delta_sum"
        ]
        == report[
            "external_terminal_bridge_replacement_frontier_external_confidence_delta_sum"
        ]
    )
    assert len(selectors) == 121
    assert len(controls) == 121
    assert len(frontier) == (
        report["external_terminal_bridge_replacement_frontier_count"]
        +
        report["external_rank_consistent_cluster_gated_native_positive_frontier_count"]
        + report[
            "external_persistent_rank_consistent_cluster_gated_recall_frontier_count"
        ]
    )
    replacement_frontier_rows = [
        row
        for row in frontier
        if row["frontier_kind"] == "terminal_bridge_replacement_frontier_v0"
    ]
    if replacement_frontier_rows:
        replacement_row = replacement_frontier_rows[0]
        assert "replacement_external_constraint_coverage_delta" in replacement_row
        assert "replacement_external_constraint_confidence_delta" in replacement_row
        assert "replacement_external_constraint_coverage_improved" in replacement_row
        assert "replacement_external_constraint_confidence_improved" in replacement_row
    assert len(row_status) == 8
    assert (
        "external_rank_consistent_cluster_gated_native_positive_frontier_count"
        in dashboard
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_recovered_event_count"
        in dashboard
    )
    assert "external_score_margin_expanded_added_event_count" in dashboard
    assert "external_boundary_field_replacement_probe_long_range_recall" in dashboard
    assert "external_macro_scale_future_preserved_long_range_recall" in dashboard
    assert (
        "external_macro_scale_future_preserved_long_range_recall_margin_vs_matched_controls"
        in dashboard
    )
    assert "external_multiscale_future_preserved_long_range_recall" in dashboard
    assert (
        "external_multiscale_future_preserved_long_range_recall_margin_vs_matched_controls"
        in dashboard
    )
    assert "external_multiscale_critical_boundary_long_range_recall" in dashboard
    assert (
        "external_multiscale_critical_boundary_long_range_recall_delta_vs_tuned_multiscale"
        in dashboard
    )
    assert "external_terminal_bridge_replacement_frontier_count" in dashboard
    assert (
        "external_terminal_bridge_replacement_frontier_external_count_delta_positive_count"
        in dashboard
    )
    assert (
        "external_score_margin_expanded_long_range_recall_delta_vs_persistent"
        in dashboard
    )
    assert "external_boundary_continuity_expanded_added_event_count" in dashboard
    assert (
        "external_boundary_continuity_expanded_long_range_recall_delta_vs_score_margin"
        in dashboard
    )
    assert "external_edge_continuity_expanded_added_event_count" in dashboard
    assert (
        "external_edge_continuity_expanded_long_range_recall_delta_vs_boundary_continuity"
        in dashboard
    )
    assert "external_pressure_release_expanded_added_event_count" in dashboard
    assert (
        "external_pressure_release_expanded_long_range_recall_delta_vs_edge_continuity"
        in dashboard
    )
    assert "external_registry_extension_expanded_added_event_count" in dashboard
    assert (
        "external_registry_extension_expanded_long_range_recall_delta_vs_pressure_release"
        in dashboard
    )
    assert "external_terminal_bridge_expanded_added_event_count" in dashboard
    assert (
        "external_terminal_bridge_expanded_long_range_recall_delta_vs_registry_extension"
        in dashboard
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count"
        in dashboard
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count_margin_vs_matched_controls"
        in dashboard
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_matched_controls"
        in dashboard
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_adversarial_controls"
        in dashboard
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_signal_seen"
        in dashboard
    )
    assert (
        "external_persistent_rank_consistent_cluster_gated_selector_score_probe_passed"
        in dashboard
    )
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
    assert payload["external_evolutionary_couplings_used"] is False
    assert payload["external_constraint_count"] == 0
    assert payload["result"] == "no_external_data_built"
    assert payload["constraints"] == []
    assert len(manifest) == 8
    assert len(rows) == 8
    assert {row["external_coupling_status"] for row in rows} == {
        "not_attempted_no_acquisition_pipeline"
    }
    assert all(row["duplicate_count_dropped"] == "0" for row in rows)


def test_external_probe_result_classification() -> None:
    assert (
        classify_external_probe_result(
            available_rows=0,
            external_real_beats_physical=True,
            external_real_beats_matched_controls=True,
        )
        == "no_external_data_built"
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
