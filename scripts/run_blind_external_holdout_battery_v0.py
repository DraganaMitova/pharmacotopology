from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.artifact_io import write_csv_rows  # noqa: E402
from pharmacotopology.folding_coupling_negative_controls import (  # noqa: E402
    EXTERNAL_COUPLING_CONTROL_NAMES,
    generate_external_coupling_negative_controls,
)
from pharmacotopology.folding_evolutionary_constraints import (  # noqa: E402
    CouplingDataset,
)
from pharmacotopology.folding_external_coupling_importer import (  # noqa: E402
    import_external_coupling_dataset,
)
from pharmacotopology.folding_external_coupling_trace_loop import (  # noqa: E402
    TraceLoopRun,
    _build_multiscale_physical_contexts,
    _run_multiscale_phase_aligned_critical_boundary_selector,
    _run_multiscale_phase_aligned_external_novelty_boundary_selector,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


BATTERY_REPORT_KIND = "blind_external_holdout_battery_v0"
BATTERY_CERTIFICATE_KIND = "blind_external_holdout_battery_certificate_v0"
DEFAULT_TARGET_MANIFEST = Path("data/blind_holdout_manifest_v0.locked.json")
DEFAULT_COUPLING_DIR = Path("data/blind_external_couplings_v0")
DEFAULT_OUTPUT_DIR = Path(
    "first_contact_clean_pharmacotopology_layer_run/blind_external_holdout_v0"
)
DEFAULT_SELECTOR = "external_multiscale_phase_aligned_external_novelty_boundary"
FRONTIER_SELECTOR = "external_multiscale_phase_aligned_external_novelty_frontier"
PRECISION_BOUNDARY_SELECTOR = (
    "external_multiscale_phase_aligned_critical_boundary"
)
SATURATION_BOUNDARY_SELECTOR = DEFAULT_SELECTOR
PASS_TARGET_WIN_RATE_MIN = 0.70


@dataclass(frozen=True)
class ExactContactMetrics:
    exact_contact_precision_top_L: float
    exact_contact_precision_L_over_2: float
    exact_long_range_contact_precision: float
    exact_long_range_contact_recall: float
    region_width_penalty: float
    region_pair_over_native_hit_ratio: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(value, 6)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return f"sha256:{digest.hexdigest()}"


def _current_commit_hash() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return ""
    return result.stdout.strip()


def _load_json(path: Path) -> Mapping[str, object]:
    parsed = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, Mapping):
        raise ValueError(f"{path} must contain a JSON object")
    return parsed


def _resolve_path(raw: object, *, base_dir: Path) -> Path:
    path = Path(str(raw))
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def _target_coupling_file(
    target: Mapping[str, object],
    *,
    coupling_dir: Path,
) -> Path:
    for field_name in (
        "external_coupling_file",
        "external_coupling_filename",
        "coupling_file",
    ):
        raw = target.get(field_name)
        if raw:
            path = Path(str(raw))
            if path.is_absolute():
                return path
            repo_relative = REPO_ROOT / path
            if repo_relative.exists():
                return repo_relative
            return coupling_dir / path
    target_id = str(target["target_id"])
    return coupling_dir / f"{target_id}.locked.json"


def _constraint_precision(
    constraints: Sequence[object],
    native_pairs: set[tuple[int, int]],
) -> float:
    if not constraints:
        return 0.0
    supported = sum(
        1
        for constraint in constraints
        if getattr(constraint, "pair")() in native_pairs
    )
    return _rounded(supported / len(constraints))


def _exact_contact_metrics(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
    run: TraceLoopRun,
) -> ExactContactMetrics:
    constraints_by_row = dataset.constraints_by_row_id()
    top_l_precisions: list[float] = []
    top_half_precisions: list[float] = []
    long_precisions: list[float] = []
    long_recalls: list[float] = []
    for row in rows:
        native_pairs = set(row.native_contact_pairs())
        native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
        constraints = tuple(constraints_by_row.get(row.row_id, ()))
        ranked = tuple(
            sorted(
                constraints,
                key=lambda constraint: (
                    constraint.rank or 10**9,
                    -constraint.confidence,
                    constraint.i,
                    constraint.j,
                ),
            )
        )
        top_l = ranked[: row.sequence_length]
        half_count = max(1, row.sequence_length // 2)
        top_half = ranked[:half_count]
        long_constraints = tuple(
            constraint
            for constraint in ranked
            if constraint.sequence_separation >= 24
        )
        supported_long = {
            constraint.pair()
            for constraint in long_constraints
            if constraint.pair() in native_long
        }
        top_l_precisions.append(_constraint_precision(top_l, native_pairs))
        top_half_precisions.append(_constraint_precision(top_half, native_pairs))
        long_precisions.append(_constraint_precision(long_constraints, native_pairs))
        long_recalls.append(
            _rounded(len(supported_long) / len(native_long))
            if native_long
            else 1.0
        )

    possible_region_pair_count = sum(
        len(event.candidate_region_pairs())
        for event in run.selected_events
    )
    native_hit_count = sum(
        event.native_contact_count_after_scoring
        for event in run.selected_events
    )
    cluster_precision = run.metric.contact_cluster_precision
    return ExactContactMetrics(
        exact_contact_precision_top_L=_rounded(mean(top_l_precisions))
        if top_l_precisions
        else 0.0,
        exact_contact_precision_L_over_2=_rounded(mean(top_half_precisions))
        if top_half_precisions
        else 0.0,
        exact_long_range_contact_precision=_rounded(mean(long_precisions))
        if long_precisions
        else 0.0,
        exact_long_range_contact_recall=_rounded(mean(long_recalls))
        if long_recalls
        else 0.0,
        region_width_penalty=_rounded(1.0 - cluster_precision),
        region_pair_over_native_hit_ratio=_rounded(
            possible_region_pair_count / native_hit_count
        )
        if native_hit_count
        else 0.0,
    )


def _metric_row(
    *,
    target_id: str,
    run: TraceLoopRun,
    exact_metrics: ExactContactMetrics,
    benchmark_file: Path,
    external_coupling_file: Path,
    benchmark_file_sha256: str,
    external_coupling_file_sha256: str,
    source_accessions: Sequence[str],
    control_name: str = "",
) -> dict[str, object]:
    metric = run.metric
    false_event_count = round(
        metric.false_nucleus_rate * metric.selected_event_count
    )
    row = {
        "target_id": target_id,
        "selector_name": run.selector_name,
        "control_name": control_name,
        "control_kind": run.control_kind,
        "benchmark_file": str(benchmark_file),
        "external_coupling_file": str(external_coupling_file),
        "benchmark_file_sha256": benchmark_file_sha256,
        "external_coupling_file_sha256": external_coupling_file_sha256,
        "source_accessions": ";".join(source_accessions),
        "constraint_count": run.constraint_count,
        "selected_event_count": metric.selected_event_count,
        "false_event_count": false_event_count,
        "false_nucleus_rate": metric.false_nucleus_rate,
        "cluster_precision": metric.contact_cluster_precision,
        "long_range_contact_recall": metric.long_range_contact_recall,
        "coupling_constraint_recall": metric.coupling_constraint_recall,
        "real_vs_decoy_coupling_enrichment_ratio": (
            metric.real_vs_decoy_coupling_enrichment_ratio
        ),
        "coordinate_truth_used_to_build_constraints": (
            metric.coordinate_truth_used_to_build_constraints
        ),
        "native_truth_used_before_coupling_selection": (
            metric.native_truth_used_before_coupling_selection
        ),
        "raw_sequence_exposed": metric.raw_sequence_exposed,
    }
    row.update(exact_metrics.to_dict())
    return row


def _run_real_frontier(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
    physical_contexts: Mapping[int, object],
    control_kind: str,
) -> tuple[TraceLoopRun, TraceLoopRun]:
    precision = _run_multiscale_phase_aligned_critical_boundary_selector(
        rows=rows,
        dataset=dataset,
        selector_name=PRECISION_BOUNDARY_SELECTOR,
        control_kind=control_kind,
        physical_contexts=physical_contexts,
    )
    saturation = _run_multiscale_phase_aligned_external_novelty_boundary_selector(
        rows=rows,
        dataset=dataset,
        selector_name=SATURATION_BOUNDARY_SELECTOR,
        control_kind=control_kind,
        physical_contexts=physical_contexts,
    )
    return precision, saturation


def _run_control_frontier(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
    physical_contexts: Mapping[int, object],
    control_name: str,
    control_kind: str,
) -> tuple[TraceLoopRun, TraceLoopRun]:
    precision = _run_multiscale_phase_aligned_critical_boundary_selector(
        rows=rows,
        dataset=dataset,
        selector_name=f"{control_name}_precision_boundary",
        control_kind=control_kind,
        physical_contexts=physical_contexts,
    )
    saturation = _run_multiscale_phase_aligned_external_novelty_boundary_selector(
        rows=rows,
        dataset=dataset,
        selector_name=f"{control_name}_saturation_boundary",
        control_kind=control_kind,
        physical_contexts=physical_contexts,
    )
    return precision, saturation


def _dominates(left: TraceLoopRun, right: TraceLoopRun) -> bool:
    left_metric = left.metric
    right_metric = right.metric
    return (
        left_metric.false_nucleus_rate <= right_metric.false_nucleus_rate
        and left_metric.long_range_contact_recall
        >= right_metric.long_range_contact_recall
        and left_metric.contact_cluster_precision
        > right_metric.contact_cluster_precision
    )


def _frontier_dominates(
    left_runs: Sequence[TraceLoopRun],
    right_runs: Sequence[TraceLoopRun],
) -> bool:
    return all(
        any(_dominates(left, right) for left in left_runs)
        for right in right_runs
    )


def _frontier_best_metric(
    runs: Sequence[TraceLoopRun],
    field_name: str,
) -> float:
    return _rounded(max(float(getattr(run.metric, field_name)) for run in runs))


def _write_json(path: Path, payload: Mapping[str, object]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _control_names(raw: Sequence[str]) -> tuple[str, ...]:
    if not raw or raw == ("all",):
        return tuple(EXTERNAL_COUPLING_CONTROL_NAMES)
    names = tuple(raw)
    unknown = sorted(set(names).difference(EXTERNAL_COUPLING_CONTROL_NAMES))
    if unknown:
        raise ValueError(f"unknown control name(s): {', '.join(unknown)}")
    return names


def _mean_field(rows: Sequence[Mapping[str, object]], field_name: str) -> float:
    values = [float(row[field_name]) for row in rows]
    return _rounded(mean(values)) if values else 0.0


def _max_field(rows: Sequence[Mapping[str, object]], field_name: str) -> float:
    values = [float(row[field_name]) for row in rows]
    return _rounded(max(values)) if values else 0.0


def _battery_classification(
    *,
    target_rows: Sequence[Mapping[str, object]],
    control_rows: Sequence[Mapping[str, object]],
    failure_rows: Sequence[Mapping[str, object]],
    oracle_taint_violation_count: int,
    target_manifest_frozen: bool,
    selector_frozen_after_manifest: bool,
) -> str:
    if (
        oracle_taint_violation_count
        or not target_manifest_frozen
        or not selector_frozen_after_manifest
    ):
        return "folding_claim_forbidden"
    if failure_rows:
        return "blind_batch_failed"
    if not target_rows or not control_rows:
        return "blind_batch_inconclusive"

    frontier_win_rate = _mean_field(target_rows, "frontier_control_win")
    scaffold_signal = frontier_win_rate >= PASS_TARGET_WIN_RATE_MIN
    exact_signal = (
        scaffold_signal
        and _mean_field(target_rows, "exact_contact_precision_top_L")
        > _max_field(control_rows, "exact_contact_precision_top_L")
        and _mean_field(target_rows, "exact_long_range_contact_precision")
        > _max_field(control_rows, "exact_long_range_contact_precision")
    )
    if exact_signal:
        return "blind_batch_exact_contact_signal_confirmed"
    if scaffold_signal:
        return "blind_batch_scaffold_signal_confirmed"
    return "blind_batch_failed"


def run_blind_external_holdout_battery_v0(
    *,
    target_manifest: Path,
    coupling_dir: Path,
    output_dir: Path,
    selector_name: str,
    control_names: Sequence[str],
    target_ids: Sequence[str] = (),
) -> tuple[Path, Path, Path, Path, Path, Path]:
    if selector_name not in (DEFAULT_SELECTOR, FRONTIER_SELECTOR):
        raise ValueError(f"unsupported frozen selector: {selector_name}")

    manifest = _load_json(target_manifest)
    targets = manifest.get("targets", [])
    if not isinstance(targets, list) or not targets:
        raise ValueError("target manifest must include a non-empty targets list")
    requested_target_ids = {target_id.strip() for target_id in target_ids if target_id.strip()}
    if requested_target_ids:
        targets = [
            target
            for target in targets
            if isinstance(target, Mapping)
            and str(target.get("target_id", "")).strip() in requested_target_ids
        ]
        if not targets:
            raise ValueError(
                "target-id filter did not match any manifest targets: "
                + ", ".join(sorted(requested_target_ids))
            )
    controls_to_run = _control_names(tuple(control_names))
    target_manifest_sha256 = _sha256_file(target_manifest)
    current_commit = _current_commit_hash()
    selector_frozen_after_manifest = bool(
        manifest.get("selector_frozen_after_manifest", False)
    )

    target_rows: list[dict[str, object]] = []
    control_rows: list[dict[str, object]] = []
    selected_event_rows: list[dict[str, object]] = []
    failure_rows: list[dict[str, object]] = []
    coordinate_taint_violation_count = 0
    native_taint_violation_count = 0
    oracle_taint_violation_count = 0

    for raw_target in targets:
        if not isinstance(raw_target, Mapping):
            failure_rows.append({"target_id": "", "failure_kind": "invalid_target"})
            continue
        target_id = str(raw_target.get("target_id", "")).strip()
        if not target_id:
            failure_rows.append({"target_id": "", "failure_kind": "missing_target_id"})
            continue
        benchmark_file = _resolve_path(
            raw_target.get("benchmark_file", ""),
            base_dir=REPO_ROOT,
        )
        external_coupling_file = _target_coupling_file(
            raw_target,
            coupling_dir=coupling_dir,
        )
        if not external_coupling_file.is_absolute():
            external_coupling_file = (REPO_ROOT / external_coupling_file).resolve()

        try:
            rows = load_real_coordinate_visual_rows(benchmark_file)
            import_result = import_external_coupling_dataset(
                rows=rows,
                external_coupling_file=external_coupling_file,
            )
            if import_result.dataset.coordinate_truth_tainted:
                coordinate_taint_violation_count += 1
            if import_result.dataset.native_truth_tainted:
                native_taint_violation_count += 1
            if import_result.dataset.oracle_constraint_control:
                oracle_taint_violation_count += 1
            if not import_result.dataset.external_evolutionary_couplings_used:
                failure_rows.append(
                    {
                        "target_id": target_id,
                        "failure_kind": "external_couplings_not_used",
                    }
                )
                continue
            physical_contexts = _build_multiscale_physical_contexts(rows)
            real_precision_run, real_saturation_run = _run_real_frontier(
                rows=rows,
                dataset=import_result.dataset,
                control_kind="external_real_blind_holdout",
                physical_contexts=physical_contexts,
            )
            controls = generate_external_coupling_negative_controls(
                rows=rows,
                dataset=import_result.dataset,
            )
            benchmark_file_sha256 = _sha256_file(benchmark_file)
            coupling_file_sha256 = _sha256_file(external_coupling_file)
            source_accessions = tuple(row.source_accession for row in rows)
            real_exact_metrics = _exact_contact_metrics(
                rows=rows,
                dataset=import_result.dataset,
                run=real_saturation_run,
            )
            real_precision_exact_metrics = _exact_contact_metrics(
                rows=rows,
                dataset=import_result.dataset,
                run=real_precision_run,
            )
            target_row = _metric_row(
                target_id=target_id,
                run=real_saturation_run,
                exact_metrics=real_exact_metrics,
                benchmark_file=benchmark_file,
                external_coupling_file=external_coupling_file,
                benchmark_file_sha256=benchmark_file_sha256,
                external_coupling_file_sha256=coupling_file_sha256,
                source_accessions=source_accessions,
            )
            target_row["frontier_selector_name"] = FRONTIER_SELECTOR
            target_row["precision_boundary_selected_event_count"] = (
                real_precision_run.metric.selected_event_count
            )
            target_row["precision_boundary_false_nucleus_rate"] = (
                real_precision_run.metric.false_nucleus_rate
            )
            target_row["precision_boundary_cluster_precision"] = (
                real_precision_run.metric.contact_cluster_precision
            )
            target_row["precision_boundary_long_range_contact_recall"] = (
                real_precision_run.metric.long_range_contact_recall
            )
            target_row["precision_boundary_exact_contact_precision_top_L"] = (
                real_precision_exact_metrics.exact_contact_precision_top_L
            )
            target_row["precision_boundary_region_width_penalty"] = (
                real_precision_exact_metrics.region_width_penalty
            )
            row_control_rows: list[dict[str, object]] = []
            control_frontier_win_flags: list[int] = []
            row_controls_to_run = tuple(
                control_name
                for control_name in controls_to_run
                if not (
                    control_name == "external_cross_row_swapped"
                    and len(rows) < 2
                )
            )
            for control_name in row_controls_to_run:
                control = controls[control_name]
                control_precision_run, control_saturation_run = _run_control_frontier(
                    rows=rows,
                    dataset=control.dataset,
                    control_name=control_name,
                    control_kind=control.control_kind,
                    physical_contexts=physical_contexts,
                )
                control_frontier_win_flags.append(
                    int(
                        _frontier_dominates(
                            (real_precision_run, real_saturation_run),
                            (control_precision_run, control_saturation_run),
                        )
                    )
                )
                for boundary_name, control_run in (
                    ("precision_boundary", control_precision_run),
                    ("saturation_boundary", control_saturation_run),
                ):
                    control_exact_metrics = _exact_contact_metrics(
                        rows=rows,
                        dataset=control.dataset,
                        run=control_run,
                    )
                    control_row = _metric_row(
                        target_id=target_id,
                        run=control_run,
                        exact_metrics=control_exact_metrics,
                        benchmark_file=benchmark_file,
                        external_coupling_file=external_coupling_file,
                        benchmark_file_sha256=benchmark_file_sha256,
                        external_coupling_file_sha256=coupling_file_sha256,
                        source_accessions=source_accessions,
                        control_name=control_name,
                    )
                    control_row["frontier_member"] = boundary_name
                    row_control_rows.append(control_row)
                    control_rows.append(control_row)

            target_row["max_control_false_nucleus_rate"] = _max_field(
                row_control_rows,
                "false_nucleus_rate",
            )
            target_row["max_control_long_range_contact_recall"] = _max_field(
                row_control_rows,
                "long_range_contact_recall",
            )
            target_row["max_control_cluster_precision"] = _max_field(
                row_control_rows,
                "cluster_precision",
            )
            target_row["frontier_control_count"] = len(control_frontier_win_flags)
            target_row["frontier_beaten_control_count"] = sum(
                control_frontier_win_flags
            )
            target_row["frontier_control_win"] = int(
                bool(control_frontier_win_flags)
                and all(control_frontier_win_flags)
            )
            target_row["individual_control_win"] = target_row[
                "frontier_control_win"
            ]
            target_rows.append(target_row)
            for real_run in (real_precision_run, real_saturation_run):
                for selected_row in real_run.selected_rows:
                    selected_event_rows.append(
                        {
                            "target_id": target_id,
                            "frontier_selector_name": FRONTIER_SELECTOR,
                            "benchmark_file_sha256": benchmark_file_sha256,
                            "external_coupling_file_sha256": coupling_file_sha256,
                            **selected_row,
                        }
                    )
        except Exception as exc:  # noqa: BLE001
            failure_rows.append(
                {
                    "target_id": target_id,
                    "failure_kind": exc.__class__.__name__,
                    "failure_message": str(exc),
                    "benchmark_file": str(benchmark_file),
                    "external_coupling_file": str(external_coupling_file),
                }
            )

    classification = _battery_classification(
        target_rows=target_rows,
        control_rows=control_rows,
        failure_rows=failure_rows,
        oracle_taint_violation_count=oracle_taint_violation_count,
        target_manifest_frozen=bool(manifest.get("manifest_frozen", False)),
        selector_frozen_after_manifest=selector_frozen_after_manifest,
    )
    hard_gates = {
        "coordinate_truth_used_to_build_constraints": (
            coordinate_taint_violation_count > 0
        ),
        "native_truth_used_before_coupling_selection": (
            native_taint_violation_count > 0
        ),
        "oracle_constraint_control": oracle_taint_violation_count > 0,
        "external_evolutionary_couplings_used": not failure_rows,
        "target_manifest_frozen": bool(manifest.get("manifest_frozen", False)),
        "selector_frozen_after_manifest": selector_frozen_after_manifest,
    }
    hard_gates["all_hard_gates_passed"] = (
        not hard_gates["coordinate_truth_used_to_build_constraints"]
        and not hard_gates["native_truth_used_before_coupling_selection"]
        and not hard_gates["oracle_constraint_control"]
        and bool(hard_gates["external_evolutionary_couplings_used"])
        and bool(hard_gates["target_manifest_frozen"])
        and bool(hard_gates["selector_frozen_after_manifest"])
    )
    report = {
        "report_kind": BATTERY_REPORT_KIND,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "commit_hash": current_commit,
        "selector_name": FRONTIER_SELECTOR,
        "saturation_selector_name": SATURATION_BOUNDARY_SELECTOR,
        "precision_boundary_selector_name": PRECISION_BOUNDARY_SELECTOR,
        "target_manifest": str(target_manifest),
        "target_manifest_sha256": target_manifest_sha256,
        "target_manifest_declared_kind": str(manifest.get("manifest_kind", "")),
        "target_manifest_frozen": bool(manifest.get("manifest_frozen", False)),
        "selector_frozen_after_manifest": selector_frozen_after_manifest,
        "control_names": tuple(controls_to_run),
        "target_id_filter": tuple(sorted(requested_target_ids)),
        "target_count": len(target_rows),
        "failure_count": len(failure_rows),
        "coordinate_taint_violation_count": coordinate_taint_violation_count,
        "native_taint_violation_count": native_taint_violation_count,
        "oracle_taint_violation_count": oracle_taint_violation_count,
        "mean_false_nucleus_rate": _mean_field(target_rows, "false_nucleus_rate"),
        "mean_long_range_contact_recall": _mean_field(
            target_rows,
            "long_range_contact_recall",
        ),
        "mean_cluster_precision": _mean_field(target_rows, "cluster_precision"),
        "mean_exact_contact_precision_top_L": _mean_field(
            target_rows,
            "exact_contact_precision_top_L",
        ),
        "mean_exact_contact_precision_L_over_2": _mean_field(
            target_rows,
            "exact_contact_precision_L_over_2",
        ),
        "mean_exact_long_range_contact_precision": _mean_field(
            target_rows,
            "exact_long_range_contact_precision",
        ),
        "mean_region_width_penalty": _mean_field(
            target_rows,
            "region_width_penalty",
        ),
        "frontier_control_win_rate": _mean_field(
            target_rows,
            "frontier_control_win",
        ),
        "individual_control_win_rate": _mean_field(
            target_rows,
            "individual_control_win",
        ),
        "max_control_false_nucleus_rate": _max_field(
            control_rows,
            "false_nucleus_rate",
        ),
        "max_control_long_range_contact_recall": _max_field(
            control_rows,
            "long_range_contact_recall",
        ),
        "max_control_cluster_precision": _max_field(
            control_rows,
            "cluster_precision",
        ),
        "max_control_exact_contact_precision_top_L": _max_field(
            control_rows,
            "exact_contact_precision_top_L",
        ),
        "max_control_exact_long_range_contact_precision": _max_field(
            control_rows,
            "exact_long_range_contact_precision",
        ),
        "hard_gates": hard_gates,
        "final_classification": classification,
        "folding_problem_solved": False,
        "claim_allowed": classification
        in (
            "blind_batch_scaffold_signal_confirmed",
            "blind_batch_exact_contact_signal_confirmed",
        )
        and bool(hard_gates["all_hard_gates_passed"]),
    }
    certificate = {
        "certificate_kind": BATTERY_CERTIFICATE_KIND,
        "report_kind": BATTERY_REPORT_KIND,
        "commit_hash": current_commit,
        "selector_name": FRONTIER_SELECTOR,
        "saturation_selector_name": SATURATION_BOUNDARY_SELECTOR,
        "precision_boundary_selector_name": PRECISION_BOUNDARY_SELECTOR,
        "target_manifest_sha256": target_manifest_sha256,
        "final_classification": classification,
        "folding_problem_solved": False,
        "claim_allowed": report["claim_allowed"],
        "claim_scope": (
            "long_range_folding_nucleus_region_recovery_only_not_full_atomic_folding"
        ),
        "hard_gates": report["hard_gates"],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = _write_json(output_dir / "blind_external_holdout_report.json", report)
    rows_path = write_csv_rows(
        target_rows,
        output_dir / "blind_external_holdout_rows.csv",
    )
    selected_events_path = write_csv_rows(
        selected_event_rows,
        output_dir / "blind_external_holdout_selected_events.csv",
    )
    controls_path = write_csv_rows(
        control_rows,
        output_dir / "blind_external_holdout_controls.csv",
    )
    failures_path = write_csv_rows(
        failure_rows,
        output_dir / "blind_external_holdout_failures.csv",
    )
    certificate_path = _write_json(
        output_dir / "blind_external_holdout_certificate.json",
        certificate,
    )
    return (
        report_path,
        rows_path,
        selected_events_path,
        controls_path,
        failures_path,
        certificate_path,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run a frozen blind external-coupling holdout battery with hard "
            "taint gates and matched controls."
        )
    )
    parser.add_argument("--target-manifest", default=str(DEFAULT_TARGET_MANIFEST))
    parser.add_argument("--coupling-dir", default=str(DEFAULT_COUPLING_DIR))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--selector", default=DEFAULT_SELECTOR)
    parser.add_argument(
        "--target-id",
        action="append",
        default=[],
        help=(
            "Optional manifest target_id to run. Repeat for several. Omit to "
            "run every target."
        ),
    )
    parser.add_argument(
        "--control-name",
        action="append",
        choices=tuple(EXTERNAL_COUPLING_CONTROL_NAMES) + ("all",),
        default=[],
        help=(
            "Matched control to run. Repeat for several. Omit or pass 'all' "
            "to run every matched control."
        ),
    )
    args = parser.parse_args()
    for output in run_blind_external_holdout_battery_v0(
        target_manifest=Path(args.target_manifest),
        coupling_dir=Path(args.coupling_dir),
        output_dir=Path(args.output_dir),
        selector_name=args.selector,
        control_names=tuple(args.control_name),
        target_ids=tuple(args.target_id),
    ):
        print(output)


if __name__ == "__main__":
    main()
