from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from dataclasses import asdict, dataclass
from math import sqrt
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from pharmacotopology.folding_contact_topology import (
    PREDICTOR_INPUT_BOUNDARY,
    ContactTopologyPrediction,
    predict_contact_topology,
    sha256_sequence,
)
from pharmacotopology.folding_energy_landscape import (
    EnergyLandscapePacket,
    build_energy_landscape,
    render_curve_svg,
)
from pharmacotopology.folding_native_contact_eval import (
    ContactMetricPacket,
    ContactPair,
    contact_map_hash,
    evaluate_contact_prediction,
    normalized_contact_pairs,
)
from pharmacotopology.folding_reference_loader import load_folding_reference_dataset
from pharmacotopology.folding_structure_benchmark import AA3_TO_1
from pharmacotopology.folding_topology import normalize_sequence
from pharmacotopology.folding_visual_trajectory import (
    render_coarse_grain_svg,
    render_contact_map_svg,
    render_contact_overlay_svg,
    render_trajectory_html,
    write_visual_file,
)


REAL_COORDINATE_VISUAL_BENCHMARK_KIND = "real_coordinate_visual_contact_benchmark_v1"
REAL_COORDINATE_VISUAL_SPLIT = "real_coordinate_visual_8"
REAL_COORDINATE_NATIVE_KIND = "rcsb_ca_coordinate_native_contacts_v1"
REAL_COORDINATE_VISUAL_SIGNATURE_KIND = (
    "sequence_only_against_coordinate_derived_native_contacts"
)
REAL_COORDINATE_VISUAL_CERTIFICATE_KIND = (
    "real_coordinate_visual_contact_safety_certificate"
)

CONTACT_CUTOFF_ANGSTROM = 8.0
MIN_SEQUENCE_SEPARATION = 3

ROOT_OUTPUT_NAMES = (
    "real_coordinate_visual_8_report.json",
    "real_coordinate_visual_8_rows.csv",
    "real_coordinate_visual_8_contact_metrics.csv",
    "real_coordinate_visual_8_native_contact_summary.csv",
    "real_coordinate_visual_8_dashboard.html",
    "real_coordinate_visual_8_certificate.json",
)

PER_ROW_VISUAL_NAMES = (
    "native_coordinate_trace.svg",
    "native_contact_map.svg",
    "predicted_contact_map.svg",
    "contact_map_overlay.svg",
    "folding_trajectory.html",
    "energy_curve.svg",
    "contact_closure_curve.svg",
    "coarse_grain_final.svg",
)


@dataclass(frozen=True)
class CoordinatePoint:
    sequence_index: int
    residue_number: int
    insertion_code: str
    x: float
    y: float
    z: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RealCoordinateVisualRow:
    row_id: str
    source_id: str
    source_accession: str
    source_database: str
    reference_structure_source: str
    reference_fold_class: str
    sequence: str
    sequence_sha256: str
    sequence_length: int
    coordinate_source_url: str
    coordinate_points: tuple[CoordinatePoint, ...]
    coordinate_trace_hash: str
    coordinate_residue_count: int
    coordinate_coverage: float
    native_contact_map_hash: str
    truth_axes: dict[str, str]

    def native_contact_pairs(self) -> tuple[ContactPair, ...]:
        return coordinate_native_contact_pairs(self.coordinate_points)

    def to_safe_dict(self) -> dict[str, object]:
        return {
            "row_id": self.row_id,
            "source_id": self.source_id,
            "source_accession": self.source_accession,
            "source_database": self.source_database,
            "reference_structure_source": self.reference_structure_source,
            "reference_fold_class": self.reference_fold_class,
            "sequence_sha256": self.sequence_sha256,
            "sequence_length": self.sequence_length,
            "coordinate_source_url": self.coordinate_source_url,
            "coordinate_trace_hash": self.coordinate_trace_hash,
            "coordinate_residue_count": self.coordinate_residue_count,
            "coordinate_coverage": self.coordinate_coverage,
            "native_contact_map_hash": self.native_contact_map_hash,
            "truth_axes": dict(self.truth_axes),
        }


@dataclass(frozen=True)
class RealCoordinateVisualPacket:
    row: RealCoordinateVisualRow
    native_contact_pairs: tuple[ContactPair, ...]
    prediction: ContactTopologyPrediction
    metrics: ContactMetricPacket
    energy: EnergyLandscapePacket
    visible_partial_success: bool
    failure_cohort: str
    visual_paths: Mapping[str, str]

    def safe_row(self) -> dict[str, object]:
        axes = self.row.truth_axes
        return {
            "row_id": self.row.row_id,
            "source_id": self.row.source_id,
            "source_accession": self.row.source_accession,
            "source_database": self.row.source_database,
            "sequence_hash": self.row.sequence_sha256,
            "sequence_length": self.row.sequence_length,
            "reference_fold_class": self.row.reference_fold_class,
            "truth_secondary_structure_axis": axes.get(
                "secondary_structure_axis",
                "weak_or_unknown",
            ),
            "truth_architecture_axis": axes.get("architecture_axis", "unknown"),
            "truth_order_axis": axes.get("order_axis", "unknown"),
            "truth_environment_axis": axes.get("environment_axis", "unknown"),
            "coordinate_residue_count": self.row.coordinate_residue_count,
            "coordinate_coverage": self.row.coordinate_coverage,
            "coordinate_trace_hash": self.row.coordinate_trace_hash,
            "native_contact_derivation_kind": REAL_COORDINATE_NATIVE_KIND,
            "native_contact_count": self.metrics.native_contact_count,
            "predicted_contact_count": self.metrics.predicted_contact_count,
            "true_positive_contacts": self.metrics.true_positive_contacts,
            "false_positive_contacts": self.metrics.false_positive_contacts,
            "false_negative_contacts": self.metrics.false_negative_contacts,
            "contact_map_f1": self.metrics.contact_map_f1,
            "native_contact_precision": self.metrics.native_contact_precision,
            "native_contact_recall": self.metrics.native_contact_recall,
            "long_range_contact_recall": self.metrics.long_range_contact_recall,
            "short_range_contact_recall": self.metrics.short_range_contact_recall,
            "false_contact_rate": self.metrics.false_contact_rate,
            "visible_partial_success": self.visible_partial_success,
            "failure_cohort": self.failure_cohort,
            "native_contact_map_hash": self.row.native_contact_map_hash,
            "predicted_contact_map_hash": self.prediction.predicted_contact_map_hash,
            "native_truth_used_before_prediction": (
                self.prediction.native_truth_used_before_prediction
            ),
            "coordinate_truth_used_before_prediction": False,
            "raw_sequence_exposed": self.prediction.raw_sequence_exposed,
            "repair_heuristic_applied": False,
            "mechanism_discovery_claim_allowed": False,
            "global_folding_claim_allowed": False,
            "folding_problem_solved": False,
            "visual_dir": f"real_coordinate_visuals/{self.row.row_id}",
            **{
                f"visual_{name.replace('.', '_')}": path
                for name, path in self.visual_paths.items()
            },
        }

    def contact_metric_row(self) -> dict[str, object]:
        return {
            "row_id": self.row.row_id,
            "source_accession": self.row.source_accession,
            "sequence_hash": self.row.sequence_sha256,
            "coordinate_residue_count": self.row.coordinate_residue_count,
            "coordinate_coverage": self.row.coordinate_coverage,
            "native_contact_count": self.metrics.native_contact_count,
            "predicted_contact_count": self.metrics.predicted_contact_count,
            "true_positive_contacts": self.metrics.true_positive_contacts,
            "false_positive_contacts": self.metrics.false_positive_contacts,
            "false_negative_contacts": self.metrics.false_negative_contacts,
            "native_contact_recall": self.metrics.native_contact_recall,
            "native_contact_precision": self.metrics.native_contact_precision,
            "contact_map_f1": self.metrics.contact_map_f1,
            "long_range_contact_recall": self.metrics.long_range_contact_recall,
            "short_range_contact_recall": self.metrics.short_range_contact_recall,
            "false_contact_rate": self.metrics.false_contact_rate,
            "visible_partial_success": self.visible_partial_success,
            "failure_cohort": self.failure_cohort,
        }


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _distance(left: CoordinatePoint, right: CoordinatePoint) -> float:
    return sqrt(
        (left.x - right.x) ** 2
        + (left.y - right.y) ** 2
        + (left.z - right.z) ** 2
    )


def coordinate_native_contact_pairs(
    points: Sequence[CoordinatePoint],
    *,
    contact_cutoff_angstrom: float = CONTACT_CUTOFF_ANGSTROM,
    minimum_sequence_separation: int = MIN_SEQUENCE_SEPARATION,
) -> tuple[ContactPair, ...]:
    pairs: list[tuple[int, int]] = []
    ordered = tuple(sorted(points, key=lambda point: point.sequence_index))
    for left_index, left in enumerate(ordered):
        for right in ordered[left_index + 1 :]:
            if right.sequence_index - left.sequence_index < minimum_sequence_separation:
                continue
            if _distance(left, right) <= contact_cutoff_angstrom:
                pairs.append((left.sequence_index, right.sequence_index))
    return normalized_contact_pairs(pairs)


def coordinate_trace_hash(points: Iterable[CoordinatePoint]) -> str:
    normalized = sorted(points, key=lambda point: point.sequence_index)
    encoded = ";".join(
        (
            f"{point.sequence_index}:{point.residue_number}:"
            f"{point.insertion_code}:{point.x:.3f}:{point.y:.3f}:{point.z:.3f}"
        )
        for point in normalized
    )
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _load_json(path: Path) -> Mapping[str, Any]:
    parsed = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, Mapping):
        raise ValueError(f"{path} must contain a JSON object")
    return parsed


def _coordinate_points_from_raw(raw_points: object) -> tuple[CoordinatePoint, ...]:
    if not isinstance(raw_points, list):
        raise ValueError("coordinate_points must be a list")
    points: list[CoordinatePoint] = []
    for index, raw in enumerate(raw_points, start=1):
        if not isinstance(raw, Mapping):
            raise ValueError(f"coordinate_points[{index}] must be an object")
        points.append(
            CoordinatePoint(
                sequence_index=int(raw["sequence_index"]),
                residue_number=int(raw["residue_number"]),
                insertion_code=str(raw.get("insertion_code", "")),
                x=float(raw["x"]),
                y=float(raw["y"]),
                z=float(raw["z"]),
            )
        )
    return tuple(points)


def load_real_coordinate_visual_rows(path: Path) -> tuple[RealCoordinateVisualRow, ...]:
    data = _load_json(path)
    if data.get("benchmark_kind") != REAL_COORDINATE_VISUAL_BENCHMARK_KIND:
        raise ValueError("unexpected real-coordinate visual benchmark kind")
    raw_rows = data.get("references")
    if not isinstance(raw_rows, list):
        raise ValueError("real-coordinate benchmark must include references")
    rows: list[RealCoordinateVisualRow] = []
    for index, raw in enumerate(raw_rows, start=1):
        if not isinstance(raw, Mapping):
            raise ValueError(f"row[{index}] must be an object")
        if "native_contact_pairs" in raw:
            raise ValueError(f"row[{index}].native_contact_pairs_must_not_be_stored")
        sequence = normalize_sequence(str(raw.get("sequence", "")))
        actual_sha = sha256_sequence(sequence)
        if str(raw.get("sequence_sha256", "")) != actual_sha:
            raise ValueError(f"row[{index}].sequence_sha256_mismatch")
        sequence_length = int(raw.get("sequence_length", 0))
        if sequence_length != len(sequence):
            raise ValueError(f"row[{index}].sequence_length_mismatch")
        points = _coordinate_points_from_raw(raw.get("coordinate_points"))
        if int(raw.get("coordinate_residue_count", 0)) != len(points):
            raise ValueError(f"row[{index}].coordinate_residue_count_mismatch")
        actual_trace_hash = coordinate_trace_hash(points)
        if str(raw.get("coordinate_trace_hash", "")) != actual_trace_hash:
            raise ValueError(f"row[{index}].coordinate_trace_hash_mismatch")
        native_hash = contact_map_hash(coordinate_native_contact_pairs(points))
        if str(raw.get("native_contact_map_hash", "")) != native_hash:
            raise ValueError(f"row[{index}].native_contact_map_hash_mismatch")
        truth_axes = raw.get("truth_axes", {})
        if not isinstance(truth_axes, Mapping):
            raise ValueError(f"row[{index}].truth_axes_must_be_object")
        rows.append(
            RealCoordinateVisualRow(
                row_id=str(raw["row_id"]),
                source_id=str(raw["source_id"]),
                source_accession=str(raw["source_accession"]),
                source_database=str(raw["source_database"]),
                reference_structure_source=str(raw["reference_structure_source"]),
                reference_fold_class=str(raw["reference_fold_class"]),
                sequence=sequence,
                sequence_sha256=actual_sha,
                sequence_length=sequence_length,
                coordinate_source_url=str(raw["coordinate_source_url"]),
                coordinate_points=points,
                coordinate_trace_hash=actual_trace_hash,
                coordinate_residue_count=len(points),
                coordinate_coverage=_rounded(len(points) / max(sequence_length, 1)),
                native_contact_map_hash=native_hash,
                truth_axes={str(key): str(value) for key, value in truth_axes.items()},
            )
        )
    return tuple(rows)


def validate_real_coordinate_visual_lock(
    rows: Sequence[RealCoordinateVisualRow],
) -> dict[str, object]:
    row_ids = [row.row_id for row in rows]
    sequence_hashes = [row.sequence_sha256 for row in rows]
    contact_counts = [len(row.native_contact_pairs()) for row in rows]
    violations: list[str] = []
    if len(rows) != 8:
        violations.append("real_coordinate_visual_row_count_not_8")
    if len(set(row_ids)) != len(row_ids):
        violations.append("duplicate_row_id")
    if len(set(sequence_hashes)) != len(sequence_hashes):
        violations.append("duplicate_sequence_sha256")
    if any(row.coordinate_residue_count <= 0 for row in rows):
        violations.append("empty_coordinate_trace")
    if any(count <= 0 for count in contact_counts):
        violations.append("empty_coordinate_native_contact_map")
    if any(row.coordinate_coverage < 0.80 for row in rows):
        violations.append("coordinate_coverage_below_0_80")
    return {
        "real_coordinate_visual_row_count": len(rows),
        "coordinate_backed_row_count": sum(
            1 for row in rows if row.coordinate_residue_count > 0
        ),
        "unique_sequence_count": len(set(sequence_hashes)),
        "native_contact_map_count": len(contact_counts),
        "min_native_contact_count": min(contact_counts) if contact_counts else 0,
        "max_native_contact_count": max(contact_counts) if contact_counts else 0,
        "real_coordinate_visual_lock_valid": not violations,
        "real_coordinate_visual_lock_violations": tuple(violations),
    }


def _visual_paths(row_id: str) -> dict[str, str]:
    return {
        name: f"real_coordinate_visuals/{row_id}/{name}"
        for name in PER_ROW_VISUAL_NAMES
    }


def _visible_partial_success(metrics: ContactMetricPacket) -> bool:
    if metrics.contact_map_f1 >= 0.08:
        return True
    return (
        metrics.native_contact_recall >= 0.04
        and metrics.native_contact_precision >= 0.10
    )


def _failure_cohort(
    row: RealCoordinateVisualRow,
    metrics: ContactMetricPacket,
    visible_partial_success: bool,
) -> str:
    if visible_partial_success:
        return "coordinate_visible_partial_success"
    axes = row.truth_axes
    if axes.get("environment_axis") == "membrane_like":
        return "coordinate_membrane_gap"
    if axes.get("architecture_axis") == "multidomain_or_segmented":
        return "coordinate_architecture_gap"
    if (
        axes.get("secondary_structure_axis") == "beta_rich"
        and metrics.long_range_contact_recall < 0.10
    ):
        return "coordinate_beta_long_range_gap"
    if metrics.true_positive_contacts == 0:
        return "coordinate_no_native_overlap"
    if metrics.false_contact_rate >= 0.90:
        return "coordinate_false_contact_overprediction"
    return "coordinate_low_recall_gap"


def real_coordinate_visual_packets(
    rows: Sequence[RealCoordinateVisualRow],
) -> list[RealCoordinateVisualPacket]:
    packets: list[RealCoordinateVisualPacket] = []
    for row in rows:
        prediction = predict_contact_topology(row.sequence, row_id=row.row_id)
        native_pairs = row.native_contact_pairs()
        metrics = evaluate_contact_prediction(
            native_pairs=native_pairs,
            predicted_pairs=prediction.predicted_contact_pairs,
        )
        energy = build_energy_landscape(prediction.candidates)
        partial_success = _visible_partial_success(metrics)
        packets.append(
            RealCoordinateVisualPacket(
                row=row,
                native_contact_pairs=native_pairs,
                prediction=prediction,
                metrics=metrics,
                energy=energy,
                visible_partial_success=partial_success,
                failure_cohort=_failure_cohort(row, metrics, partial_success),
                visual_paths=_visual_paths(row.row_id),
            )
        )
    return packets


def safe_real_coordinate_visual_rows(
    packets: Sequence[RealCoordinateVisualPacket],
) -> list[dict[str, object]]:
    return [packet.safe_row() for packet in packets]


def contact_metric_rows(
    packets: Sequence[RealCoordinateVisualPacket],
) -> list[dict[str, object]]:
    return [packet.contact_metric_row() for packet in packets]


def native_contact_summary_rows(
    packets: Sequence[RealCoordinateVisualPacket],
) -> list[dict[str, object]]:
    rows = []
    for packet in packets:
        native = packet.native_contact_pairs
        long_count = sum(1 for left, right in native if right - left >= 24)
        short_count = sum(1 for left, right in native if right - left <= 12)
        rows.append(
            {
                "row_id": packet.row.row_id,
                "source_accession": packet.row.source_accession,
                "coordinate_residue_count": packet.row.coordinate_residue_count,
                "coordinate_coverage": packet.row.coordinate_coverage,
                "native_contact_count": len(native),
                "native_long_range_contact_count": long_count,
                "native_short_range_contact_count": short_count,
                "native_contact_derivation_kind": REAL_COORDINATE_NATIVE_KIND,
                "contact_cutoff_angstrom": CONTACT_CUTOFF_ANGSTROM,
                "minimum_sequence_separation": MIN_SEQUENCE_SEPARATION,
                "native_contact_map_hash": packet.row.native_contact_map_hash,
                "coordinate_trace_hash": packet.row.coordinate_trace_hash,
            }
        )
    return rows


def build_real_coordinate_visual_report(
    *,
    packets: Sequence[RealCoordinateVisualPacket],
    source_benchmark_file: Path,
    lock_validation: Mapping[str, object],
) -> dict[str, object]:
    rows = safe_real_coordinate_visual_rows(packets)
    f1_values = [float(row["contact_map_f1"]) for row in rows]
    precision_values = [float(row["native_contact_precision"]) for row in rows]
    recall_values = [float(row["native_contact_recall"]) for row in rows]
    partial_success_count = sum(
        1 for row in rows if bool(row["visible_partial_success"])
    )
    failure_count = len(rows) - partial_success_count
    native_truth_flags = [
        bool(row["native_truth_used_before_prediction"]) for row in rows
    ]
    coordinate_truth_flags = [
        bool(row["coordinate_truth_used_before_prediction"]) for row in rows
    ]
    raw_sequence_flags = [bool(row["raw_sequence_exposed"]) for row in rows]
    output_artifact_count = len(ROOT_OUTPUT_NAMES) + len(rows) * len(PER_ROW_VISUAL_NAMES)
    failure_counts = Counter(str(row["failure_cohort"]) for row in rows)
    return {
        "benchmark_kind": REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
        "real_coordinate_visual_signature_kind": (
            REAL_COORDINATE_VISUAL_SIGNATURE_KIND
        ),
        "holdout_split": REAL_COORDINATE_VISUAL_SPLIT,
        "source_benchmark_file": str(source_benchmark_file),
        "benchmark_size": len(rows),
        "lock_validation": dict(lock_validation),
        "prediction_input_boundary": PREDICTOR_INPUT_BOUNDARY,
        "truth_scoring_boundary": (
            "coordinates_and_coordinate_native_contacts_used_only_after_sequence_prediction"
        ),
        "native_contact_derivation_kind": REAL_COORDINATE_NATIVE_KIND,
        "native_contact_derivation_rule": (
            "C-alpha distance <= 8.0 angstrom with sequence separation >= 3"
        ),
        "real_coordinate_native_contacts_extracted": True,
        "toy_locked_contact_targets_used": False,
        "coarse_ca_only": True,
        "full_atomic_folding_available": False,
        "coordinate_backed_row_count": len(rows),
        "coordinate_contact_map_count": len(rows),
        "visual_artifacts_generated_for_rows": len(rows),
        "visual_files_per_row": len(PER_ROW_VISUAL_NAMES),
        "visual_artifacts_generated_count": output_artifact_count,
        "contact_map_f1_computed_count": len(f1_values),
        "mean_contact_map_f1": (
            _rounded(sum(f1_values) / len(f1_values)) if f1_values else 0.0
        ),
        "max_contact_map_f1": max(f1_values) if f1_values else 0.0,
        "mean_native_contact_precision": (
            _rounded(sum(precision_values) / len(precision_values))
            if precision_values
            else 0.0
        ),
        "mean_native_contact_recall": (
            _rounded(sum(recall_values) / len(recall_values))
            if recall_values
            else 0.0
        ),
        "visible_partial_success_count": partial_success_count,
        "visible_failure_count": failure_count,
        "failures_visualized": failure_count > 0,
        "failure_cohort_count": len(failure_counts),
        "failure_cohorts": dict(sorted(failure_counts.items())),
        "native_truth_used_before_prediction": any(native_truth_flags),
        "coordinate_truth_used_before_prediction": any(coordinate_truth_flags),
        "raw_sequence_exposed": any(raw_sequence_flags),
        "repair_heuristic_applied": False,
        "contact_repair_overfit_risk_separate": True,
        "blind_visual_contact_benchmark": True,
        "coordinate_native_scoring_claim_allowed": True,
        "mechanism_discovery_claim_allowed": False,
        "mechanism_discovery_claim_created": False,
        "global_folding_claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "claim_allowed": False,
        "boundary_statement": (
            "This benchmark extracts coarse C-alpha native contact maps from "
            "locked RCSB coordinate traces, then scores sequence-only contact "
            "hypotheses after prediction. It upgrades the proof target beyond "
            "toy locked contacts, but it is still not atomistic folding, "
            "mechanism discovery, or a solved folding engine."
        ),
        "rows": rows,
    }


def build_real_coordinate_visual_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": REAL_COORDINATE_VISUAL_CERTIFICATE_KIND,
        "benchmark_kind": report["benchmark_kind"],
        "real_coordinate_visual_signature_kind": report[
            "real_coordinate_visual_signature_kind"
        ],
        "holdout_split": report["holdout_split"],
        "benchmark_size": report["benchmark_size"],
        "real_coordinate_native_contacts_extracted": report[
            "real_coordinate_native_contacts_extracted"
        ],
        "toy_locked_contact_targets_used": report["toy_locked_contact_targets_used"],
        "coarse_ca_only": report["coarse_ca_only"],
        "full_atomic_folding_available": report["full_atomic_folding_available"],
        "blind_visual_contact_benchmark": report["blind_visual_contact_benchmark"],
        "repair_heuristic_applied": report["repair_heuristic_applied"],
        "contact_map_f1_computed_count": report["contact_map_f1_computed_count"],
        "native_truth_used_before_prediction": report[
            "native_truth_used_before_prediction"
        ],
        "coordinate_truth_used_before_prediction": report[
            "coordinate_truth_used_before_prediction"
        ],
        "raw_sequence_exposed": report["raw_sequence_exposed"],
        "mechanism_discovery_claim_allowed": report[
            "mechanism_discovery_claim_allowed"
        ],
        "global_folding_claim_allowed": report["global_folding_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "output_artifacts": tuple(output_names),
    }


def write_real_coordinate_visual_outputs(
    *,
    report: Mapping[str, object],
    packets: Sequence[RealCoordinateVisualPacket],
    report_path: Path,
    rows_path: Path,
    contact_metrics_path: Path,
    native_contact_summary_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
    visuals_root: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(safe_real_coordinate_visual_rows(packets), rows_path)
    _write_csv_rows(contact_metric_rows(packets), contact_metrics_path)
    _write_csv_rows(native_contact_summary_rows(packets), native_contact_summary_path)
    dashboard_path.write_text(
        render_real_coordinate_visual_dashboard(report),
        encoding="utf-8",
    )
    certificate = build_real_coordinate_visual_certificate(
        report,
        output_names=ROOT_OUTPUT_NAMES,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    for packet in packets:
        _write_packet_visuals(packet, visuals_root)
    return (
        report_path,
        rows_path,
        contact_metrics_path,
        native_contact_summary_path,
        dashboard_path,
        certificate_path,
    )


def run_real_coordinate_visual_benchmark(
    *,
    benchmark_file: Path,
    report_path: Path,
    rows_path: Path,
    contact_metrics_path: Path,
    native_contact_summary_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
    visuals_root: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    lock_validation = validate_real_coordinate_visual_lock(rows)
    packets = real_coordinate_visual_packets(rows)
    report = build_real_coordinate_visual_report(
        packets=packets,
        source_benchmark_file=benchmark_file,
        lock_validation=lock_validation,
    )
    return write_real_coordinate_visual_outputs(
        report=report,
        packets=packets,
        report_path=report_path,
        rows_path=rows_path,
        contact_metrics_path=contact_metrics_path,
        native_contact_summary_path=native_contact_summary_path,
        dashboard_path=dashboard_path,
        certificate_path=certificate_path,
        visuals_root=visuals_root,
    )


def parse_pdb_ca_coordinate_points(
    pdb_text: str,
    *,
    chain_id: str,
) -> tuple[CoordinatePoint, ...]:
    points: list[CoordinatePoint] = []
    seen: set[tuple[str, int, str]] = set()
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL") and points:
            break
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        alternate = line[16:17].strip()
        chain = line[21:22].strip()
        if atom_name != "CA" or chain != chain_id or alternate not in ("", "A"):
            continue
        residue_name = line[17:20].strip()
        if residue_name not in AA3_TO_1:
            continue
        try:
            residue_number = int(line[22:26])
            insertion_code = line[26:27].strip()
            key = (chain, residue_number, insertion_code)
            if key in seen:
                continue
            seen.add(key)
            points.append(
                CoordinatePoint(
                    sequence_index=len(points) + 1,
                    residue_number=residue_number,
                    insertion_code=insertion_code,
                    x=round(float(line[30:38]), 3),
                    y=round(float(line[38:46]), 3),
                    z=round(float(line[46:54]), 3),
                )
            )
        except ValueError:
            continue
    return tuple(points)


def _axis_for_reference_class(reference_fold_class: str) -> dict[str, str]:
    if reference_fold_class == "alpha_rich":
        secondary = "alpha_rich"
        architecture = "compact_single_domain"
    elif reference_fold_class == "beta_rich":
        secondary = "beta_rich"
        architecture = "compact_single_domain"
    elif reference_fold_class == "alpha_beta_mixed":
        secondary = "alpha_beta_mixed"
        architecture = "compact_single_domain"
    elif reference_fold_class == "multidomain_boundary":
        secondary = "alpha_beta_mixed"
        architecture = "multidomain_or_segmented"
    else:
        secondary = "weak_or_unknown"
        architecture = "unknown"
    return {
        "secondary_structure_axis": secondary,
        "architecture_axis": architecture,
        "order_axis": "ordered",
        "environment_axis": "soluble_like",
    }


def build_real_coordinate_visual_lock_payload(
    *,
    source_benchmark_file: Path,
    pdb_dir: Path,
) -> dict[str, object]:
    dataset = load_folding_reference_dataset(source_benchmark_file, require_external=True)
    rows: list[dict[str, object]] = []
    for reference in dataset.references:
        accession = reference.source_accession
        if ":" not in accession or not accession[:4].isalnum():
            continue
        pdb_id, chain_id = accession.split(":", 1)
        pdb_id = pdb_id.upper()
        pdb_path = pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            continue
        points = parse_pdb_ca_coordinate_points(
            pdb_path.read_text(encoding="utf-8"),
            chain_id=chain_id,
        )
        if not points:
            continue
        sequence = normalize_sequence(reference.sequence)
        native_pairs = coordinate_native_contact_pairs(points)
        index = len(rows) + 1
        rows.append(
            {
                "row_id": f"coord_{index:03d}_{reference.protein_id}",
                "source_id": reference.protein_id,
                "source_accession": accession,
                "source_database": reference.source_database,
                "reference_structure_source": reference.reference_structure_source,
                "reference_fold_class": reference.reference_fold_class,
                "sequence": sequence,
                "sequence_sha256": sha256_sequence(sequence),
                "sequence_length": len(sequence),
                "coordinate_source_url": (
                    f"https://files.rcsb.org/download/{pdb_id}.pdb"
                ),
                "coordinate_points": [point.to_dict() for point in points],
                "coordinate_trace_hash": coordinate_trace_hash(points),
                "coordinate_residue_count": len(points),
                "coordinate_coverage": _rounded(len(points) / max(len(sequence), 1)),
                "native_contact_map_hash": contact_map_hash(native_pairs),
                "truth_axes": _axis_for_reference_class(reference.reference_fold_class),
            }
        )
        if len(rows) == 8:
            break
    return {
        "benchmark_kind": REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
        "holdout_split": REAL_COORDINATE_VISUAL_SPLIT,
        "source_benchmark_file": str(source_benchmark_file),
        "source_coordinate_database": "RCSB_PDB",
        "coordinate_fixture_kind": "locked_c_alpha_coordinate_trace_v1",
        "native_contact_derivation_kind": REAL_COORDINATE_NATIVE_KIND,
        "native_contact_derivation_rule": (
            "C-alpha distance <= 8.0 angstrom with sequence separation >= 3"
        ),
        "contact_cutoff_angstrom": CONTACT_CUTOFF_ANGSTROM,
        "minimum_sequence_separation": MIN_SEQUENCE_SEPARATION,
        "benchmark_size": len(rows),
        "locked_after_generation": True,
        "blind_prediction_before_coordinate_scoring": True,
        "toy_locked_contact_targets_used": False,
        "coarse_ca_only": True,
        "full_atomic_folding_available": False,
        "mechanism_discovery_claim_allowed": False,
        "global_folding_claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "references": rows,
    }


def write_real_coordinate_visual_lock(
    *,
    source_benchmark_file: Path,
    pdb_dir: Path,
    output_path: Path,
) -> Path:
    payload = build_real_coordinate_visual_lock_payload(
        source_benchmark_file=source_benchmark_file,
        pdb_dir=pdb_dir,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path


def _write_packet_visuals(
    packet: RealCoordinateVisualPacket,
    visuals_root: Path,
) -> None:
    row_dir = visuals_root / packet.row.row_id
    write_visual_file(
        row_dir / "native_coordinate_trace.svg",
        render_coordinate_trace_svg(packet.row),
    )
    write_visual_file(
        row_dir / "native_contact_map.svg",
        render_contact_map_svg(
            row_id=packet.row.row_id,
            sequence_length=packet.row.sequence_length,
            contact_pairs=packet.native_contact_pairs,
            title="Coordinate-derived C-alpha native contacts",
            color="#294f9b",
        ),
    )
    write_visual_file(
        row_dir / "predicted_contact_map.svg",
        render_contact_map_svg(
            row_id=packet.row.row_id,
            sequence_length=packet.row.sequence_length,
            contact_pairs=packet.prediction.predicted_contact_pairs,
            title="Sequence-only predicted contact candidates",
            color="#c44b3a",
        ),
    )
    write_visual_file(
        row_dir / "contact_map_overlay.svg",
        render_contact_overlay_svg(
            row_id=packet.row.row_id,
            sequence_length=packet.row.sequence_length,
            native_pairs=packet.native_contact_pairs,
            predicted_pairs=packet.prediction.predicted_contact_pairs,
        ),
    )
    write_visual_file(
        row_dir / "folding_trajectory.html",
        render_trajectory_html(
            row_id=packet.row.row_id,
            sequence_length=packet.row.sequence_length,
            candidates=packet.prediction.candidates,
            energy=packet.energy,
        ),
    )
    write_visual_file(
        row_dir / "energy_curve.svg",
        render_curve_svg(
            packet.energy.energy_values,
            title="Coarse energy descent curve",
            y_label="relative energy",
        ),
    )
    write_visual_file(
        row_dir / "contact_closure_curve.svg",
        render_curve_svg(
            tuple(float(value) for value in packet.energy.contact_counts),
            title="Contact closure curve",
            y_label="active contacts",
        ),
    )
    write_visual_file(
        row_dir / "coarse_grain_final.svg",
        render_coarse_grain_svg(
            row_id=packet.row.row_id,
            sequence_length=packet.row.sequence_length,
            contacts=packet.prediction.predicted_contact_pairs,
        ),
    )


def render_coordinate_trace_svg(
    row: RealCoordinateVisualRow,
    *,
    width: int = 760,
    height: int = 260,
) -> str:
    points = row.coordinate_points
    if not points:
        return ""
    min_x = min(point.x for point in points)
    max_x = max(point.x for point in points)
    min_y = min(point.y for point in points)
    max_y = max(point.y for point in points)
    span_x = max(max_x - min_x, 1.0)
    span_y = max(max_y - min_y, 1.0)
    margin = 36

    def scaled(point: CoordinatePoint) -> tuple[float, float]:
        x = margin + (point.x - min_x) / span_x * (width - margin * 2)
        y = margin + (point.y - min_y) / span_y * (height - margin * 2)
        return x, y

    path_points = " ".join(
        f"{x:.2f},{y:.2f}" for x, y in (scaled(point) for point in points)
    )
    native_contacts = coordinate_native_contact_pairs(points)
    contact_lines = []
    by_index = {point.sequence_index: point for point in points}
    for left, right in native_contacts[:42]:
        if left not in by_index or right not in by_index:
            continue
        lx, ly = scaled(by_index[left])
        rx, ry = scaled(by_index[right])
        contact_lines.append(
            f'<line x1="{lx:.2f}" y1="{ly:.2f}" x2="{rx:.2f}" y2="{ry:.2f}" '
            'stroke="#326e8f" stroke-width="1" stroke-opacity="0.22"/>'
        )
    bead_step = max(1, len(points) // 45)
    beads = []
    for point in points[::bead_step]:
        x, y = scaled(point)
        beads.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3" fill="#22302d"/>')
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="native coordinate trace">
  <rect width="{width}" height="{height}" fill="#f8faf7"/>
  <text x="34" y="24" font-family="Arial, sans-serif" font-size="15" fill="#22302d">Native C-alpha Coordinate Trace</text>
  {''.join(contact_lines)}
  <polyline points="{path_points}" fill="none" stroke="#7d8f86" stroke-width="2.6"/>
  {''.join(beads)}
  <text x="34" y="{height - 12}" font-family="Arial, sans-serif" font-size="11" fill="#58635e">row: {_escape(row.row_id)} | projected coordinates, not full atomic folding</text>
</svg>
"""


def _write_csv_rows(
    rows: Sequence[Mapping[str, object]],
    path: Path,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        if not rows:
            return path
        fieldnames = list(rows[0])
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    return path


def _escape(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _metric_cards(report: Mapping[str, object]) -> str:
    labels = (
        "benchmark_size",
        "coordinate_backed_row_count",
        "real_coordinate_native_contacts_extracted",
        "toy_locked_contact_targets_used",
        "mean_contact_map_f1",
        "visible_partial_success_count",
        "visible_failure_count",
        "repair_heuristic_applied",
        "mechanism_discovery_claim_allowed",
        "folding_problem_solved",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _rule_cards() -> str:
    rules = (
        (
            "Real Coordinate Native Contacts",
            "Native maps are derived from locked RCSB C-alpha coordinate traces.",
        ),
        (
            "Prediction Is Blind To Coordinates",
            "Coordinates and native contacts are used only after sequence-only prediction.",
        ),
        (
            "C-alpha Coarse, Not Full Atomic Folding",
            "The benchmark scores contact topology, not atomistic pathways or dynamics.",
        ),
        (
            "No Repair Heuristic Applied",
            "This surface is a baseline proof-target upgrade, not another repair pass.",
        ),
        (
            "Global Folding Claim Remains Locked",
            "The benchmark can expose failures but cannot claim folding is solved.",
        ),
    )
    return "".join(
        "<div class=\"rule\">"
        f"<h3>{_escape(title)}</h3><p>{_escape(body)}</p>"
        "</div>"
        for title, body in rules
    )


def _dashboard_preview_rows(report: Mapping[str, object]) -> str:
    rows = report.get("rows", [])
    if not isinstance(rows, Sequence):
        return ""
    body = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        overlay = _escape(row["visual_contact_map_overlay_svg"])
        trace = _escape(row["visual_native_coordinate_trace_svg"])
        body.append(
            "<tr>"
            f"<td>{_escape(row['row_id'])}</td>"
            f"<td>{_escape(row['source_accession'])}</td>"
            f"<td>{_escape(row['reference_fold_class'])}</td>"
            f"<td>{_escape(row['native_contact_count'])}</td>"
            f"<td>{_escape(row['contact_map_f1'])}</td>"
            f"<td>{_escape(row['failure_cohort'])}</td>"
            f"<td><a href=\"{overlay}\">overlay</a> | "
            f"<a href=\"{trace}\">trace</a></td>"
            "</tr>"
        )
    return (
        "<section><h2>Coordinate Visual Rows</h2>"
        "<table><thead><tr>"
        "<th>row</th><th>source</th><th>class</th><th>native contacts</th>"
        "<th>contact_map_f1</th><th>cohort</th><th>visuals</th>"
        "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


def _visual_grid(report: Mapping[str, object]) -> str:
    rows = report.get("rows", [])
    if not isinstance(rows, Sequence):
        return ""
    cards = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        cards.append(
            "<article class=\"visual-card\">"
            f"<h3>{_escape(row['row_id'])}</h3>"
            f"<img src=\"{_escape(row['visual_contact_map_overlay_svg'])}\" "
            f"alt=\"coordinate contact overlay for {_escape(row['row_id'])}\">"
            f"<p>F1: {_escape(row['contact_map_f1'])} | "
            f"{_escape(row['failure_cohort'])}</p>"
            "</article>"
        )
    return (
        "<section><h2>Coordinate Contact Overlays</h2>"
        "<div class=\"visual-grid\">"
        + "".join(cards)
        + "</div></section>"
    )


def render_real_coordinate_visual_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Real Coordinate Visual Contact Benchmark</title>
  <style>
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f5f6f1;
      color: #202623;
    }}
    header {{
      padding: 34px;
      background: #24302c;
      color: #f6f7f2;
    }}
    main {{
      max-width: 1220px;
      margin: 0 auto;
      padding: 24px;
    }}
    h1, h2, h3 {{
      margin: 0 0 12px;
      letter-spacing: 0;
    }}
    section {{
      margin: 24px 0;
    }}
    .metrics, .rules, .visual-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
    }}
    .metric, .rule, .visual-card {{
      background: #ffffff;
      border: 1px solid #d4ddd6;
      border-radius: 6px;
      padding: 14px;
    }}
    .metric span {{
      display: block;
      color: #58635e;
      font-size: 12px;
      overflow-wrap: anywhere;
    }}
    .metric strong {{
      display: block;
      margin-top: 8px;
      font-size: 20px;
      overflow-wrap: anywhere;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: #ffffff;
      border: 1px solid #d4ddd6;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid #e3e8e3;
      text-align: left;
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }}
    .visual-card img {{
      width: 100%;
      height: auto;
      border: 1px solid #d4ddd6;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Real Coordinate Visual Contact Benchmark</h1>
    <p>Sequence-only contact hypotheses are scored against C-alpha native contact maps extracted from locked RCSB coordinate traces.</p>
  </header>
  <main>
    <section class="metrics">{_metric_cards(report)}</section>
    <section><h2>Boundary Rules</h2><div class="rules">{_rule_cards()}</div></section>
    {_dashboard_preview_rows(report)}
    {_visual_grid(report)}
    <section>
      <h2>Claim Boundary</h2>
      <p>{_escape(report.get("boundary_statement", ""))}</p>
    </section>
  </main>
</body>
</html>
"""
