from __future__ import annotations

import csv
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_contact_gap_analysis import (
    CONTACT_GAP_ANALYSIS_SIGNATURE_KIND,
    ContactGapPacket,
    analyze_contact_gap,
    gap_rows_from_packets,
)
from pharmacotopology.folding_contact_topology import (
    VISUAL_MECHANISM_BENCHMARK_KIND,
    VisualMechanismRow,
    load_visual_mechanism_rows,
    predict_contact_topology,
    validate_visual_mechanism_lock,
)
from pharmacotopology.folding_energy_landscape import (
    EnergyLandscapePacket,
    build_energy_landscape,
    render_curve_svg,
)
from pharmacotopology.folding_long_range_contact_repair import (
    LONG_RANGE_CONTACT_REPAIR_SIGNATURE_KIND,
    REPAIR_INPUT_BOUNDARY,
    ContactRepairPacket,
    repair_contact_topology,
)
from pharmacotopology.folding_native_contact_eval import (
    ContactMetricPacket,
    evaluate_contact_prediction,
)
from pharmacotopology.folding_visual_trajectory import (
    render_coarse_grain_svg,
    render_contact_map_svg,
    render_contact_overlay_svg,
    render_trajectory_html,
    write_visual_file,
)


CONTACT_TOPOLOGY_REPAIR_BENCHMARK_KIND = "contact_topology_repair_benchmark_v1"
CONTACT_TOPOLOGY_REPAIR_SIGNATURE_KIND = "native_gap_guided_repair_scoring_v1"
CONTACT_TOPOLOGY_REPAIR_CERTIFICATE_KIND = (
    "contact_topology_repair_safety_certificate"
)

ROOT_OUTPUT_NAMES = (
    "contact_topology_repair_12_report.json",
    "contact_topology_repair_12_rows.csv",
    "contact_topology_repair_12_gap_analysis.csv",
    "contact_topology_repair_12_failure_cohorts.csv",
    "contact_topology_repair_12_dashboard.html",
    "contact_topology_repair_12_certificate.json",
)

PER_ROW_REPAIR_VISUAL_NAMES = (
    "native_contact_map.svg",
    "repaired_contact_map.svg",
    "repair_contact_map_overlay.svg",
    "repair_folding_trajectory.html",
    "repair_energy_curve.svg",
    "repair_contact_closure_curve.svg",
    "repair_coarse_grain_final.svg",
)


@dataclass(frozen=True)
class ContactTopologyRepairBenchmarkPacket:
    row: VisualMechanismRow
    repair: ContactRepairPacket
    baseline_metrics: ContactMetricPacket
    repaired_metrics: ContactMetricPacket
    baseline_energy: EnergyLandscapePacket
    repaired_energy: EnergyLandscapePacket
    baseline_gap: ContactGapPacket
    repaired_gap: ContactGapPacket
    baseline_visible_partial_success: bool
    repaired_visible_partial_success: bool
    visual_paths: Mapping[str, str]

    def safe_row(self) -> dict[str, object]:
        axes = self.row.truth_axes
        beta_baseline = self.baseline_gap.beta_pairing_contact_recall
        beta_repaired = self.repaired_gap.beta_pairing_contact_recall
        return {
            "row_id": self.row.row_id,
            "source_id": self.row.source_id,
            "sequence_hash": self.row.sequence_sha256,
            "sequence_length": self.row.length,
            "truth_secondary_structure_axis": axes.get(
                "secondary_structure_axis",
                "weak_or_unknown",
            ),
            "truth_architecture_axis": axes.get("architecture_axis", "unknown"),
            "truth_order_axis": axes.get("order_axis", "unknown"),
            "truth_environment_axis": axes.get("environment_axis", "unknown"),
            "baseline_predicted_contact_count": (
                self.baseline_metrics.predicted_contact_count
            ),
            "repaired_predicted_contact_count": (
                self.repaired_metrics.predicted_contact_count
            ),
            "baseline_contact_map_f1": self.baseline_metrics.contact_map_f1,
            "repaired_contact_map_f1": self.repaired_metrics.contact_map_f1,
            "contact_map_f1_delta": round(
                self.repaired_metrics.contact_map_f1
                - self.baseline_metrics.contact_map_f1,
                6,
            ),
            "baseline_long_range_contact_recall": (
                self.baseline_metrics.long_range_contact_recall
            ),
            "repaired_long_range_contact_recall": (
                self.repaired_metrics.long_range_contact_recall
            ),
            "long_range_contact_recall_delta": round(
                self.repaired_metrics.long_range_contact_recall
                - self.baseline_metrics.long_range_contact_recall,
                6,
            ),
            "baseline_short_range_contact_recall": (
                self.baseline_metrics.short_range_contact_recall
            ),
            "repaired_short_range_contact_recall": (
                self.repaired_metrics.short_range_contact_recall
            ),
            "baseline_beta_pairing_contact_recall": (
                "" if beta_baseline is None else beta_baseline
            ),
            "repaired_beta_pairing_contact_recall": (
                "" if beta_repaired is None else beta_repaired
            ),
            "baseline_native_cluster_miss_count": (
                self.baseline_gap.native_cluster_miss_count
            ),
            "repaired_native_cluster_miss_count": (
                self.repaired_gap.native_cluster_miss_count
            ),
            "baseline_false_contact_cluster_count": (
                self.baseline_gap.false_contact_cluster_count
            ),
            "repaired_false_contact_cluster_count": (
                self.repaired_gap.false_contact_cluster_count
            ),
            "baseline_premature_compaction": (
                self.baseline_gap.premature_compaction
            ),
            "repaired_premature_compaction": (
                self.repaired_gap.premature_compaction
            ),
            "baseline_failure_mechanism": self.baseline_gap.failure_mechanism,
            "repaired_failure_mechanism": self.repaired_gap.failure_mechanism,
            "baseline_closure_timing_label": (
                self.baseline_gap.closure_timing_label
            ),
            "repaired_closure_timing_label": (
                self.repaired_gap.closure_timing_label
            ),
            "baseline_visible_partial_success": (
                self.baseline_visible_partial_success
            ),
            "repaired_visible_partial_success": (
                self.repaired_visible_partial_success
            ),
            "compact_anchor_candidate_count": (
                self.repair.compact_anchor_candidate_count
            ),
            "beta_pairing_candidate_count": (
                self.repair.beta_pairing_candidate_count
            ),
            "local_overclosure_trimmed_count": (
                self.repair.local_overclosure_trimmed_count
            ),
            "repair_candidate_count": self.repair.repair_candidate_count,
            "repair_notes": ";".join(self.repair.repair_notes),
            "native_truth_used_before_prediction": (
                self.repair.baseline_prediction.native_truth_used_before_prediction
            ),
            "native_truth_used_before_repair": (
                self.repair.native_truth_used_before_repair
            ),
            "raw_sequence_exposed": self.repair.raw_sequence_exposed,
            "global_folding_claim_allowed": False,
            "folding_problem_solved": False,
            "visual_dir": f"contact_repair_visuals/{self.row.row_id}",
            **{
                f"visual_{name.replace('.', '_')}": path
                for name, path in self.visual_paths.items()
            },
        }


def contact_topology_repair_packets(
    rows: Sequence[VisualMechanismRow],
) -> list[ContactTopologyRepairBenchmarkPacket]:
    packets: list[ContactTopologyRepairBenchmarkPacket] = []
    for row in rows:
        baseline_prediction = predict_contact_topology(row.sequence, row_id=row.row_id)
        repair = repair_contact_topology(
            row.sequence,
            baseline_prediction=baseline_prediction,
        )
        baseline_metrics = evaluate_contact_prediction(
            native_pairs=row.native_contact_pairs,
            predicted_pairs=baseline_prediction.predicted_contact_pairs,
        )
        repaired_metrics = evaluate_contact_prediction(
            native_pairs=row.native_contact_pairs,
            predicted_pairs=repair.repaired_prediction.predicted_contact_pairs,
        )
        baseline_energy = build_energy_landscape(baseline_prediction.candidates)
        repaired_energy = build_energy_landscape(repair.repaired_prediction.candidates)
        baseline_gap = analyze_contact_gap(
            row=row,
            prediction=baseline_prediction,
            metrics=baseline_metrics,
            energy=baseline_energy,
            visible_partial_success=False,
        )
        baseline_visible = visible_partial_success(
            metrics=baseline_metrics,
            gap=baseline_gap,
        )
        baseline_gap = analyze_contact_gap(
            row=row,
            prediction=baseline_prediction,
            metrics=baseline_metrics,
            energy=baseline_energy,
            visible_partial_success=baseline_visible,
        )
        repaired_gap = analyze_contact_gap(
            row=row,
            prediction=repair.repaired_prediction,
            metrics=repaired_metrics,
            energy=repaired_energy,
            visible_partial_success=False,
        )
        repaired_visible = visible_partial_success(
            metrics=repaired_metrics,
            gap=repaired_gap,
        )
        repaired_gap = analyze_contact_gap(
            row=row,
            prediction=repair.repaired_prediction,
            metrics=repaired_metrics,
            energy=repaired_energy,
            visible_partial_success=repaired_visible,
        )
        packets.append(
            ContactTopologyRepairBenchmarkPacket(
                row=row,
                repair=repair,
                baseline_metrics=baseline_metrics,
                repaired_metrics=repaired_metrics,
                baseline_energy=baseline_energy,
                repaired_energy=repaired_energy,
                baseline_gap=baseline_gap,
                repaired_gap=repaired_gap,
                baseline_visible_partial_success=baseline_visible,
                repaired_visible_partial_success=repaired_visible,
                visual_paths=_repair_visual_paths(row.row_id),
            )
        )
    return packets


def visible_partial_success(
    *,
    metrics: ContactMetricPacket,
    gap: ContactGapPacket,
) -> bool:
    if metrics.contact_map_f1 >= 0.12:
        return True
    if (
        metrics.short_range_contact_recall >= 0.20
        and metrics.native_contact_precision >= 0.08
    ):
        return True
    return (
        gap.beta_pairing_contact_recall is not None
        and gap.beta_pairing_contact_recall >= 0.25
        and metrics.native_contact_recall >= 0.20
        and metrics.native_contact_precision >= 0.07
    )


def safe_repair_rows(
    packets: Sequence[ContactTopologyRepairBenchmarkPacket],
) -> list[dict[str, object]]:
    return [packet.safe_row() for packet in packets]


def repair_gap_rows(
    packets: Sequence[ContactTopologyRepairBenchmarkPacket],
) -> list[dict[str, object]]:
    rows = []
    for packet in packets:
        baseline = packet.baseline_gap.to_safe_dict()
        repaired = packet.repaired_gap.to_safe_dict()
        row = {
            "row_id": packet.row.row_id,
            "baseline_failure_mechanism": baseline["failure_mechanism"],
            "repaired_failure_mechanism": repaired["failure_mechanism"],
            "baseline_native_cluster_miss_count": baseline[
                "native_cluster_miss_count"
            ],
            "repaired_native_cluster_miss_count": repaired[
                "native_cluster_miss_count"
            ],
            "baseline_false_contact_cluster_count": baseline[
                "false_contact_cluster_count"
            ],
            "repaired_false_contact_cluster_count": repaired[
                "false_contact_cluster_count"
            ],
            "baseline_missed_native_contact_clusters": baseline[
                "missed_native_contact_clusters"
            ],
            "repaired_missed_native_contact_clusters": repaired[
                "missed_native_contact_clusters"
            ],
            "baseline_false_predicted_contact_clusters": baseline[
                "false_predicted_contact_clusters"
            ],
            "repaired_false_predicted_contact_clusters": repaired[
                "false_predicted_contact_clusters"
            ],
            "baseline_closure_timing_label": baseline["closure_timing_label"],
            "repaired_closure_timing_label": repaired["closure_timing_label"],
            "baseline_premature_compaction": baseline["premature_compaction"],
            "repaired_premature_compaction": repaired["premature_compaction"],
            "baseline_beta_pairing_contact_recall": baseline[
                "beta_pairing_contact_recall"
            ],
            "repaired_beta_pairing_contact_recall": repaired[
                "beta_pairing_contact_recall"
            ],
        }
        rows.append(row)
    return rows


def repair_failure_cohort_rows(
    packets: Sequence[ContactTopologyRepairBenchmarkPacket],
) -> list[dict[str, object]]:
    counts = Counter(packet.repaired_gap.failure_mechanism for packet in packets)
    baseline_counts = Counter(packet.baseline_gap.failure_mechanism for packet in packets)
    cohorts = sorted(set(counts) | set(baseline_counts))
    return [
        {
            "failure_mechanism": cohort,
            "baseline_row_count": baseline_counts.get(cohort, 0),
            "repaired_row_count": counts.get(cohort, 0),
            "repaired_row_ids": ";".join(
                packet.row.row_id
                for packet in packets
                if packet.repaired_gap.failure_mechanism == cohort
            ),
            "failures_visualized": cohort != "visible_partial_success",
        }
        for cohort in cohorts
    ]


def build_contact_topology_repair_report(
    *,
    packets: Sequence[ContactTopologyRepairBenchmarkPacket],
    source_benchmark_file: Path,
    lock_validation: Mapping[str, object],
) -> dict[str, object]:
    baseline_success = sum(
        1 for packet in packets if packet.baseline_visible_partial_success
    )
    repaired_success = sum(
        1 for packet in packets if packet.repaired_visible_partial_success
    )
    baseline_long_mean = _mean(
        packet.baseline_metrics.long_range_contact_recall for packet in packets
    )
    repaired_long_mean = _mean(
        packet.repaired_metrics.long_range_contact_recall for packet in packets
    )
    baseline_beta_mean = _mean(
        packet.baseline_gap.beta_pairing_contact_recall
        for packet in packets
        if packet.baseline_gap.beta_pairing_contact_recall is not None
    )
    repaired_beta_mean = _mean(
        packet.repaired_gap.beta_pairing_contact_recall
        for packet in packets
        if packet.repaired_gap.beta_pairing_contact_recall is not None
    )
    rows = safe_repair_rows(packets)
    return {
        "benchmark_kind": CONTACT_TOPOLOGY_REPAIR_BENCHMARK_KIND,
        "source_visual_mechanism_benchmark_kind": VISUAL_MECHANISM_BENCHMARK_KIND,
        "contact_topology_repair_signature_kind": (
            CONTACT_TOPOLOGY_REPAIR_SIGNATURE_KIND
        ),
        "contact_gap_analysis_signature_kind": CONTACT_GAP_ANALYSIS_SIGNATURE_KIND,
        "repair_signature_kind": LONG_RANGE_CONTACT_REPAIR_SIGNATURE_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "benchmark_size": len(packets),
        "lock_validation": dict(lock_validation),
        "predictor_input_boundary": "sequence_only_no_native_contacts_no_truth_axes",
        "repair_input_boundary": REPAIR_INPUT_BOUNDARY,
        "truth_scoring_boundary": (
            "native_contacts_and_truth_axes_used_only_after_baseline_and_repair"
        ),
        "visual_artifacts_generated_for_rows": len(packets),
        "visual_artifacts_generated_count": (
            len(ROOT_OUTPUT_NAMES) + len(packets) * len(PER_ROW_REPAIR_VISUAL_NAMES)
        ),
        "contact_map_f1_computed_count": len(packets),
        "baseline_visible_partial_success_count": baseline_success,
        "visible_partial_success_count": repaired_success,
        "visible_partial_success_delta": repaired_success - baseline_success,
        "baseline_visible_failure_count": len(packets) - baseline_success,
        "visible_failure_count": len(packets) - repaired_success,
        "baseline_mean_contact_map_f1": _mean(
            packet.baseline_metrics.contact_map_f1 for packet in packets
        ),
        "repaired_mean_contact_map_f1": _mean(
            packet.repaired_metrics.contact_map_f1 for packet in packets
        ),
        "mean_contact_map_f1_delta": round(
            _mean(packet.repaired_metrics.contact_map_f1 for packet in packets)
            - _mean(packet.baseline_metrics.contact_map_f1 for packet in packets),
            6,
        ),
        "baseline_mean_long_range_contact_recall": baseline_long_mean,
        "repaired_mean_long_range_contact_recall": repaired_long_mean,
        "long_range_contact_recall_delta": round(
            repaired_long_mean - baseline_long_mean,
            6,
        ),
        "baseline_mean_beta_pairing_contact_recall": baseline_beta_mean,
        "repaired_mean_beta_pairing_contact_recall": repaired_beta_mean,
        "beta_pairing_contact_recall_delta": round(
            repaired_beta_mean - baseline_beta_mean,
            6,
        ),
        "premature_compaction_count": sum(
            1 for packet in packets if packet.repaired_gap.premature_compaction
        ),
        "false_contact_cluster_count": sum(
            packet.repaired_gap.false_contact_cluster_count for packet in packets
        ),
        "native_cluster_miss_count": sum(
            packet.repaired_gap.native_cluster_miss_count for packet in packets
        ),
        "visual_failure_cohort_count": len(
            {
                packet.repaired_gap.failure_mechanism
                for packet in packets
                if packet.repaired_gap.failure_mechanism != "visible_partial_success"
            }
        ),
        "repair_candidate_count": sum(
            packet.repair.repair_candidate_count for packet in packets
        ),
        "beta_pairing_candidate_count": sum(
            packet.repair.beta_pairing_candidate_count for packet in packets
        ),
        "compact_anchor_candidate_count": sum(
            packet.repair.compact_anchor_candidate_count for packet in packets
        ),
        "local_overclosure_trimmed_count": sum(
            packet.repair.local_overclosure_trimmed_count for packet in packets
        ),
        "native_truth_used_before_prediction": any(
            packet.repair.baseline_prediction.native_truth_used_before_prediction
            for packet in packets
        ),
        "native_truth_used_before_repair": any(
            packet.repair.native_truth_used_before_repair for packet in packets
        ),
        "raw_sequence_exposed": any(packet.repair.raw_sequence_exposed for packet in packets),
        "global_folding_claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "claim_allowed": False,
        "failure_cohorts": {
            row["failure_mechanism"]: row["repaired_row_count"]
            for row in repair_failure_cohort_rows(packets)
        },
        "boundary_statement": (
            "This layer repairs sequence-only contact candidates and analyzes "
            "native gaps after prediction. It improves visual mechanism evidence "
            "without using native truth during prediction or claiming folding is solved."
        ),
        "rows": rows,
    }


def build_contact_topology_repair_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": CONTACT_TOPOLOGY_REPAIR_CERTIFICATE_KIND,
        "benchmark_kind": report["benchmark_kind"],
        "benchmark_size": report["benchmark_size"],
        "contact_map_f1_computed_count": report["contact_map_f1_computed_count"],
        "visual_artifacts_generated_for_rows": report[
            "visual_artifacts_generated_for_rows"
        ],
        "baseline_visible_partial_success_count": report[
            "baseline_visible_partial_success_count"
        ],
        "visible_partial_success_count": report["visible_partial_success_count"],
        "visible_failure_count": report["visible_failure_count"],
        "long_range_contact_recall_delta": report[
            "long_range_contact_recall_delta"
        ],
        "beta_pairing_contact_recall_delta": report[
            "beta_pairing_contact_recall_delta"
        ],
        "native_truth_used_before_prediction": report[
            "native_truth_used_before_prediction"
        ],
        "native_truth_used_before_repair": report[
            "native_truth_used_before_repair"
        ],
        "raw_sequence_exposed": report["raw_sequence_exposed"],
        "global_folding_claim_allowed": report[
            "global_folding_claim_allowed"
        ],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "output_artifacts": tuple(output_names),
    }


def write_contact_topology_repair_outputs(
    *,
    report: Mapping[str, object],
    packets: Sequence[ContactTopologyRepairBenchmarkPacket],
    report_path: Path,
    rows_path: Path,
    gap_analysis_path: Path,
    failure_cohorts_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
    visuals_root: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(safe_repair_rows(packets), rows_path)
    _write_csv_rows(repair_gap_rows(packets), gap_analysis_path)
    _write_csv_rows(repair_failure_cohort_rows(packets), failure_cohorts_path)
    dashboard_path.write_text(render_contact_topology_repair_dashboard(report), encoding="utf-8")
    certificate = build_contact_topology_repair_certificate(
        report,
        output_names=ROOT_OUTPUT_NAMES,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    for packet in packets:
        _write_repair_packet_visuals(packet, visuals_root)
    return (
        report_path,
        rows_path,
        gap_analysis_path,
        failure_cohorts_path,
        dashboard_path,
        certificate_path,
    )


def run_contact_topology_repair_benchmark(
    *,
    benchmark_file: Path,
    report_path: Path,
    rows_path: Path,
    gap_analysis_path: Path,
    failure_cohorts_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
    visuals_root: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    rows = load_visual_mechanism_rows(benchmark_file)
    lock_validation = validate_visual_mechanism_lock(rows)
    packets = contact_topology_repair_packets(rows)
    report = build_contact_topology_repair_report(
        packets=packets,
        source_benchmark_file=benchmark_file,
        lock_validation=lock_validation,
    )
    return write_contact_topology_repair_outputs(
        report=report,
        packets=packets,
        report_path=report_path,
        rows_path=rows_path,
        gap_analysis_path=gap_analysis_path,
        failure_cohorts_path=failure_cohorts_path,
        dashboard_path=dashboard_path,
        certificate_path=certificate_path,
        visuals_root=visuals_root,
    )


def _repair_visual_paths(row_id: str) -> dict[str, str]:
    return {
        name: f"contact_repair_visuals/{row_id}/{name}"
        for name in PER_ROW_REPAIR_VISUAL_NAMES
    }


def _write_repair_packet_visuals(
    packet: ContactTopologyRepairBenchmarkPacket,
    visuals_root: Path,
) -> None:
    row_dir = visuals_root / packet.row.row_id
    repaired_prediction = packet.repair.repaired_prediction
    write_visual_file(
        row_dir / "native_contact_map.svg",
        render_contact_map_svg(
            row_id=packet.row.row_id,
            sequence_length=packet.row.length,
            contact_pairs=packet.row.native_contact_pairs,
            title="Locked coarse native contact target",
            color="#294f9b",
        ),
    )
    write_visual_file(
        row_dir / "repaired_contact_map.svg",
        render_contact_map_svg(
            row_id=packet.row.row_id,
            sequence_length=packet.row.length,
            contact_pairs=repaired_prediction.predicted_contact_pairs,
            title="Repaired sequence-only contact candidates",
            color="#c44b3a",
        ),
    )
    write_visual_file(
        row_dir / "repair_contact_map_overlay.svg",
        render_contact_overlay_svg(
            row_id=packet.row.row_id,
            sequence_length=packet.row.length,
            native_pairs=packet.row.native_contact_pairs,
            predicted_pairs=repaired_prediction.predicted_contact_pairs,
        ),
    )
    write_visual_file(
        row_dir / "repair_folding_trajectory.html",
        render_trajectory_html(
            row_id=packet.row.row_id,
            sequence_length=packet.row.length,
            candidates=repaired_prediction.candidates,
            energy=packet.repaired_energy,
        ),
    )
    write_visual_file(
        row_dir / "repair_energy_curve.svg",
        render_curve_svg(
            packet.repaired_energy.energy_values,
            title="Repaired coarse energy descent curve",
            y_label="relative energy",
        ),
    )
    write_visual_file(
        row_dir / "repair_contact_closure_curve.svg",
        render_curve_svg(
            tuple(float(value) for value in packet.repaired_energy.contact_counts),
            title="Repaired contact closure curve",
            y_label="active contacts",
        ),
    )
    write_visual_file(
        row_dir / "repair_coarse_grain_final.svg",
        render_coarse_grain_svg(
            row_id=packet.row.row_id,
            sequence_length=packet.row.length,
            contacts=repaired_prediction.predicted_contact_pairs,
            title="Repaired coarse-grain final collapse",
        ),
    )


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


def _mean(values: Sequence[float] | object) -> float:
    value_list = [float(value) for value in values]  # type: ignore[arg-type]
    if not value_list:
        return 0.0
    return round(sum(value_list) / len(value_list), 6)


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
        "baseline_visible_partial_success_count",
        "visible_partial_success_count",
        "visible_failure_count",
        "long_range_contact_recall_delta",
        "beta_pairing_contact_recall_delta",
        "premature_compaction_count",
        "native_truth_used_before_repair",
        "raw_sequence_exposed",
        "folding_problem_solved",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _repair_rows_table(report: Mapping[str, object]) -> str:
    rows = report.get("rows", [])
    if not isinstance(rows, Sequence):
        return ""
    body = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        overlay = _escape(row["visual_repair_contact_map_overlay_svg"])
        trajectory = _escape(row["visual_repair_folding_trajectory_html"])
        body.append(
            "<tr>"
            f"<td>{_escape(row['row_id'])}</td>"
            f"<td>{_escape(row['baseline_contact_map_f1'])}</td>"
            f"<td>{_escape(row['repaired_contact_map_f1'])}</td>"
            f"<td>{_escape(row['repaired_long_range_contact_recall'])}</td>"
            f"<td>{_escape(row['repaired_beta_pairing_contact_recall'])}</td>"
            f"<td>{_escape(row['repaired_failure_mechanism'])}</td>"
            f"<td><a href=\"{overlay}\">overlay</a> | "
            f"<a href=\"{trajectory}\">trajectory</a></td>"
            "</tr>"
        )
    return (
        "<section><h2>Native-Gap Repair Rows</h2>"
        "<table><thead><tr>"
        "<th>row</th><th>baseline F1</th><th>repaired F1</th>"
        "<th>long recall</th><th>beta recall</th><th>mechanism</th><th>visuals</th>"
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
            f"<img src=\"{_escape(row['visual_repair_contact_map_overlay_svg'])}\" "
            f"alt=\"repair contact overlay for {_escape(row['row_id'])}\">"
            f"<p>F1: {_escape(row['repaired_contact_map_f1'])} | "
            f"{_escape(row['repaired_failure_mechanism'])}</p>"
            "</article>"
        )
    return (
        "<section><h2>Repaired Contact Map Overlays</h2>"
        "<div class=\"visual-grid\">"
        + "".join(cards)
        + "</div></section>"
    )


def render_contact_topology_repair_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Contact Topology Repair Workbench</title>
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
      color: #59655f;
      font-size: 12px;
      overflow-wrap: anywhere;
    }}
    .metric strong {{
      display: block;
      margin-top: 8px;
      font-size: 24px;
    }}
    .visual-card img {{
      width: 100%;
      aspect-ratio: 1 / 1;
      object-fit: contain;
      background: #f8faf7;
      border: 1px solid #dce4de;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: #ffffff;
      border: 1px solid #d4ddd6;
      border-radius: 6px;
      overflow: hidden;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid #e6ede7;
      text-align: left;
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }}
    th {{
      background: #e8eee8;
      color: #35443e;
    }}
    a {{
      color: #245d62;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Contact Topology Repair Workbench</h1>
    <p>Native-gap analysis explains failures after prediction; repairs remain sequence-only before scoring.</p>
    <div class="metrics">{_metric_cards(report)}</div>
  </header>
  <main>
    <section>
      <h2>Boundary Rules</h2>
      <div class="rules">
        <div class="rule"><strong>Native Gap Analysis Happens After Prediction</strong><br>Missed and false clusters are scored only after sequence-only contact generation.</div>
        <div class="rule"><strong>Long-Range Repair Is Sequence-Only</strong><br>Compact anchors and beta registry candidates use residue composition and length, not native contacts.</div>
        <div class="rule"><strong>Failures Stay Visible</strong><br>Disorder over-collapse, membrane mis-topology, and false-contact clusters remain in the dashboard.</div>
        <div class="rule"><strong>Global Folding Claim Remains Locked</strong><br>This is contact-topology repair, not a solved folding engine.</div>
      </div>
    </section>
    {_visual_grid(report)}
    {_repair_rows_table(report)}
  </main>
</body>
</html>
"""
