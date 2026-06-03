from __future__ import annotations

import csv
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_contact_topology import (
    VISUAL_MECHANISM_BENCHMARK_KIND,
    VISUAL_MECHANISM_SPLIT,
    ContactTopologyPrediction,
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


VISUAL_MECHANISM_SIGNATURE_KIND = "coarse_contact_visual_mechanism_workbench"
VISUAL_MECHANISM_CERTIFICATE_KIND = "visual_folding_mechanism_safety_certificate"

ROOT_OUTPUT_NAMES = (
    "visual_mechanism_12_report.json",
    "visual_mechanism_12_rows.csv",
    "visual_mechanism_12_contact_metrics.csv",
    "visual_mechanism_12_failure_cohorts.csv",
    "visual_mechanism_12_dashboard.html",
    "visual_mechanism_12_certificate.json",
)

PER_ROW_VISUAL_NAMES = (
    "native_contact_map.svg",
    "predicted_contact_map.svg",
    "contact_map_overlay.svg",
    "folding_trajectory.html",
    "energy_curve.svg",
    "contact_closure_curve.svg",
    "coarse_grain_final.svg",
)


@dataclass(frozen=True)
class VisualMechanismPacket:
    row: VisualMechanismRow
    prediction: ContactTopologyPrediction
    metrics: ContactMetricPacket
    energy: EnergyLandscapePacket
    failure_cohort: str
    visible_partial_success: bool
    visual_paths: Mapping[str, str]

    def safe_row(self) -> dict[str, object]:
        axes = self.row.truth_axes
        return {
            "row_id": self.row.row_id,
            "source_id": self.row.source_id,
            "source_kind": self.row.source_kind,
            "sequence_hash": self.row.sequence_sha256,
            "sequence_length": self.row.length,
            "native_scope": self.row.native_scope,
            "mechanism_expected_difficulty": self.row.mechanism_expected_difficulty,
            "truth_secondary_structure_axis": axes.get(
                "secondary_structure_axis",
                "weak_or_unknown",
            ),
            "truth_architecture_axis": axes.get("architecture_axis", "unknown"),
            "truth_order_axis": axes.get("order_axis", "unknown"),
            "truth_environment_axis": axes.get("environment_axis", "unknown"),
            "predicted_contact_count": self.metrics.predicted_contact_count,
            "native_contact_count": self.metrics.native_contact_count,
            "contact_map_f1": self.metrics.contact_map_f1,
            "native_contact_precision": self.metrics.native_contact_precision,
            "native_contact_recall": self.metrics.native_contact_recall,
            "long_range_contact_recall": self.metrics.long_range_contact_recall,
            "short_range_contact_recall": self.metrics.short_range_contact_recall,
            "false_contact_rate": self.metrics.false_contact_rate,
            "trajectory_contact_gain_monotonicity": (
                self.energy.trajectory_contact_gain_monotonicity
            ),
            "energy_descent_consistency": self.energy.energy_descent_consistency,
            "collapse_stability_score": self.energy.collapse_stability_score,
            "visible_partial_success": self.visible_partial_success,
            "failure_cohort": self.failure_cohort,
            "native_contact_map_hash": self.row.native_contact_map_hash,
            "predicted_contact_map_hash": self.prediction.predicted_contact_map_hash,
            "native_truth_used_before_prediction": (
                self.prediction.native_truth_used_before_prediction
            ),
            "raw_sequence_exposed": self.prediction.raw_sequence_exposed,
            "global_folding_claim_allowed": False,
            "folding_problem_solved": False,
            "visual_dir": f"visuals/{self.row.row_id}",
            **{
                f"visual_{name.replace('.', '_')}": path
                for name, path in self.visual_paths.items()
            },
        }

    def contact_metric_row(self) -> dict[str, object]:
        return {
            "row_id": self.row.row_id,
            "source_id": self.row.source_id,
            "sequence_hash": self.row.sequence_sha256,
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


def _visual_paths(row_id: str) -> dict[str, str]:
    return {
        name: f"visuals/{row_id}/{name}"
        for name in PER_ROW_VISUAL_NAMES
    }


def _visible_partial_success(
    metrics: ContactMetricPacket,
    energy: EnergyLandscapePacket,
) -> bool:
    if metrics.contact_map_f1 >= 0.12:
        return True
    return (
        metrics.short_range_contact_recall >= 0.20
        and metrics.native_contact_precision >= 0.08
        and energy.energy_descent_consistency >= 0.80
    )


def _failure_cohort(
    row: VisualMechanismRow,
    metrics: ContactMetricPacket,
    visible_partial_success: bool,
) -> str:
    axes = row.truth_axes
    if visible_partial_success:
        return "visible_partial_success"
    if axes.get("order_axis") == "disordered_flexible":
        return "disorder_control_overclosure_failure"
    if axes.get("environment_axis") == "membrane_like":
        return "membrane_contact_topology_failure"
    if axes.get("architecture_axis") == "multidomain_or_segmented":
        return "architecture_boundary_contact_failure"
    if (
        axes.get("secondary_structure_axis") == "beta_rich"
        and metrics.long_range_contact_recall < 0.25
    ):
        return "beta_long_range_pairing_failure"
    if metrics.contact_map_f1 == 0.0:
        return "no_native_contact_overlap"
    if metrics.false_contact_rate >= 0.85:
        return "false_contact_overprediction"
    return "low_recall_partial_failure"


def visual_mechanism_packets(
    rows: Sequence[VisualMechanismRow],
) -> list[VisualMechanismPacket]:
    packets: list[VisualMechanismPacket] = []
    for row in rows:
        prediction = predict_contact_topology(row.sequence, row_id=row.row_id)
        metrics = evaluate_contact_prediction(
            native_pairs=row.native_contact_pairs,
            predicted_pairs=prediction.predicted_contact_pairs,
        )
        energy = build_energy_landscape(prediction.candidates)
        partial_success = _visible_partial_success(metrics, energy)
        packets.append(
            VisualMechanismPacket(
                row=row,
                prediction=prediction,
                metrics=metrics,
                energy=energy,
                failure_cohort=_failure_cohort(row, metrics, partial_success),
                visible_partial_success=partial_success,
                visual_paths=_visual_paths(row.row_id),
            )
        )
    return packets


def safe_visual_mechanism_rows(
    packets: Sequence[VisualMechanismPacket],
) -> list[dict[str, object]]:
    return [packet.safe_row() for packet in packets]


def contact_metric_rows(
    packets: Sequence[VisualMechanismPacket],
) -> list[dict[str, object]]:
    return [packet.contact_metric_row() for packet in packets]


def failure_cohort_rows(
    packets: Sequence[VisualMechanismPacket],
) -> list[dict[str, object]]:
    counts = Counter(packet.failure_cohort for packet in packets)
    return [
        {
            "failure_cohort": cohort,
            "row_count": counts[cohort],
            "row_ids": ";".join(
                packet.row.row_id
                for packet in packets
                if packet.failure_cohort == cohort
            ),
            "failures_visualized": cohort != "visible_partial_success",
        }
        for cohort in sorted(counts)
    ]


def build_visual_mechanism_report(
    *,
    packets: Sequence[VisualMechanismPacket],
    source_benchmark_file: Path,
    lock_validation: Mapping[str, object],
) -> dict[str, object]:
    rows = safe_visual_mechanism_rows(packets)
    f1_values = [float(row["contact_map_f1"]) for row in rows]
    partial_success_count = sum(
        1 for row in rows if bool(row["visible_partial_success"])
    )
    failure_count = len(rows) - partial_success_count
    output_artifact_count = len(ROOT_OUTPUT_NAMES) + len(rows) * len(PER_ROW_VISUAL_NAMES)
    native_truth_flags = [
        bool(row["native_truth_used_before_prediction"]) for row in rows
    ]
    raw_sequence_flags = [bool(row["raw_sequence_exposed"]) for row in rows]
    return {
        "benchmark_kind": VISUAL_MECHANISM_BENCHMARK_KIND,
        "visual_mechanism_signature_kind": VISUAL_MECHANISM_SIGNATURE_KIND,
        "holdout_split": VISUAL_MECHANISM_SPLIT,
        "source_benchmark_file": str(source_benchmark_file),
        "benchmark_size": len(rows),
        "lock_validation": dict(lock_validation),
        "predictor_input_boundary": "sequence_only_no_native_contacts_no_truth_axes",
        "truth_scoring_boundary": (
            "native_contacts_and_truth_axes_used_only_after_contact_prediction"
        ),
        "visual_artifacts_generated_for_rows": len(rows),
        "visual_artifacts_generated_count": output_artifact_count,
        "visual_files_per_row": len(PER_ROW_VISUAL_NAMES),
        "contact_map_f1_computed_count": len(f1_values),
        "mean_contact_map_f1": (
            _rounded(sum(f1_values) / len(f1_values)) if f1_values else 0.0
        ),
        "max_contact_map_f1": max(f1_values) if f1_values else 0.0,
        "visible_partial_success_count": partial_success_count,
        "visible_failure_count": failure_count,
        "failures_visualized": failure_count > 0,
        "failure_cohort_count": len(set(row["failure_cohort"] for row in rows)),
        "failure_cohorts": {
            row["failure_cohort"]: row["row_count"]
            for row in failure_cohort_rows(packets)
        },
        "native_truth_used_before_prediction": any(native_truth_flags),
        "raw_sequence_exposed": any(raw_sequence_flags),
        "global_folding_claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "claim_allowed": False,
        "boundary_statement": (
            "This workbench renders sequence-only contact candidates against "
            "locked coarse native-contact targets after prediction. It creates "
            "visual mechanism evidence, not a solved folding engine."
        ),
        "rows": rows,
    }


def build_visual_mechanism_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": VISUAL_MECHANISM_CERTIFICATE_KIND,
        "benchmark_kind": report["benchmark_kind"],
        "visual_mechanism_signature_kind": report[
            "visual_mechanism_signature_kind"
        ],
        "holdout_split": report["holdout_split"],
        "benchmark_size": report["benchmark_size"],
        "visual_artifacts_generated_for_rows": report[
            "visual_artifacts_generated_for_rows"
        ],
        "contact_map_f1_computed_count": report["contact_map_f1_computed_count"],
        "visible_partial_success_count": report["visible_partial_success_count"],
        "visible_failure_count": report["visible_failure_count"],
        "failures_visualized": report["failures_visualized"],
        "native_truth_used_before_prediction": report[
            "native_truth_used_before_prediction"
        ],
        "raw_sequence_exposed": report["raw_sequence_exposed"],
        "global_folding_claim_allowed": report[
            "global_folding_claim_allowed"
        ],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "output_artifacts": tuple(output_names),
    }


def write_visual_mechanism_outputs(
    *,
    report: Mapping[str, object],
    packets: Sequence[VisualMechanismPacket],
    report_path: Path,
    rows_path: Path,
    contact_metrics_path: Path,
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
    _write_csv_rows(safe_visual_mechanism_rows(packets), rows_path)
    _write_csv_rows(contact_metric_rows(packets), contact_metrics_path)
    _write_csv_rows(failure_cohort_rows(packets), failure_cohorts_path)
    dashboard_path.write_text(render_visual_mechanism_dashboard(report), encoding="utf-8")
    certificate = build_visual_mechanism_certificate(
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
        failure_cohorts_path,
        dashboard_path,
        certificate_path,
    )


def run_visual_mechanism_benchmark(
    *,
    benchmark_file: Path,
    report_path: Path,
    rows_path: Path,
    contact_metrics_path: Path,
    failure_cohorts_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
    visuals_root: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    rows = load_visual_mechanism_rows(benchmark_file)
    lock_validation = validate_visual_mechanism_lock(rows)
    packets = visual_mechanism_packets(rows)
    report = build_visual_mechanism_report(
        packets=packets,
        source_benchmark_file=benchmark_file,
        lock_validation=lock_validation,
    )
    return write_visual_mechanism_outputs(
        report=report,
        packets=packets,
        report_path=report_path,
        rows_path=rows_path,
        contact_metrics_path=contact_metrics_path,
        failure_cohorts_path=failure_cohorts_path,
        dashboard_path=dashboard_path,
        certificate_path=certificate_path,
        visuals_root=visuals_root,
    )


def _write_packet_visuals(packet: VisualMechanismPacket, visuals_root: Path) -> None:
    row_dir = visuals_root / packet.row.row_id
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
        row_dir / "predicted_contact_map.svg",
        render_contact_map_svg(
            row_id=packet.row.row_id,
            sequence_length=packet.row.length,
            contact_pairs=packet.prediction.predicted_contact_pairs,
            title="Sequence-only predicted contact candidates",
            color="#c44b3a",
        ),
    )
    write_visual_file(
        row_dir / "contact_map_overlay.svg",
        render_contact_overlay_svg(
            row_id=packet.row.row_id,
            sequence_length=packet.row.length,
            native_pairs=packet.row.native_contact_pairs,
            predicted_pairs=packet.prediction.predicted_contact_pairs,
        ),
    )
    write_visual_file(
        row_dir / "folding_trajectory.html",
        render_trajectory_html(
            row_id=packet.row.row_id,
            sequence_length=packet.row.length,
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
            sequence_length=packet.row.length,
            contacts=packet.prediction.predicted_contact_pairs,
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
        "visual_artifacts_generated_for_rows",
        "contact_map_f1_computed_count",
        "visible_partial_success_count",
        "visible_failure_count",
        "mean_contact_map_f1",
        "native_truth_used_before_prediction",
        "raw_sequence_exposed",
        "global_folding_claim_allowed",
        "folding_problem_solved",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
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
        trajectory = _escape(row["visual_folding_trajectory_html"])
        body.append(
            "<tr>"
            f"<td>{_escape(row['row_id'])}</td>"
            f"<td>{_escape(row['truth_secondary_structure_axis'])}</td>"
            f"<td>{_escape(row['truth_architecture_axis'])}</td>"
            f"<td>{_escape(row['contact_map_f1'])}</td>"
            f"<td>{_escape(row['failure_cohort'])}</td>"
            f"<td><a href=\"{overlay}\">overlay</a> | "
            f"<a href=\"{trajectory}\">trajectory</a></td>"
            "</tr>"
        )
    return (
        "<section><h2>Visual Rows</h2>"
        "<table><thead><tr>"
        "<th>row</th><th>secondary</th><th>architecture</th>"
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
            f"alt=\"contact overlay for {_escape(row['row_id'])}\">"
            f"<p>F1: {_escape(row['contact_map_f1'])} | "
            f"{_escape(row['failure_cohort'])}</p>"
            "</article>"
        )
    return (
        "<section><h2>Contact Map Overlays</h2>"
        "<div class=\"visual-grid\">"
        + "".join(cards)
        + "</div></section>"
    )


def render_visual_mechanism_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Visual Folding Mechanism Workbench</title>
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
    <h1>Visual Folding Mechanism Workbench</h1>
    <p>This Is A Mechanism Visualization, Not A Solved Folding Engine</p>
    <div class="metrics">{_metric_cards(report)}</div>
  </header>
  <main>
    <section>
      <h2>Boundary Rules</h2>
      <div class="rules">
        <div class="rule"><strong>Contact Maps Are The First Proof Target</strong><br>Scoring starts with coarse native-contact overlap, not a global fold label.</div>
        <div class="rule"><strong>Trajectory Is Coarse-Grained</strong><br>The path is a visual reconstruction of candidate closure, not atomistic dynamics.</div>
        <div class="rule"><strong>Native Contacts Are Used Only After Prediction</strong><br>Prediction is sequence-only; locked contact targets enter only during scoring and rendering.</div>
        <div class="rule"><strong>Failures Are Visual Evidence</strong><br>Low overlap and overclosure are kept visible instead of hidden behind a single label.</div>
        <div class="rule"><strong>Global Folding Claim Remains Locked</strong><br>The workbench refuses solved-folding, clinical, molecule, and protein-design claims.</div>
      </div>
    </section>
    {_visual_grid(report)}
    {_dashboard_preview_rows(report)}
  </main>
</body>
</html>
"""
