from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping, Sequence


VISUAL_MECHANISM_AUDIT_KIND = "visual_mechanism_claim_audit_v1"
VISUAL_MECHANISM_AUDIT_CERTIFICATE_KIND = "visual_mechanism_audit_certificate"

ROOT_OUTPUT_NAMES = (
    "visual_mechanism_audit_report.json",
    "visual_mechanism_audit_rows.csv",
    "visual_mechanism_audit_overfit_risks.csv",
    "visual_mechanism_audit_dashboard.html",
    "visual_mechanism_audit_certificate.json",
)


def load_json_object(path: Path) -> dict[str, object]:
    parsed = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return parsed


def build_visual_mechanism_audit_report(
    *,
    baseline_report_path: Path,
    repair_report_path: Path,
    source_benchmark_file: Path,
) -> dict[str, object]:
    baseline_report = load_json_object(baseline_report_path)
    repair_report = load_json_object(repair_report_path)
    rows = audit_rows(
        baseline_rows=_safe_rows(baseline_report),
        repair_rows=_safe_rows(repair_report),
    )
    overfit_rows = overfit_risk_rows(rows)
    return {
        "audit_kind": VISUAL_MECHANISM_AUDIT_KIND,
        "baseline_report_file": str(baseline_report_path),
        "repair_report_file": str(repair_report_path),
        "source_benchmark_file": str(source_benchmark_file),
        "visual_12_is_toy_benchmark": True,
        "visual_12_scope": "toy_coarse_internal_contact_map_benchmark",
        "coarse_native_contacts_only": True,
        "full_physical_native_contacts_available": False,
        "baseline_visual_rows_rendered": baseline_report.get(
            "visual_artifacts_generated_for_rows",
            0,
        ),
        "repaired_visual_rows_rendered": repair_report.get(
            "visual_artifacts_generated_for_rows",
            0,
        ),
        "baseline_visible_partial_success_count": baseline_report.get(
            "visible_partial_success_count",
            0,
        ),
        "repaired_visible_partial_success_count": repair_report.get(
            "visible_partial_success_count",
            0,
        ),
        "visible_partial_success_delta": repair_report.get(
            "visible_partial_success_delta",
            0,
        ),
        "baseline_visible_failure_count": baseline_report.get(
            "visible_failure_count",
            0,
        ),
        "repaired_visible_failure_count": repair_report.get(
            "visible_failure_count",
            0,
        ),
        "baseline_mean_contact_map_f1": baseline_report.get(
            "mean_contact_map_f1",
            0,
        ),
        "repaired_mean_contact_map_f1": repair_report.get(
            "repaired_mean_contact_map_f1",
            0,
        ),
        "long_range_contact_recall_delta": repair_report.get(
            "long_range_contact_recall_delta",
            0,
        ),
        "beta_pairing_contact_recall_delta": repair_report.get(
            "beta_pairing_contact_recall_delta",
            0,
        ),
        "hardcoded_beta_registry_pair_templates_detected": True,
        "hardcoded_beta_registry_pattern_families": (
            "compact_beta_registry_centers",
            "fixed_cross_sheet_anchors",
        ),
        "contact_repair_overfit_risk_reported": True,
        "contact_repair_overfit_risk_level": "high_for_visual_12_benchmark",
        "overfit_risk_row_count": len(overfit_rows),
        "beta_template_success_gain_row_count": sum(
            1
            for row in rows
            if row["baseline_visible_partial_success"] is False
            and row["repaired_visible_partial_success"] is True
            and int(row["beta_pairing_candidate_count"]) > 0
        ),
        "mechanism_discovery_claim_allowed": False,
        "mechanism_discovery_claim_created": False,
        "folding_problem_solved": False,
        "global_folding_claim_allowed": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "raw_sequence_exposed": False,
        "native_truth_used_before_prediction": False,
        "native_truth_used_before_repair": False,
        "artifact_reproducible": True,
        "clean_archive_required": True,
        "clean_archive_command": (
            "git archive --format=zip --output /private/tmp/pharmacotopology-clean.zip HEAD"
        ),
        "finder_zip_allowed": False,
        "claim_boundary": (
            "The visual and repaired 12-row artifacts are toy, coarse, internal "
            "contact-map benchmarks. They support inspection of contact "
            "hypotheses and failure modes, not mechanism discovery or solved "
            "protein folding."
        ),
        "rows": rows,
        "overfit_risk_rows": overfit_rows,
    }


def audit_rows(
    *,
    baseline_rows: Sequence[Mapping[str, object]],
    repair_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    baseline_by_id = {str(row["row_id"]): row for row in baseline_rows}
    output = []
    for repair_row in repair_rows:
        row_id = str(repair_row["row_id"])
        baseline_row = baseline_by_id[row_id]
        beta_candidate_count = int(repair_row.get("beta_pairing_candidate_count", 0))
        compact_candidate_count = int(
            repair_row.get("compact_anchor_candidate_count", 0)
        )
        baseline_success = _as_bool(
            baseline_row.get("visible_partial_success", False)
        )
        repaired_success = _as_bool(
            repair_row.get("repaired_visible_partial_success", False)
        )
        risk_reason = _overfit_risk_reason(
            beta_candidate_count=beta_candidate_count,
            compact_candidate_count=compact_candidate_count,
            baseline_success=baseline_success,
            repaired_success=repaired_success,
            repaired_failure_mechanism=str(
                repair_row.get("repaired_failure_mechanism", "")
            ),
        )
        output.append(
            {
                "row_id": row_id,
                "truth_secondary_structure_axis": repair_row[
                    "truth_secondary_structure_axis"
                ],
                "truth_architecture_axis": repair_row["truth_architecture_axis"],
                "truth_order_axis": repair_row["truth_order_axis"],
                "truth_environment_axis": repair_row["truth_environment_axis"],
                "baseline_contact_map_f1": baseline_row["contact_map_f1"],
                "repaired_contact_map_f1": repair_row["repaired_contact_map_f1"],
                "contact_map_f1_delta": repair_row["contact_map_f1_delta"],
                "baseline_visible_partial_success": baseline_success,
                "repaired_visible_partial_success": repaired_success,
                "baseline_failure_cohort": baseline_row["failure_cohort"],
                "repaired_failure_mechanism": repair_row[
                    "repaired_failure_mechanism"
                ],
                "beta_pairing_candidate_count": beta_candidate_count,
                "compact_anchor_candidate_count": compact_candidate_count,
                "local_overclosure_trimmed_count": int(
                    repair_row.get("local_overclosure_trimmed_count", 0)
                ),
                "toy_benchmark_claim_only": True,
                "mechanism_discovery_claim_allowed": False,
                "overfit_risk_flag": bool(risk_reason),
                "overfit_risk_reason": risk_reason,
            }
        )
    return output


def overfit_risk_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    keys = (
        "row_id",
        "truth_secondary_structure_axis",
        "baseline_contact_map_f1",
        "repaired_contact_map_f1",
        "contact_map_f1_delta",
        "baseline_visible_partial_success",
        "repaired_visible_partial_success",
        "beta_pairing_candidate_count",
        "compact_anchor_candidate_count",
        "repaired_failure_mechanism",
        "overfit_risk_reason",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if bool(row["overfit_risk_flag"])
    ]


def build_visual_mechanism_audit_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": VISUAL_MECHANISM_AUDIT_CERTIFICATE_KIND,
        "audit_kind": report["audit_kind"],
        "visual_12_is_toy_benchmark": report["visual_12_is_toy_benchmark"],
        "coarse_native_contacts_only": report["coarse_native_contacts_only"],
        "contact_repair_overfit_risk_reported": report[
            "contact_repair_overfit_risk_reported"
        ],
        "mechanism_discovery_claim_allowed": report[
            "mechanism_discovery_claim_allowed"
        ],
        "folding_problem_solved": report["folding_problem_solved"],
        "global_folding_claim_allowed": report["global_folding_claim_allowed"],
        "artifact_reproducible": report["artifact_reproducible"],
        "clean_archive_required": report["clean_archive_required"],
        "finder_zip_allowed": report["finder_zip_allowed"],
        "raw_sequence_exposed": report["raw_sequence_exposed"],
        "native_truth_used_before_prediction": report[
            "native_truth_used_before_prediction"
        ],
        "native_truth_used_before_repair": report[
            "native_truth_used_before_repair"
        ],
        "output_artifacts": tuple(output_names),
    }


def write_visual_mechanism_audit_outputs(
    *,
    report: Mapping[str, object],
    report_path: Path,
    rows_path: Path,
    overfit_risks_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(report["rows"], rows_path)  # type: ignore[arg-type]
    _write_csv_rows(
        report["overfit_risk_rows"],  # type: ignore[arg-type]
        overfit_risks_path,
    )
    dashboard_path.write_text(render_visual_mechanism_audit_dashboard(report), encoding="utf-8")
    certificate = build_visual_mechanism_audit_certificate(
        report,
        output_names=ROOT_OUTPUT_NAMES,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        rows_path,
        overfit_risks_path,
        dashboard_path,
        certificate_path,
    )


def _safe_rows(report: Mapping[str, object]) -> Sequence[Mapping[str, object]]:
    rows = report.get("rows", ())
    if not isinstance(rows, Sequence):
        raise ValueError("report rows must be a sequence")
    return tuple(row for row in rows if isinstance(row, Mapping))


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).lower() == "true"


def _overfit_risk_reason(
    *,
    beta_candidate_count: int,
    compact_candidate_count: int,
    baseline_success: bool,
    repaired_success: bool,
    repaired_failure_mechanism: str,
) -> str:
    if beta_candidate_count and not baseline_success and repaired_success:
        return "beta_registry_template_created_success_gain_on_toy_benchmark"
    if beta_candidate_count and repaired_failure_mechanism in {
        "bad_beta_pairing",
        "premature_compaction",
    }:
        return "beta_registry_template_still_drives_failure_or_compaction"
    if compact_candidate_count and not baseline_success and repaired_success:
        return "compact_anchor_template_created_success_gain_on_toy_benchmark"
    return ""


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
        "visual_12_is_toy_benchmark",
        "baseline_visible_partial_success_count",
        "repaired_visible_partial_success_count",
        "overfit_risk_row_count",
        "contact_repair_overfit_risk_reported",
        "mechanism_discovery_claim_allowed",
        "folding_problem_solved",
        "artifact_reproducible",
        "clean_archive_required",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _audit_rows_table(report: Mapping[str, object]) -> str:
    rows = report.get("rows", [])
    if not isinstance(rows, Sequence):
        return ""
    body = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        body.append(
            "<tr>"
            f"<td>{_escape(row['row_id'])}</td>"
            f"<td>{_escape(row['baseline_contact_map_f1'])}</td>"
            f"<td>{_escape(row['repaired_contact_map_f1'])}</td>"
            f"<td>{_escape(row['repaired_failure_mechanism'])}</td>"
            f"<td>{_escape(row['overfit_risk_reason'])}</td>"
            "</tr>"
        )
    return (
        "<section><h2>Audit Rows</h2>"
        "<table><thead><tr>"
        "<th>row</th><th>baseline F1</th><th>repaired F1</th>"
        "<th>repair outcome</th><th>overfit risk</th>"
        "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


def render_visual_mechanism_audit_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Visual Mechanism Audit</title>
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
      max-width: 1180px;
      margin: 0 auto;
      padding: 24px;
    }}
    h1, h2 {{
      margin: 0 0 12px;
      letter-spacing: 0;
    }}
    section {{
      margin: 24px 0;
    }}
    .metrics, .rules {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
    }}
    .metric, .rule {{
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
  </style>
</head>
<body>
  <header>
    <h1>Visual Mechanism Audit</h1>
    <p>The 12-row visual benchmark is toy, coarse, internal, and not evidence that folding is solved.</p>
    <div class="metrics">{_metric_cards(report)}</div>
  </header>
  <main>
    <section>
      <h2>Freeze Rules</h2>
      <div class="rules">
        <div class="rule"><strong>Toy Benchmark</strong><br>The 12-row visual set is internal and coarse, with locked contact targets rather than full physical truth.</div>
        <div class="rule"><strong>Overfit Risk Reported</strong><br>Hardcoded beta registry patterns improved the tiny benchmark and must not be treated as mechanism discovery.</div>
        <div class="rule"><strong>No Discovery Claim</strong><br>Mechanism discovery, global folding, clinical, molecule, and protein-design claims remain refused.</div>
        <div class="rule"><strong>Archive Rule</strong><br>Use git archive from HEAD for sharing; Finder zips are not accepted as reproducibility evidence.</div>
      </div>
    </section>
    {_audit_rows_table(report)}
  </main>
</body>
</html>
"""
