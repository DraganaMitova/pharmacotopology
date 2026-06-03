from __future__ import annotations

import argparse
import html
import json
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_structure_features import (  # noqa: E402
    FOLD_CLASS_TO_STRATUM,
    similarity_distribution,
    topology_signature_delta,
)
from pharmacotopology.folding_topology import FOLDING_TOPOLOGY_DIMENSIONS  # noqa: E402


DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_500_report.json"
)
DEFAULT_CSV_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_500_rows.csv"
)
DEFAULT_OUTPUT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/real_folding_500_dashboard.html"
)


def _read_json(path: Path) -> dict[str, Any]:
    parsed = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    return parsed


def _escape(value: object) -> str:
    return html.escape(str(value), quote=True)


def _comparisons(report: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        row
        for row in report.get("comparisons", [])
        if isinstance(row, dict)
    ]


def _confusion_matrix(rows: Sequence[Mapping[str, Any]]) -> dict[str, dict[str, int]]:
    matrix: dict[str, dict[str, int]] = {}
    for row in rows:
        actual = str(row.get("reference_fold_class", "unknown"))
        predicted = str(row.get("predicted_fold_class", "unknown"))
        matrix.setdefault(actual, {})
        matrix[actual][predicted] = matrix[actual].get(predicted, 0) + 1
    return {key: dict(sorted(value.items())) for key, value in sorted(matrix.items())}


def _per_class_accuracy(rows: Sequence[Mapping[str, Any]]) -> dict[str, float]:
    totals: dict[str, int] = {}
    matches: dict[str, int] = {}
    for row in rows:
        actual = str(row.get("reference_fold_class", "unknown"))
        totals[actual] = totals.get(actual, 0) + 1
        if bool(row.get("fold_class_match")):
            matches[actual] = matches.get(actual, 0) + 1
    return {
        fold_class: round(matches.get(fold_class, 0) / total, 6)
        for fold_class, total in sorted(totals.items())
    }


def _similarity_distribution_from_rows(
    rows: Sequence[Mapping[str, Any]],
) -> dict[str, int]:
    class _ComparisonProxy:
        def __init__(self, similarity: float) -> None:
            self.contact_map_similarity = similarity

    return similarity_distribution(
        tuple(
            _ComparisonProxy(float(row.get("contact_map_similarity", 0.0)))
            for row in rows
        )
    )


def _worst_rows(
    rows: Sequence[Mapping[str, Any]],
    *,
    limit: int = 20,
) -> list[Mapping[str, Any]]:
    return sorted(
        rows,
        key=lambda row: (
            bool(row.get("fold_class_match")),
            float(row.get("contact_map_similarity", 0.0)),
        ),
    )[:limit]


def _metric_cards(report: Mapping[str, Any]) -> str:
    metrics = [
        (
            "Benchmark size",
            report.get("benchmark_size", report.get("comparisons_reviewed", 0)),
        ),
        ("External rows", report.get("external_rows", 0)),
        (
            "Accuracy",
            report.get("accuracy", report.get("fold_class_match_rate", 0.0)),
        ),
        ("Fold match rate", report.get("fold_class_match_rate", 0.0)),
        ("Mean similarity", report.get("mean_contact_map_similarity", 0.0)),
        ("Perfect matches", report.get("perfect_matches", 0)),
        ("Mismatches", report.get("mismatches", 0)),
        ("Locked", report.get("locked_after_generation", False)),
    ]
    return "\n".join(
        (
            "<div class=\"metric\">"
            f"<span>{_escape(label)}</span>"
            f"<strong>{_escape(value)}</strong>"
            "</div>"
        )
        for label, value in metrics
    )


def _render_confusion_matrix(rows: Sequence[Mapping[str, Any]]) -> str:
    matrix = _confusion_matrix(rows)
    predicted_classes = sorted(
        {
            predicted
            for predictions in matrix.values()
            for predicted in predictions
        }
    )
    if not matrix:
        return "<p class=\"empty\">No benchmark rows loaded.</p>"
    header = "".join(f"<th>{_escape(item)}</th>" for item in predicted_classes)
    body_rows = []
    for actual, predictions in matrix.items():
        cells = "".join(
            f"<td>{predictions.get(predicted, 0)}</td>"
            for predicted in predicted_classes
        )
        body_rows.append(f"<tr><th>{_escape(actual)}</th>{cells}</tr>")
    return (
        "<table><thead><tr><th>Actual \\ Predicted</th>"
        + header
        + "</tr></thead><tbody>"
        + "".join(body_rows)
        + "</tbody></table>"
    )


def _render_distribution(rows: Sequence[Mapping[str, Any]]) -> str:
    distribution = _similarity_distribution_from_rows(rows)
    total = max(sum(distribution.values()), 1)
    bars = []
    for bucket, count in distribution.items():
        width = round((count / total) * 100, 3)
        bars.append(
            "<div class=\"bar-row\">"
            f"<span>{_escape(bucket)}</span>"
            "<div class=\"bar-track\">"
            f"<div class=\"bar\" style=\"width:{width}%\"></div>"
            "</div>"
            f"<strong>{count}</strong>"
            "</div>"
        )
    return "".join(bars)


def _render_accuracy(rows: Sequence[Mapping[str, Any]]) -> str:
    accuracy = _per_class_accuracy(rows)
    if not accuracy:
        return "<p class=\"empty\">No per-class accuracy yet.</p>"
    return "".join(
        (
            "<div class=\"bar-row\">"
            f"<span>{_escape(fold_class)}</span>"
            "<div class=\"bar-track\">"
            f"<div class=\"bar accent\" style=\"width:{round(value * 100, 3)}%\"></div>"
            "</div>"
            f"<strong>{value:.3f}</strong>"
            "</div>"
        )
        for fold_class, value in accuracy.items()
    )


def _render_failure_table(rows: Sequence[Mapping[str, Any]]) -> str:
    worst = _worst_rows(rows)
    if not worst:
        return "<p class=\"empty\">No benchmark rows loaded.</p>"
    body = []
    for row in worst:
        body.append(
            "<tr>"
            f"<td>{_escape(row.get('protein_id', ''))}</td>"
            f"<td>{_escape(row.get('reference_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('predicted_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('contact_map_similarity', 0.0))}</td>"
            f"<td>{_escape(row.get('failure_reason', ''))}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr>"
        "<th>protein_id</th><th>actual</th><th>predicted</th>"
        "<th>similarity</th><th>failure_reason</th>"
        "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table>"
    )


def _point(cx: float, cy: float, radius: float, index: int, count: int, value: float) -> str:
    from math import cos, pi, sin

    angle = (-pi / 2) + (2 * pi * index / count)
    return f"{cx + cos(angle) * radius * value:.3f},{cy + sin(angle) * radius * value:.3f}"


def _render_radar(rows: Sequence[Mapping[str, Any]]) -> str:
    if not rows:
        return "<p class=\"empty\">No radar overlay yet.</p>"
    row = _worst_rows(rows, limit=1)[0]
    predicted = row.get("predicted_topology_signature", {})
    reference = row.get("reference_topology_signature", {})
    if not isinstance(predicted, Mapping) or not isinstance(reference, Mapping):
        return "<p class=\"empty\">Missing signature data.</p>"
    dimensions = FOLDING_TOPOLOGY_DIMENSIONS
    cx = cy = 170.0
    radius = 120.0
    predicted_points = " ".join(
        _point(cx, cy, radius, index, len(dimensions), float(predicted.get(dim, 0.0)))
        for index, dim in enumerate(dimensions)
    )
    reference_points = " ".join(
        _point(cx, cy, radius, index, len(dimensions), float(reference.get(dim, 0.0)))
        for index, dim in enumerate(dimensions)
    )
    axes = []
    for index, dimension in enumerate(dimensions):
        end = _point(cx, cy, radius, index, len(dimensions), 1.0)
        label = _point(cx, cy, radius + 26, index, len(dimensions), 1.0)
        x, y = label.split(",")
        axes.append(
            f"<line x1=\"{cx}\" y1=\"{cy}\" x2=\"{end.split(',')[0]}\" y2=\"{end.split(',')[1]}\" />"
            f"<text x=\"{x}\" y=\"{y}\">{_escape(dimension.replace('_', ' '))}</text>"
        )
    delta = topology_signature_delta(row)
    largest_delta = max(delta.items(), key=lambda item: abs(item[1]))
    return (
        "<div class=\"radar-wrap\">"
        "<svg viewBox=\"0 0 340 340\" role=\"img\" aria-label=\"Radar overlay\">"
        "<circle cx=\"170\" cy=\"170\" r=\"120\"></circle>"
        "<circle cx=\"170\" cy=\"170\" r=\"80\"></circle>"
        "<circle cx=\"170\" cy=\"170\" r=\"40\"></circle>"
        + "".join(axes)
        + f"<polygon class=\"reference\" points=\"{reference_points}\"></polygon>"
        + f"<polygon class=\"predicted\" points=\"{predicted_points}\"></polygon>"
        + "</svg>"
        "<div>"
        f"<p><strong>{_escape(row.get('protein_id', ''))}</strong></p>"
        "<p class=\"legend\"><span class=\"swatch ref\"></span>reference "
        "<span class=\"swatch pred\"></span>predicted</p>"
        f"<p>Largest signature delta: {_escape(largest_delta[0])} = {_escape(largest_delta[1])}</p>"
        "</div></div>"
    )


def _render_certificate(report: Mapping[str, Any], report_path: Path, csv_path: Path) -> str:
    certificate = report.get("lock_certificate", {})
    validation = report.get("reference_dataset_validation", {})
    if not isinstance(certificate, Mapping):
        certificate = {}
    if not isinstance(validation, Mapping):
        validation = {}
    blockers = certificate.get("lock_blockers", [])
    if not isinstance(blockers, list):
        blockers = []
    blocker_items = "".join(f"<li>{_escape(item)}</li>" for item in blockers) or "<li>none recorded</li>"
    return (
        "<dl class=\"certificate\">"
        f"<dt>report</dt><dd>{_escape(report_path)}</dd>"
        f"<dt>csv</dt><dd>{_escape(csv_path)}</dd>"
        f"<dt>dataset hash</dt><dd>{_escape(certificate.get('dataset_hash', 'not attached'))}</dd>"
        f"<dt>commit hash</dt><dd>{_escape(certificate.get('recipe_commit_hash', 'not attached'))}</dd>"
        f"<dt>locked_after_generation</dt><dd>{_escape(report.get('locked_after_generation', False))}</dd>"
        f"<dt>no_retuning_flag</dt><dd>{_escape(report.get('no_retuning_flag', False))}</dd>"
        f"<dt>external rows</dt><dd>{_escape(validation.get('external_reference_count', 0))}</dd>"
        "</dl>"
        "<h3>Lock Blockers</h3><ul class=\"blockers\">"
        + blocker_items
        + "</ul>"
    )


def render_folding_benchmark_dashboard(
    report_path: Path = DEFAULT_REPORT_PATH,
    csv_path: Path = DEFAULT_CSV_PATH,
    output_path: Path = DEFAULT_OUTPUT_PATH,
) -> Path:
    report = _read_json(report_path)
    rows = _comparisons(report)
    source_labels = ", ".join(str(item) for item in report.get("benchmark_sources", []))
    strata = ", ".join(
        f"{key} -> {value}"
        for key, value in sorted(FOLD_CLASS_TO_STRATUM.items())
    )
    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Folding Benchmark Dashboard</title>
<style>
:root {{
  color-scheme: light;
  --ink: #18201d;
  --muted: #5d6762;
  --line: #cad6d0;
  --panel: #f7faf8;
  --accent: #0d766e;
  --warn: #9b3d2e;
  --blue: #285f9f;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  color: var(--ink);
  background: #ffffff;
}}
header {{
  padding: 28px clamp(18px, 4vw, 54px);
  border-bottom: 1px solid var(--line);
  background: #eef5f2;
}}
h1 {{ margin: 0 0 8px; font-size: clamp(28px, 4vw, 48px); letter-spacing: 0; }}
h2 {{ margin: 0 0 16px; font-size: 22px; letter-spacing: 0; }}
h3 {{ margin: 18px 0 8px; font-size: 16px; letter-spacing: 0; }}
p {{ max-width: 980px; line-height: 1.55; }}
main {{ padding: 24px clamp(18px, 4vw, 54px) 54px; }}
section {{ margin: 0 0 34px; }}
.warning {{
  display: inline-block;
  padding: 8px 10px;
  border: 1px solid var(--warn);
  color: var(--warn);
  font-weight: 700;
  background: #fff8f6;
}}
.metrics {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
  gap: 10px;
  margin-top: 18px;
}}
.metric {{
  border: 1px solid var(--line);
  background: var(--panel);
  padding: 12px;
  min-height: 82px;
}}
.metric span {{ display: block; color: var(--muted); font-size: 13px; }}
.metric strong {{ display: block; margin-top: 8px; font-size: 22px; }}
table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 14px;
}}
th, td {{
  border: 1px solid var(--line);
  padding: 8px;
  text-align: left;
  vertical-align: top;
}}
th {{ background: #eef5f2; }}
.bar-row {{
  display: grid;
  grid-template-columns: minmax(150px, 260px) 1fr 70px;
  gap: 10px;
  align-items: center;
  margin: 10px 0;
}}
.bar-track {{ height: 16px; border: 1px solid var(--line); background: #fff; }}
.bar {{ height: 100%; background: var(--blue); }}
.bar.accent {{ background: var(--accent); }}
.grid {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
  gap: 26px;
}}
.radar-wrap {{ display: grid; grid-template-columns: 360px 1fr; gap: 20px; align-items: center; }}
svg {{ max-width: 100%; }}
svg circle, svg line {{ stroke: var(--line); fill: none; }}
svg text {{ font-size: 8px; fill: var(--muted); text-anchor: middle; }}
polygon.reference {{ fill: rgba(13, 118, 110, 0.20); stroke: var(--accent); stroke-width: 2; }}
polygon.predicted {{ fill: rgba(40, 95, 159, 0.20); stroke: var(--blue); stroke-width: 2; }}
.swatch {{ display: inline-block; width: 12px; height: 12px; margin: 0 4px 0 10px; vertical-align: -1px; }}
.swatch.ref {{ background: var(--accent); }}
.swatch.pred {{ background: var(--blue); }}
.certificate {{
  display: grid;
  grid-template-columns: minmax(160px, 260px) 1fr;
  border: 1px solid var(--line);
}}
.certificate dt, .certificate dd {{ margin: 0; padding: 8px; border-bottom: 1px solid var(--line); }}
.certificate dt {{ background: #eef5f2; font-weight: 700; }}
.empty {{ color: var(--muted); font-style: italic; }}
.blockers {{ margin-top: 0; }}
@media (max-width: 760px) {{
  .radar-wrap {{ grid-template-columns: 1fr; }}
  .bar-row {{ grid-template-columns: 1fr; }}
  .certificate {{ grid-template-columns: 1fr; }}
}}
</style>
</head>
<body>
<header>
  <p class="warning">NOT A FOLDING SOLUTION / NOT A DRUG-DESIGN TOOL</p>
  <h1>Folding Benchmark Dashboard</h1>
  <p>This page visualizes external-alignment evidence and failures. It does not
  claim that protein folding has been solved, and it should become more useful
  as real external rows replace placeholder shells.</p>
  <div class="metrics">{_metric_cards(report)}</div>
</header>
<main>
  <section>
    <h2>Locked Benchmark Certificate</h2>
    {_render_certificate(report, report_path, csv_path)}
  </section>
  <section>
    <h2>Confusion Matrix</h2>
    {_render_confusion_matrix(rows)}
  </section>
  <section class="grid">
    <div>
      <h2>Similarity Distribution</h2>
      {_render_distribution(rows)}
    </div>
    <div>
      <h2>Per-Class Accuracy</h2>
      {_render_accuracy(rows)}
    </div>
  </section>
  <section>
    <h2>Failure Table</h2>
    {_render_failure_table(rows)}
  </section>
  <section>
    <h2>Radar Overlay</h2>
    {_render_radar(rows)}
  </section>
  <section>
    <h2>Benchmark Source Notes</h2>
    <p>Sources: {_escape(source_labels or "not attached")}</p>
    <p>Fold-class strata: {_escape(strata)}</p>
  </section>
</main>
</body>
</html>
"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html_text, encoding="utf-8")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render a static folding benchmark proof dashboard."
    )
    parser.add_argument(
        "--report",
        default=str(DEFAULT_REPORT_PATH),
        help="Benchmark JSON report path.",
    )
    parser.add_argument(
        "--csv",
        default=str(DEFAULT_CSV_PATH),
        help="Benchmark CSV path for certificate display.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_PATH),
        help="Output HTML dashboard path.",
    )
    args = parser.parse_args()
    print(
        render_folding_benchmark_dashboard(
            Path(args.report),
            Path(args.csv),
            Path(args.output),
        )
    )


if __name__ == "__main__":
    main()
