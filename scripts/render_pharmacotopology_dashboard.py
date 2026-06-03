from __future__ import annotations

import argparse
import json
import sys
from html import escape
from pathlib import Path
from typing import Any, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.layer import (  # noqa: E402
    DEFAULT_MECHANISM_VECTORS,
    DEFAULT_NORMAL_BOUNDED_PROFILE,
    DEFAULT_TOPOLOGY_PROFILES,
    build_pharmacotopology_review,
)


DEFAULT_INPUT_PATH = Path("first_contact_clean_pharmacotopology_layer_run/memory.jsonl")
DEFAULT_OUTPUT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/pharmacotopology_dashboard.html"
)

TOPOLOGY_DIMENSIONS: tuple[str, ...] = (
    "salience_amplification",
    "recurrence_overbinding",
    "symbolic_closure_pressure",
    "threat_propagation",
    "falsification_weakness",
    "boundary_instability",
    "agency_confusion",
    "sensory_intrusion",
    "cognitive_fragmentation",
    "negative_shutdown",
    "sleep_instability",
    "collapse_cost",
)

REQUIRED_REVIEW_KEY = "Φ.review"


def _json_rows(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parsed = json.loads(line)
        if isinstance(parsed, dict):
            rows.append(parsed)
    return rows


def read_phi_packet(memory_path: Path) -> dict[str, Any]:
    for row in _json_rows(memory_path):
        content = row.get("content")
        if not isinstance(content, str):
            continue
        packet = json.loads(content)
        if isinstance(packet, dict) and REQUIRED_REVIEW_KEY in packet:
            return packet
    raise ValueError(f"No {REQUIRED_REVIEW_KEY} packet found in {memory_path}")


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _fmt(value: object) -> str:
    if isinstance(value, (int, float)):
        return f"{float(value):.3f}"
    return str(value)


def _cell_class(delta: float) -> str:
    if delta < 0:
        return "reduce"
    if delta > 0:
        return "worsen"
    return "neutral"


def _cell_label(delta: float) -> str:
    if delta < 0:
        return f"reduce {_fmt(delta)}"
    if delta > 0:
        return f"worsen +{_fmt(delta)}"
    return "neutral 0.000"


def _topology_bar_svg(value: float) -> str:
    bounded = _clamp(value, 0.0, 1.0)
    width = round(bounded * 220, 2)
    return (
        '<svg class="topology-bar" viewBox="0 0 260 18" role="img" '
        f'aria-label="topology value {_fmt(value)}">'
        '<rect class="bar-track" x="0" y="4" width="220" height="10" rx="2" />'
        f'<rect class="bar-fill" x="0" y="4" width="{width}" height="10" rx="2" />'
        f'<text x="232" y="13">{escape(_fmt(value))}</text>'
        "</svg>"
    )


def _score_bar_svg(value: float, *, signed: bool, css_class: str) -> str:
    if signed:
        bounded = _clamp(value, -1.0, 1.0)
        midpoint = 110.0
        span = 100.0
        bar_width = abs(bounded) * span
        x = midpoint if bounded >= 0 else midpoint - bar_width
        return (
            '<svg class="score-bar signed" viewBox="0 0 260 18" role="img" '
            f'aria-label="score {_fmt(value)}">'
            '<rect class="bar-track" x="10" y="4" width="200" height="10" rx="2" />'
            '<line class="zero-line" x1="110" y1="2" x2="110" y2="16" />'
            f'<rect class="{css_class}" x="{x:.2f}" y="4" '
            f'width="{bar_width:.2f}" height="10" rx="2" />'
            f'<text x="224" y="13">{escape(_fmt(value))}</text>'
            "</svg>"
        )

    bounded = _clamp(value, 0.0, 1.0)
    bar_width = bounded * 200.0
    return (
        '<svg class="score-bar" viewBox="0 0 260 18" role="img" '
        f'aria-label="score {_fmt(value)}">'
        '<rect class="bar-track" x="10" y="4" width="200" height="10" rx="2" />'
        f'<rect class="{css_class}" x="10" y="4" '
        f'width="{bar_width:.2f}" height="10" rx="2" />'
        f'<text x="224" y="13">{escape(_fmt(value))}</text>'
        "</svg>"
    )


def _source_dimensions(review: Mapping[str, Any]) -> dict[str, float]:
    source = review["Φ.source_profile"]
    dimensions = dict(source["dimensions"])
    dimensions["collapse_cost"] = 0.0
    return {key: float(dimensions.get(key, 0.0)) for key in TOPOLOGY_DIMENSIONS}


def _mechanism_lookup(review: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        str(vector["mechanism_id"]): vector
        for vector in review.get("Φ.mechanism_vectors", [])
        if isinstance(vector, dict)
    }


def _result_lookup(review: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        str(result["mechanism_id"]): result
        for result in review.get("Φ.results", [])
        if isinstance(result, dict)
    }


def _ranking(review: Mapping[str, Any]) -> list[dict[str, Any]]:
    ranking = review.get("Φ.ranking", [])
    return [row for row in ranking if isinstance(row, dict)]


def _html_header(packet: Mapping[str, Any]) -> str:
    schema = escape(str(packet.get("κ", "κ.unavailable")))
    layer = escape(str(packet.get("χ", "χ.clean.pharmacotopology_layer")))
    return f"""
    <header class="hero">
      <div>
        <p class="kicker">KNOT topology simulation</p>
        <h1>KNOT Φ.review Pharmacotopology Dashboard</h1>
        <p class="layer-line">Layer: {layer} · Schema: {schema}</p>
      </div>
      <div class="safety-badge">SIMULATION ONLY &mdash; NOT MEDICAL ADVICE</div>
    </header>
    """


def _topology_state_map(review: Mapping[str, Any]) -> str:
    dimensions = _source_dimensions(review)
    rows = []
    for dimension in TOPOLOGY_DIMENSIONS:
        value = dimensions[dimension]
        rows.append(
            "<div class=\"topology-row\">"
            f"<span>{escape(dimension)}</span>"
            f"{_topology_bar_svg(value)}"
            "</div>"
        )
    return f"""
    <section>
      <div class="section-heading">
        <h2>Topology State Map</h2>
        <p>Baseline pressure profile. Collapse cost starts at zero until a perturbation is applied.</p>
      </div>
      <div class="topology-grid">
        {''.join(rows)}
      </div>
    </section>
    """


def _ranking_chart(review: Mapping[str, Any]) -> str:
    results = _result_lookup(review)
    rows = []
    for item in _ranking(review):
        mechanism_key = str(item["mechanism_id"])
        mechanism_id = escape(mechanism_key)
        result = results.get(mechanism_key, {})
        pathology = float(item["pathology_reduction_score"])
        collapse = float(item["collapse_cost_score"])
        net = float(item["net_topology_health_score"])
        evidence_weight = float(
            item.get("evidence_weight", result.get("evidence_weight", 0.0))
        )
        uncertainty_radius = float(
            item.get("uncertainty_radius", result.get("uncertainty_radius", 0.0))
        )
        evidence_label = escape(
            str(
                item.get(
                    "evidence_readiness_label",
                    result.get("evidence_readiness_label", "unknown"),
                )
            )
        )
        net_interval = result.get("net_topology_health_interval", {})
        net_lower = float(net_interval.get("lower", net))
        net_upper = float(net_interval.get("upper", net))
        fit_label = escape(str(item.get("fit_label", "")))
        rows.append(
            "<article class=\"ranking-card\">"
            f"<h3>#{int(item['rank'])} {mechanism_id}</h3>"
            f"<p class=\"fit-label\">{fit_label}</p>"
            f"<p class=\"evidence-line\">{evidence_label} · "
            f"evidence_weight {_fmt(evidence_weight)} · "
            f"uncertainty_radius {_fmt(uncertainty_radius)} · "
            f"net interval {_fmt(net_lower)} to {_fmt(net_upper)}</p>"
            "<div class=\"score-row\"><span>pathology_reduction_score <em>topology distance change, not cure</em></span>"
            f"{_score_bar_svg(pathology, signed=True, css_class='score-pathology')}</div>"
            "<div class=\"score-row\"><span>collapse_cost_score</span>"
            f"{_score_bar_svg(collapse, signed=False, css_class='score-collapse')}</div>"
            "<div class=\"score-row\"><span>net_topology_health_score</span>"
            f"{_score_bar_svg(net, signed=True, css_class='score-net')}</div>"
            "</article>"
        )
    return f"""
    <section>
      <div class="section-heading">
        <h2>Mechanism Ranking Chart</h2>
        <p>Ranking compares topology perturbation scores. Collapse cost stays visible beside every reduction score.</p>
      </div>
      <div class="ranking-list">
        {''.join(rows)}
      </div>
    </section>
    """


def _calibration_readiness(review: Mapping[str, Any]) -> str:
    calibration = review.get("Φ.calibration_readiness", {})
    status = escape(str(calibration.get("calibration_status", "unknown")))
    practical_use = escape(str(calibration.get("practical_use", "unknown")))
    clinical_allowed = str(calibration.get("clinical_use_allowed", False)).lower()
    evidence_count = int(calibration.get("evidence_backed_vectors", 0))
    vector_count = int(calibration.get("mechanism_vectors_reviewed", 0))
    mean_weight = float(calibration.get("mean_evidence_weight", 0.0))
    mean_uncertainty = float(calibration.get("mean_uncertainty_radius", 0.0))
    blockers = [
        str(item) for item in calibration.get("blockers", []) if isinstance(item, str)
    ]
    next_steps = [
        str(item)
        for item in calibration.get("next_steps", [])
        if isinstance(item, str)
    ]
    blocker_items = "".join(f"<li>{escape(item)}</li>" for item in blockers)
    next_step_items = "".join(f"<li>{escape(item)}</li>" for item in next_steps)
    return f"""
    <section>
      <div class="section-heading">
        <h2>Calibration Readiness</h2>
        <p>Practical status for bounded hypothesis comparison. Clinical use remains closed.</p>
      </div>
      <div class="readiness-grid">
        <div class="readiness-panel">
          <h3>{status}</h3>
          <p>practical_use: {practical_use}</p>
          <p>clinical_use_allowed: {escape(clinical_allowed)}</p>
          <p>evidence_backed_vectors: {evidence_count} / {vector_count}</p>
          <p>mean_evidence_weight: {_fmt(mean_weight)}</p>
          <p>mean_uncertainty_radius: {_fmt(mean_uncertainty)}</p>
        </div>
        <div class="readiness-panel">
          <h3>Blockers</h3>
          <ul>{blocker_items}</ul>
        </div>
        <div class="readiness-panel">
          <h3>Next Calibration Steps</h3>
          <ul>{next_step_items}</ul>
        </div>
      </div>
    </section>
    """


def _calibration_vector_table(review: Mapping[str, Any]) -> str:
    results = _result_lookup(review)
    rows = []
    for item in _ranking(review):
        mechanism_id = str(item["mechanism_id"])
        result = results.get(mechanism_id, {})
        net_interval = result.get("net_topology_health_interval", {})
        pathology_interval = result.get("pathology_reduction_interval", {})
        collapse_interval = result.get("collapse_cost_interval", {})
        primary_sources = result.get("primary_evidence_sources", ())
        evidence_refs = result.get("evidence_refs", ())
        calibration_blockers = result.get("calibration_blockers", ())
        source_label = "; ".join(primary_sources) if primary_sources else "none_attached"
        refs_label = "; ".join(evidence_refs) if evidence_refs else "none_attached"
        blockers_label = (
            "; ".join(calibration_blockers)
            if calibration_blockers
            else "none_attached"
        )
        rows.append(
            "<tr>"
            f"<th>{escape(mechanism_id)}</th>"
            f"<td>{escape(str(result.get('abstract_compound_class', 'unknown')))}</td>"
            f"<td>{escape(str(result.get('protein_family', 'unknown')))}</td>"
            f"<td>{escape(str(result.get('protein_mechanism_class', 'unknown')))}</td>"
            f"<td>{escape(str(result.get('protein_state_shift', 'unknown')))}</td>"
            f"<td>{escape(str(result.get('pathway_network_perturbation', 'unknown')))}</td>"
            f"<td>{escape(str(result.get('evidence_stage', 'unknown')))}</td>"
            f"<td>{escape(_fmt(item.get('evidence_weight', 0.0)))}</td>"
            f"<td>{escape(_fmt(item.get('uncertainty_radius', 0.0)))}</td>"
            f"<td>{escape(str(item.get('evidence_readiness_label', 'unknown')))}</td>"
            f"<td>{escape(source_label)}</td>"
            f"<td>{escape(refs_label)}</td>"
            f"<td>{escape(blockers_label)}</td>"
            f"<td>{escape(str(result.get('confidence_interval_kind', 'unknown')))}</td>"
            f"<td>{escape(_fmt(pathology_interval.get('lower', 0.0)))} to "
            f"{escape(_fmt(pathology_interval.get('upper', 0.0)))}</td>"
            f"<td>{escape(_fmt(collapse_interval.get('lower', 0.0)))} to "
            f"{escape(_fmt(collapse_interval.get('upper', 0.0)))}</td>"
            f"<td>{escape(_fmt(net_interval.get('lower', 0.0)))} to "
            f"{escape(_fmt(net_interval.get('upper', 0.0)))}</td>"
            "</tr>"
        )
    return f"""
    <section>
      <div class="section-heading">
        <h2>Per-Vector Calibration Table</h2>
        <p>Evidence posture and uncertainty interval for every mechanism vector.</p>
      </div>
      <div class="table-scroll">
        <table class="readiness-table">
          <thead>
            <tr>
              <th>mechanism_id</th>
              <th>abstract compound class</th>
              <th>protein family</th>
              <th>protein mechanism class</th>
              <th>protein state shift</th>
              <th>pathway/network perturbation</th>
              <th>evidence_stage</th>
              <th>evidence_weight</th>
              <th>uncertainty_radius</th>
              <th>readiness</th>
              <th>primary evidence sources</th>
              <th>evidence refs</th>
              <th>calibration blockers</th>
              <th>interval kind</th>
              <th>pathology interval</th>
              <th>collapse interval</th>
              <th>net interval</th>
            </tr>
          </thead>
          <tbody>{''.join(rows)}</tbody>
        </table>
      </div>
    </section>
    """


def _profile_comparison() -> str:
    rows = []
    for profile_key, profile in DEFAULT_TOPOLOGY_PROFILES.items():
        comparison_review = build_pharmacotopology_review(
            source=profile,
            target=DEFAULT_NORMAL_BOUNDED_PROFILE,
            mechanisms=DEFAULT_MECHANISM_VECTORS,
        )
        top = comparison_review["Φ.ranking"][0]
        net_interval = top.get("net_topology_health_interval", {})
        rows.append(
            "<tr>"
            f"<th>{escape(profile_key)}</th>"
            f"<td>{escape(str(profile.profile_id))}</td>"
            f"<td>{escape(str(top['mechanism_id']))}</td>"
            f"<td>{escape(str(top['fit_label']))}</td>"
            f"<td>{escape(_fmt(top['net_topology_health_score']))}</td>"
            f"<td>{escape(_fmt(net_interval.get('lower', 0.0)))} to "
            f"{escape(_fmt(net_interval.get('upper', 0.0)))}</td>"
            f"<td>{int(comparison_review['Φ.review']['destabilizing_mechanism_count'])}</td>"
            "</tr>"
        )
    return f"""
    <section>
      <div class="section-heading">
        <h2>Synthetic Profile Comparison</h2>
        <p>Same mechanism set, different synthetic pressure maps. This is a ranking robustness view, not a clinical comparison.</p>
      </div>
      <div class="table-scroll">
        <table class="profile-comparison">
          <thead>
            <tr>
              <th>profile key</th>
              <th>profile_id</th>
              <th>top mechanism</th>
              <th>top fit label</th>
              <th>top net score</th>
              <th>top net interval</th>
              <th>destabilizing count</th>
            </tr>
          </thead>
          <tbody>{''.join(rows)}</tbody>
        </table>
      </div>
    </section>
    """


def _heatmap(review: Mapping[str, Any]) -> str:
    mechanisms = _mechanism_lookup(review)
    header_cells = "".join(f"<th>{escape(dimension)}</th>" for dimension in TOPOLOGY_DIMENSIONS)
    rows = []
    for item in _ranking(review):
        mechanism_id = str(item["mechanism_id"])
        vector = mechanisms[mechanism_id]
        deltas = vector.get("deltas", {})
        cells = []
        for dimension in TOPOLOGY_DIMENSIONS:
            if dimension == "collapse_cost":
                delta = float(vector.get("collapse_cost", 0.0))
            else:
                delta = float(deltas.get(dimension, 0.0))
            cells.append(
                f'<td class="{_cell_class(delta)}">'
                f'<span>{escape(_cell_label(delta))}</span></td>'
            )
        rows.append(
            f"<tr><th>{escape(mechanism_id)}</th>{''.join(cells)}</tr>"
        )
    return f"""
    <section>
      <div class="section-heading">
        <h2>Mechanism Heatmap</h2>
        <p>Cell labels show perturbation direction: reduce, worsen, or neutral.</p>
      </div>
      <div class="table-scroll">
        <table class="heatmap">
          <thead><tr><th>mechanism_id</th>{header_cells}</tr></thead>
          <tbody>{''.join(rows)}</tbody>
        </table>
      </div>
    </section>
    """


def _before_after(review: Mapping[str, Any]) -> str:
    baseline = _source_dimensions(review)
    results = _result_lookup(review)
    sections = []
    for item in _ranking(review):
        mechanism_id = str(item["mechanism_id"])
        result = results[mechanism_id]
        state = dict(result.get("resulting_state", {}))
        deltas = dict(result.get("topology_delta", {}))
        state["collapse_cost"] = float(result["collapse_cost_score"])
        deltas["collapse_cost"] = float(result["collapse_cost_score"])
        rows = []
        for dimension in TOPOLOGY_DIMENSIONS:
            baseline_value = float(baseline[dimension])
            result_value = float(state.get(dimension, baseline_value))
            delta = float(deltas.get(dimension, result_value - baseline_value))
            rows.append(
                "<tr>"
                f"<th>{escape(dimension)}</th>"
                f"<td>{escape(_fmt(baseline_value))}</td>"
                f"<td>{escape(_fmt(result_value))}</td>"
                f'<td class="{_cell_class(delta)}">{escape(_cell_label(delta))}</td>'
                "</tr>"
            )
        sections.append(
            f"""
            <details class="comparison-card">
              <summary>{escape(mechanism_id)}</summary>
              <table>
                <thead>
                  <tr>
                    <th>dimension</th>
                    <th>baseline topology value</th>
                    <th>resulting topology value</th>
                    <th>delta</th>
                  </tr>
                </thead>
                <tbody>{''.join(rows)}</tbody>
              </table>
            </details>
            """
        )
    return f"""
    <section>
      <div class="section-heading">
        <h2>Before/After Comparison</h2>
        <p>Did this perturbation restore bounded review, or did it just suppress/collapse the system?</p>
      </div>
      <div class="comparison-list">
        {''.join(sections)}
      </div>
    </section>
    """


def _safety_footer() -> str:
    return """
    <footer class="safety-footer">
      <p>Φ.review does not recommend medication.</p>
      <p>Φ.review does not infer treatment.</p>
      <p>Φ.review only simulates bounded topology perturbations.</p>
      <p>Clinical claims are denied unless externally validated.</p>
    </footer>
    """


def _stylesheet() -> str:
    return """
    :root {
      --ink: #202124;
      --muted: #5f6368;
      --line: #d9dce1;
      --panel: #ffffff;
      --page: #f7f8fa;
      --reduce: #287878;
      --reduce-bg: #dcefed;
      --worsen: #a85648;
      --worsen-bg: #f5e4df;
      --neutral: #70757a;
      --neutral-bg: #eceff3;
      --pathology: #5a6d9f;
      --collapse: #9d6b45;
      --net: #3f7d5c;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      color: var(--ink);
      background: var(--page);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.45;
    }
    main { width: min(1180px, calc(100vw - 32px)); margin: 0 auto 48px; }
    .hero {
      display: flex;
      justify-content: space-between;
      align-items: flex-end;
      gap: 24px;
      padding: 32px 0 24px;
      border-bottom: 1px solid var(--line);
      margin-bottom: 28px;
    }
    h1, h2, h3, p { margin-top: 0; }
    h1 { font-size: 30px; line-height: 1.15; margin-bottom: 8px; }
    h2 { font-size: 20px; margin-bottom: 4px; }
    h3 { font-size: 15px; margin-bottom: 4px; }
    section { margin: 28px 0; }
    .kicker, .layer-line, .section-heading p, .fit-label, .evidence-line, .score-row em {
      color: var(--muted);
    }
    .kicker {
      text-transform: uppercase;
      font-size: 12px;
      letter-spacing: 0;
      margin-bottom: 8px;
      font-weight: 700;
    }
    .safety-badge {
      flex: 0 0 auto;
      border: 1px solid var(--worsen);
      background: #fff8f5;
      color: #6f2f24;
      padding: 10px 12px;
      font-weight: 800;
      font-size: 13px;
    }
    .topology-grid, .ranking-list, .comparison-list {
      display: grid;
      gap: 10px;
    }
    .topology-row, .ranking-card, .comparison-card {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 6px;
    }
    .topology-row {
      display: grid;
      grid-template-columns: minmax(210px, 1fr) 260px;
      align-items: center;
      gap: 16px;
      padding: 10px 12px;
    }
    .ranking-card { padding: 12px; }
    .evidence-line {
      margin-bottom: 10px;
      font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
      font-size: 12px;
    }
    .score-row {
      display: grid;
      grid-template-columns: minmax(260px, 1fr) 260px;
      gap: 16px;
      align-items: center;
      padding: 5px 0;
    }
    .score-row span { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 13px; }
    .score-row em { display: block; font-family: inherit; font-style: normal; font-size: 12px; }
    .bar-track { fill: #e4e7eb; }
    .bar-fill { fill: #52616f; }
    .zero-line { stroke: #6f747b; stroke-width: 1; }
    .score-pathology { fill: var(--pathology); }
    .score-collapse { fill: var(--collapse); }
    .score-net { fill: var(--net); }
    svg text { fill: var(--ink); font-size: 12px; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
    .table-scroll { overflow-x: auto; border: 1px solid var(--line); background: var(--panel); }
    table { width: 100%; border-collapse: collapse; font-size: 13px; }
    th, td {
      border-bottom: 1px solid var(--line);
      border-right: 1px solid var(--line);
      padding: 8px;
      text-align: left;
      vertical-align: top;
    }
    th { background: #f1f3f5; font-weight: 700; }
    td.reduce, td.worsen, td.neutral { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
    .reduce { background: var(--reduce-bg); color: #175b5b; }
    .worsen { background: var(--worsen-bg); color: #78392f; }
    .neutral { background: var(--neutral-bg); color: #4d5156; }
    .heatmap th:first-child { position: sticky; left: 0; z-index: 1; }
    .comparison-card { padding: 0; overflow: hidden; }
    .readiness-grid {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 10px;
    }
    .readiness-panel {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 12px;
    }
    .readiness-panel p {
      margin-bottom: 6px;
      font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
      font-size: 12px;
    }
    .readiness-panel ul {
      margin: 0;
      padding-left: 18px;
      font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
      font-size: 12px;
    }
    .readiness-panel li { margin-bottom: 6px; }
    .comparison-card summary {
      cursor: pointer;
      padding: 12px;
      font-weight: 700;
      font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
    }
    .comparison-card table { border-top: 1px solid var(--line); }
    .safety-footer {
      margin: 36px 0 0;
      padding: 18px;
      border: 1px solid var(--line);
      background: #fff;
      font-weight: 700;
    }
    .safety-footer p { margin: 4px 0; }
    @media (max-width: 760px) {
      main { width: min(100vw - 20px, 1180px); }
      .hero { display: block; }
      .safety-badge { display: inline-block; margin-top: 14px; }
      .topology-row, .score-row { grid-template-columns: 1fr; }
      .readiness-grid { grid-template-columns: 1fr; }
    }
    """


def render_dashboard_html(packet: Mapping[str, Any]) -> str:
    review = packet[REQUIRED_REVIEW_KEY]
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>KNOT Φ.review Pharmacotopology Dashboard</title>
  <style>{_stylesheet()}</style>
</head>
<body>
  <main>
    {_html_header(packet)}
    {_topology_state_map(review)}
    {_ranking_chart(review)}
    {_calibration_readiness(review)}
    {_calibration_vector_table(review)}
    {_profile_comparison()}
    {_heatmap(review)}
    {_before_after(review)}
    {_safety_footer()}
  </main>
</body>
</html>
"""


def render_dashboard(
    input_path: Path = DEFAULT_INPUT_PATH,
    output_path: Path = DEFAULT_OUTPUT_PATH,
) -> Path:
    packet = read_phi_packet(input_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_dashboard_html(packet), encoding="utf-8")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render a static visual-only KNOT Φ.review dashboard."
    )
    parser.add_argument(
        "--input",
        default=str(DEFAULT_INPUT_PATH),
        help="Path to the pharmacotopology memory.jsonl file.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_PATH),
        help="Path for the generated dashboard HTML file.",
    )
    args = parser.parse_args()

    output_path = render_dashboard(Path(args.input), Path(args.output))
    print(output_path)


if __name__ == "__main__":
    main()
