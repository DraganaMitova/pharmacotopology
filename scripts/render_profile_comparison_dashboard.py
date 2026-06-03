from __future__ import annotations

import argparse
from html import escape
from pathlib import Path
import sys
from typing import Any, Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.layer import (  # noqa: E402
    DEFAULT_MECHANISM_VECTORS,
    DEFAULT_NORMAL_BOUNDED_PROFILE,
    DEFAULT_TOPOLOGY_PROFILES,
    PATHOLOGY_DIMENSIONS,
    build_pharmacotopology_review,
)


DEFAULT_OUTPUT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/multi_profile_dashboard.html"
)


def _fmt(value: object) -> str:
    if isinstance(value, (int, float)):
        return f"{float(value):.3f}"
    return str(value)


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _score_class(value: float) -> str:
    if value >= 0.25:
        return "strong"
    if value > 0.0:
        return "positive"
    if value < 0.0:
        return "negative"
    return "neutral"


def _bar(value: float) -> str:
    bounded = _clamp(value, 0.0, 1.0)
    width = round(bounded * 150, 2)
    return (
        '<svg class="bar" viewBox="0 0 190 18" role="img" '
        f'aria-label="pressure {_fmt(value)}">'
        '<rect class="bar-track" x="0" y="4" width="150" height="10" rx="2" />'
        f'<rect class="bar-fill" x="0" y="4" width="{width}" height="10" rx="2" />'
        f'<text x="160" y="13">{escape(_fmt(value))}</text>'
        "</svg>"
    )


def _reviews() -> dict[str, dict[str, Any]]:
    return {
        profile_key: build_pharmacotopology_review(
            source=profile,
            target=DEFAULT_NORMAL_BOUNDED_PROFILE,
            mechanisms=DEFAULT_MECHANISM_VECTORS,
        )
        for profile_key, profile in DEFAULT_TOPOLOGY_PROFILES.items()
    }


def _ranking(review: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        row
        for row in review.get("Φ.ranking", [])
        if isinstance(row, dict)
    ]


def _result_lookup(review: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        str(result["mechanism_id"]): result
        for result in review.get("Φ.results", [])
        if isinstance(result, dict)
    }


def _profile_maps() -> str:
    cards = []
    for profile_key, profile in DEFAULT_TOPOLOGY_PROFILES.items():
        rows = []
        for dimension in PATHOLOGY_DIMENSIONS:
            value = float(profile.dimensions[dimension])
            rows.append(
                "<div class=\"profile-row\">"
                f"<span>{escape(dimension)}</span>"
                f"{_bar(value)}"
                "</div>"
            )
        cards.append(
            "<article class=\"profile-card\">"
            f"<h3>{escape(profile_key)}</h3>"
            f"<p>{escape(profile.description)}</p>"
            f"<div class=\"profile-rows\">{''.join(rows)}</div>"
            "</article>"
        )
    return f"""
    <section>
      <div class="section-heading">
        <h2>Baseline Topology Maps</h2>
        <p>Synthetic pressure maps only. They are not diagnoses, patient models, or medication targets.</p>
      </div>
      <div class="profile-grid">{''.join(cards)}</div>
    </section>
    """


def _combined_ranking_table(reviews: Mapping[str, Mapping[str, Any]]) -> str:
    profile_keys = tuple(reviews)
    header = "".join(
        f"<th>{escape(profile_key)}<br><span>rank / net</span></th>"
        for profile_key in profile_keys
    )
    mechanism_ids = tuple(vector.mechanism_id for vector in DEFAULT_MECHANISM_VECTORS)
    ranking_lookup = {
        profile_key: {
            str(item["mechanism_id"]): item
            for item in _ranking(review)
        }
        for profile_key, review in reviews.items()
    }
    rows = []
    for mechanism_id in mechanism_ids:
        cells = []
        for profile_key in profile_keys:
            item = ranking_lookup[profile_key][mechanism_id]
            score = float(item["net_topology_health_score"])
            cells.append(
                f'<td class="{_score_class(score)}">'
                f"#{int(item['rank'])} / {_fmt(score)}</td>"
            )
        rows.append(
            f"<tr><th>{escape(mechanism_id)}</th>{''.join(cells)}</tr>"
        )
    return f"""
    <section>
      <div class="section-heading">
        <h2>Combined Ranking Table</h2>
        <p>Each cell shows rank and net topology health score for the same mechanism under each synthetic profile.</p>
      </div>
      <div class="table-scroll">
        <table>
          <thead><tr><th>mechanism_id</th>{header}</tr></thead>
          <tbody>{''.join(rows)}</tbody>
        </table>
      </div>
    </section>
    """


def _mechanism_heatmap(reviews: Mapping[str, Mapping[str, Any]]) -> str:
    profile_keys = tuple(reviews)
    header = "".join(f"<th>{escape(profile_key)}</th>" for profile_key in profile_keys)
    rows = []
    for vector in DEFAULT_MECHANISM_VECTORS:
        cells = []
        for profile_key in profile_keys:
            ranking = {
                str(item["mechanism_id"]): item
                for item in _ranking(reviews[profile_key])
            }
            score = float(ranking[vector.mechanism_id]["net_topology_health_score"])
            cells.append(
                f'<td class="{_score_class(score)}"><span>{escape(_fmt(score))}</span></td>'
            )
        rows.append(
            f"<tr><th>{escape(vector.mechanism_id)}</th>{''.join(cells)}</tr>"
        )
    return f"""
    <section>
      <div class="section-heading">
        <h2>Mechanism/Profile Heatmap</h2>
        <p>Net score by mechanism and synthetic profile. Positive cells are modeling signals only.</p>
      </div>
      <div class="table-scroll">
        <table class="heatmap">
          <thead><tr><th>mechanism_id</th>{header}</tr></thead>
          <tbody>{''.join(rows)}</tbody>
        </table>
      </div>
    </section>
    """


def _summary(reviews: Mapping[str, Mapping[str, Any]]) -> str:
    summaries = []
    for vector in DEFAULT_MECHANISM_VECTORS:
        scores = []
        ranks = []
        top_profiles = []
        for profile_key, review in reviews.items():
            ranking = {
                str(item["mechanism_id"]): item
                for item in _ranking(review)
            }
            item = ranking[vector.mechanism_id]
            scores.append(float(item["net_topology_health_score"]))
            ranks.append(int(item["rank"]))
            if int(item["rank"]) <= 3:
                top_profiles.append(profile_key)
        mean_score = sum(scores) / len(scores)
        score_range = max(scores) - min(scores)
        rank_range = max(ranks) - min(ranks)
        summaries.append(
            {
                "mechanism_id": vector.mechanism_id,
                "mean_score": mean_score,
                "min_score": min(scores),
                "score_range": score_range,
                "rank_range": rank_range,
                "top_profiles": tuple(top_profiles),
            }
        )

    consistently_well = sorted(
        [
            item
            for item in summaries
            if item["min_score"] > 0.0 and item["rank_range"] <= 4
        ],
        key=lambda item: (float(item["mean_score"]), -float(item["score_range"])),
        reverse=True,
    )[:5]
    profile_specific = sorted(
        summaries,
        key=lambda item: (int(item["rank_range"]), float(item["score_range"])),
        reverse=True,
    )[:5]

    def _items(items: list[dict[str, Any]]) -> str:
        return "".join(
            "<li>"
            f"<strong>{escape(str(item['mechanism_id']))}</strong> "
            f"mean {_fmt(item['mean_score'])}, range {_fmt(item['score_range'])}, "
            f"rank range {int(item['rank_range'])}, "
            f"top profiles {escape(', '.join(item['top_profiles']) or 'none')}"
            "</li>"
            for item in items
        )

    return f"""
    <section>
      <div class="section-heading">
        <h2>Cross-Profile Summary</h2>
        <p>Broadly stabilizing means consistently positive under this synthetic model. Profile-specific means rank or score changes sharply.</p>
      </div>
      <div class="summary-grid">
        <article>
          <h3>Consistently Positive</h3>
          <ul>{_items(consistently_well)}</ul>
        </article>
        <article>
          <h3>Highly Profile-Specific</h3>
          <ul>{_items(profile_specific)}</ul>
        </article>
      </div>
    </section>
    """


def _stylesheet() -> str:
    return """
    :root {
      --ink: #202124;
      --muted: #5f6368;
      --line: #d9dce1;
      --panel: #ffffff;
      --page: #f7f8fa;
      --strong: #cfe9dd;
      --positive: #e2f0ea;
      --negative: #f5e4df;
      --neutral: #eceff3;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      color: var(--ink);
      background: var(--page);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.45;
    }
    main { width: min(1240px, calc(100vw - 32px)); margin: 0 auto 48px; }
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
    h3 { font-size: 15px; margin-bottom: 6px; }
    section { margin: 28px 0; }
    .kicker, .layer-line, .section-heading p, .profile-card p {
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
      border: 1px solid #a85648;
      background: #fff8f5;
      color: #6f2f24;
      padding: 10px 12px;
      font-weight: 800;
      font-size: 13px;
    }
    .profile-grid, .summary-grid {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 10px;
    }
    .profile-card, .summary-grid article {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 12px;
    }
    .profile-row {
      display: grid;
      grid-template-columns: minmax(190px, 1fr) 190px;
      gap: 10px;
      align-items: center;
      padding: 4px 0;
      font-size: 13px;
    }
    .bar-track { fill: #e4e7eb; }
    .bar-fill { fill: #52616f; }
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
    th span { color: var(--muted); font-weight: 500; }
    td.strong { background: var(--strong); color: #245c45; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
    td.positive { background: var(--positive); color: #2f604c; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
    td.negative { background: var(--negative); color: #78392f; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
    td.neutral { background: var(--neutral); color: #4d5156; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
    .summary-grid ul { margin: 0; padding-left: 18px; }
    .summary-grid li { margin-bottom: 8px; }
    footer {
      margin: 36px 0 0;
      padding: 18px;
      border: 1px solid var(--line);
      background: #fff;
      font-weight: 700;
    }
    footer p { margin: 4px 0; }
    @media (max-width: 880px) {
      main { width: min(100vw - 20px, 1240px); }
      .hero { display: block; }
      .safety-badge { display: inline-block; margin-top: 14px; }
      .profile-grid, .summary-grid { grid-template-columns: 1fr; }
      .profile-row { grid-template-columns: 1fr; }
    }
    """


def render_profile_comparison_dashboard_html() -> str:
    reviews = _reviews()
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>KNOT Pharmacotopology Multi-Profile Dashboard</title>
  <style>{_stylesheet()}</style>
</head>
<body>
  <main>
    <header class="hero">
      <div>
        <p class="kicker">KNOT topology simulation</p>
        <h1>Multi-Profile Pharmacotopology Dashboard</h1>
        <p class="layer-line">Same mechanism vectors across synthetic topology profiles</p>
      </div>
      <div class="safety-badge">SIMULATION ONLY &mdash; NOT MEDICAL ADVICE</div>
    </header>
    {_profile_maps()}
    {_combined_ranking_table(reviews)}
    {_mechanism_heatmap(reviews)}
    {_summary(reviews)}
    <footer>
      <p>Φ.review does not recommend medication.</p>
      <p>Φ.review does not infer treatment.</p>
      <p>Profile comparison is synthetic ranking review, not clinical comparison.</p>
      <p>Clinical claims are denied unless externally validated.</p>
    </footer>
  </main>
</body>
</html>
"""


def render_profile_comparison_dashboard(
    output_path: Path = DEFAULT_OUTPUT_PATH,
) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_profile_comparison_dashboard_html(), encoding="utf-8")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render a static multi-profile pharmacotopology dashboard."
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_PATH),
        help="Path for the generated multi-profile dashboard HTML file.",
    )
    args = parser.parse_args()
    output_path = render_profile_comparison_dashboard(Path(args.output))
    print(output_path)


if __name__ == "__main__":
    main()
