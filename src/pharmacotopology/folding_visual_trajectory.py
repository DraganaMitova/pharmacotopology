from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_contact_topology import ContactCandidate
from pharmacotopology.folding_energy_landscape import EnergyLandscapePacket
from pharmacotopology.folding_native_contact_eval import (
    ContactPair,
    normalized_contact_pairs,
)


def _escape(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _scaled_pair(pair: ContactPair, length: int, size: int, margin: int) -> tuple[float, float]:
    span = max(length - 1, 1)
    x = margin + (pair[0] - 1) / span * (size - margin * 2)
    y = margin + (pair[1] - 1) / span * (size - margin * 2)
    return x, y


def render_contact_map_svg(
    *,
    row_id: str,
    sequence_length: int,
    contact_pairs: Iterable[Sequence[int]],
    title: str,
    color: str,
    size: int = 420,
) -> str:
    margin = 34
    pairs = normalized_contact_pairs(contact_pairs)
    dots = []
    for pair in pairs:
        x, y = _scaled_pair(pair, sequence_length, size, margin)
        dots.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.2" fill="{color}"/>')
        dots.append(f'<circle cx="{y:.2f}" cy="{x:.2f}" r="2.2" fill="{color}" opacity="0.58"/>')
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" viewBox="0 0 {size} {size}" role="img" aria-label="{title}">
  <rect width="{size}" height="{size}" fill="#f8faf7"/>
  <text x="{margin}" y="22" font-family="Arial, sans-serif" font-size="14" fill="#22302d">{_escape(title)}</text>
  <line x1="{margin}" y1="{margin}" x2="{margin}" y2="{size - margin}" stroke="#9aa69f"/>
  <line x1="{margin}" y1="{size - margin}" x2="{size - margin}" y2="{size - margin}" stroke="#9aa69f"/>
  <line x1="{margin}" y1="{margin}" x2="{size - margin}" y2="{size - margin}" stroke="#cfd8d2" stroke-dasharray="4 4"/>
  {''.join(dots)}
  <text x="{margin}" y="{size - 8}" font-family="Arial, sans-serif" font-size="11" fill="#58635e">row: {_escape(row_id)} | contacts: {len(pairs)}</text>
</svg>
"""


def render_contact_overlay_svg(
    *,
    row_id: str,
    sequence_length: int,
    native_pairs: Iterable[Sequence[int]],
    predicted_pairs: Iterable[Sequence[int]],
    size: int = 420,
) -> str:
    native = set(normalized_contact_pairs(native_pairs))
    predicted = set(normalized_contact_pairs(predicted_pairs))
    all_pairs = sorted(native | predicted)
    margin = 34
    dots = []
    for pair in all_pairs:
        x, y = _scaled_pair(pair, sequence_length, size, margin)
        if pair in native and pair in predicted:
            color = "#218a5a"
            radius = 2.8
        elif pair in native:
            color = "#294f9b"
            radius = 2.4
        else:
            color = "#c44b3a"
            radius = 2.0
        dots.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius}" fill="{color}" opacity="0.88"/>')
        dots.append(f'<circle cx="{y:.2f}" cy="{x:.2f}" r="{radius}" fill="{color}" opacity="0.46"/>')
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" viewBox="0 0 {size} {size}" role="img" aria-label="contact map overlay">
  <rect width="{size}" height="{size}" fill="#f8faf7"/>
  <text x="{margin}" y="22" font-family="Arial, sans-serif" font-size="14" fill="#22302d">Contact map overlay</text>
  <line x1="{margin}" y1="{margin}" x2="{margin}" y2="{size - margin}" stroke="#9aa69f"/>
  <line x1="{margin}" y1="{size - margin}" x2="{size - margin}" y2="{size - margin}" stroke="#9aa69f"/>
  <line x1="{margin}" y1="{margin}" x2="{size - margin}" y2="{size - margin}" stroke="#cfd8d2" stroke-dasharray="4 4"/>
  {''.join(dots)}
  <g font-family="Arial, sans-serif" font-size="11" fill="#34423d">
    <rect x="{margin}" y="{size - 48}" width="250" height="32" fill="#ffffff" stroke="#d5ddd5"/>
    <circle cx="{margin + 12}" cy="{size - 36}" r="4" fill="#218a5a"/><text x="{margin + 22}" y="{size - 32}">native + predicted</text>
    <circle cx="{margin + 122}" cy="{size - 36}" r="4" fill="#294f9b"/><text x="{margin + 132}" y="{size - 32}">native only</text>
    <circle cx="{margin + 202}" cy="{size - 36}" r="4" fill="#c44b3a"/><text x="{margin + 212}" y="{size - 32}">predicted only</text>
  </g>
  <text x="{margin}" y="{size - 8}" font-family="Arial, sans-serif" font-size="11" fill="#58635e">row: {_escape(row_id)}</text>
</svg>
"""


def _coarse_positions(length: int, contacts: Sequence[ContactPair], step: int, step_count: int) -> list[tuple[float, float]]:
    width = 760
    height = 240
    margin = 42
    progress = step / max(step_count - 1, 1)
    positions = []
    for index in range(length):
        t = index / max(length - 1, 1)
        x = margin + t * (width - margin * 2)
        wave = math.sin(t * math.pi * 2.0) * 20 * progress
        y = height / 2 + wave
        positions.append((x, y))
    for left, right in contacts:
        li = left - 1
        ri = right - 1
        if li >= len(positions) or ri >= len(positions):
            continue
        lx, ly = positions[li]
        rx, ry = positions[ri]
        center_x = (lx + rx) / 2
        center_y = (ly + ry) / 2
        pull = 0.12 * progress
        positions[li] = (lx + (center_x - lx) * pull, ly + (center_y - ly) * pull)
        positions[ri] = (rx + (center_x - rx) * pull, ry + (center_y - ry) * pull)
    return positions


def render_coarse_grain_svg(
    *,
    row_id: str,
    sequence_length: int,
    contacts: Sequence[ContactPair],
    title: str = "Coarse-grain final collapse",
    width: int = 760,
    height: int = 260,
) -> str:
    positions = _coarse_positions(sequence_length, contacts, step=6, step_count=7)
    step = max(1, sequence_length // 60)
    points = " ".join(f"{x:.2f},{y:.2f}" for x, y in positions[::step])
    contact_lines = []
    for left, right in contacts[:36]:
        lx, ly = positions[left - 1]
        rx, ry = positions[right - 1]
        contact_lines.append(
            f'<line x1="{lx:.2f}" y1="{ly:.2f}" x2="{rx:.2f}" y2="{ry:.2f}" stroke="#25636a" stroke-opacity="0.22" stroke-width="1.2"/>'
        )
    contact_markup = "".join(contact_lines)
    bead_step = max(1, sequence_length // 36)
    beads = []
    for index in range(0, sequence_length, bead_step):
        x, y = positions[index]
        beads.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="#22302d"/>')
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{title}">
  <rect width="{width}" height="{height}" fill="#f8faf7"/>
  <text x="34" y="24" font-family="Arial, sans-serif" font-size="15" fill="#22302d">{_escape(title)}</text>
{contact_markup}
  <polyline points="{points}" fill="none" stroke="#7d8f86" stroke-width="3"/>
  {''.join(beads)}
  <text x="34" y="{height - 12}" font-family="Arial, sans-serif" font-size="11" fill="#58635e">row: {_escape(row_id)} | coarse beads, not atomic structure</text>
</svg>
"""


def render_trajectory_html(
    *,
    row_id: str,
    sequence_length: int,
    candidates: Sequence[ContactCandidate],
    energy: EnergyLandscapePacket,
) -> str:
    frames = []
    for step, contact_count in enumerate(energy.contact_counts):
        contacts = tuple(candidate.pair() for candidate in candidates[:contact_count])
        frame_svg = render_coarse_grain_svg(
            row_id=row_id,
            sequence_length=sequence_length,
            contacts=contacts,
            title=f"Closure step {step}: contacts {contact_count}",
        )
        frames.append(f"<section><h2>Step {step}</h2>{frame_svg}</section>")
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Folding trajectory {_escape(row_id)}</title>
  <style>
    body {{ margin: 0; font-family: Arial, sans-serif; background: #f6f7f2; color: #1f2523; }}
    header {{ padding: 28px; background: #24302c; color: #f6f7f2; }}
    main {{ max-width: 920px; margin: 0 auto; padding: 20px; }}
    section {{ margin: 22px 0; padding: 16px; background: #fff; border: 1px solid #d5ddd5; border-radius: 6px; }}
    h1, h2 {{ margin: 0 0 10px; letter-spacing: 0; }}
  </style>
</head>
<body>
  <header>
    <h1>Coarse-Grained Folding Trajectory</h1>
    <p>This is a visual mechanism reconstruction, not an atomic folding pathway.</p>
  </header>
  <main>
    {''.join(frames)}
  </main>
</body>
</html>
"""


def write_visual_file(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path
