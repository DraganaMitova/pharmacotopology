from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Sequence

from pharmacotopology.folding_contact_topology import ContactCandidate


@dataclass(frozen=True)
class EnergyLandscapePacket:
    step_count: int
    contact_counts: tuple[int, ...]
    energy_values: tuple[float, ...]
    trajectory_contact_gain_monotonicity: float
    energy_descent_consistency: float
    collapse_stability_score: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _chunked_counts(total_contacts: int, step_count: int) -> tuple[int, ...]:
    if step_count <= 1:
        return (total_contacts,)
    return tuple(
        round(total_contacts * step / (step_count - 1))
        for step in range(step_count)
    )


def build_energy_landscape(
    candidates: Sequence[ContactCandidate],
    *,
    step_count: int = 7,
) -> EnergyLandscapePacket:
    contact_counts = _chunked_counts(len(candidates), step_count)
    total_strength = sum(candidate.score for candidate in candidates)
    energy_values = []
    for count in contact_counts:
        active_strength = sum(candidate.score for candidate in candidates[:count])
        progress = active_strength / total_strength if total_strength else 0.0
        frustration = max(0.0, count - len(candidates) * 0.82) / max(len(candidates), 1)
        energy_values.append(round(1.0 - progress * 0.72 + frustration * 0.08, 6))

    contact_gain_steps = sum(
        1
        for left, right in zip(contact_counts, contact_counts[1:])
        if right >= left
    )
    energy_descent_steps = sum(
        1
        for left, right in zip(energy_values, energy_values[1:])
        if right <= left
    )
    denominator = max(step_count - 1, 1)
    return EnergyLandscapePacket(
        step_count=step_count,
        contact_counts=contact_counts,
        energy_values=tuple(energy_values),
        trajectory_contact_gain_monotonicity=_rounded(contact_gain_steps / denominator),
        energy_descent_consistency=_rounded(energy_descent_steps / denominator),
        collapse_stability_score=_rounded(
            (contact_gain_steps / denominator) * 0.45
            + (energy_descent_steps / denominator) * 0.45
            + (1.0 - energy_values[-1]) * 0.10
        ),
    )


def render_curve_svg(
    values: Sequence[float],
    *,
    title: str,
    y_label: str,
    width: int = 640,
    height: int = 220,
) -> str:
    margin = 34
    plot_width = width - margin * 2
    plot_height = height - margin * 2
    if not values:
        values = (0.0,)
    max_value = max(max(values), 1.0)
    min_value = min(min(values), 0.0)
    span = max(max_value - min_value, 1e-6)
    points = []
    for index, value in enumerate(values):
        x = margin + (plot_width * index / max(len(values) - 1, 1))
        y = margin + plot_height - ((value - min_value) / span) * plot_height
        points.append(f"{x:.2f},{y:.2f}")
    polyline = " ".join(points)
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{title}">
  <rect width="{width}" height="{height}" fill="#f8faf7"/>
  <text x="{margin}" y="22" font-family="Arial, sans-serif" font-size="15" fill="#22302d">{title}</text>
  <text x="{margin}" y="{height - 8}" font-family="Arial, sans-serif" font-size="11" fill="#58635e">closure step</text>
  <text x="8" y="{margin}" font-family="Arial, sans-serif" font-size="11" fill="#58635e">{y_label}</text>
  <line x1="{margin}" y1="{height - margin}" x2="{width - margin}" y2="{height - margin}" stroke="#9aa69f"/>
  <line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height - margin}" stroke="#9aa69f"/>
  <polyline points="{polyline}" fill="none" stroke="#25636a" stroke-width="3"/>
</svg>
"""
