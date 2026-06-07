from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Mapping, Sequence


PHYSICAL_TERM_ABLATION_KIND = "active_physical_term_ablation_v1"
ABLATION_NAMES = (
    "without_loop_strain",
    "without_steric_clash",
    "without_burial_gain",
    "without_unsatisfied_polar_penalty",
    "without_future_frustration",
    "without_decoy_margin",
)


@dataclass(frozen=True)
class PhysicalTermAblationRow:
    ablation_name: str
    false_nucleus_rate: float
    contact_cluster_precision: float
    long_range_contact_recall: float
    real_vs_decoy_physical_enrichment_ratio: float
    composite_selection_score: float
    delta_vs_full_selector: float
    term_interpretation: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def composite_selection_score(
    *,
    false_nucleus_rate: float,
    contact_cluster_precision: float,
    long_range_contact_recall: float,
    real_vs_decoy_physical_enrichment_ratio: float,
) -> float:
    return round(
        real_vs_decoy_physical_enrichment_ratio
        + contact_cluster_precision
        + long_range_contact_recall
        - false_nucleus_rate,
        6,
    )


def classify_ablation_rows(
    rows: Sequence[PhysicalTermAblationRow],
) -> dict[str, object]:
    if not rows:
        return {
            "best_physical_term": "",
            "worst_physical_term": "",
            "physical_terms_with_positive_ablation_effect": (),
            "physical_terms_rejected_as_noise": (),
        }
    best = min(rows, key=lambda row: row.delta_vs_full_selector)
    worst = max(rows, key=lambda row: row.delta_vs_full_selector)
    positive_terms = tuple(
        row.ablation_name.removeprefix("without_")
        for row in rows
        if row.delta_vs_full_selector < -0.005
    )
    rejected = tuple(
        row.ablation_name.removeprefix("without_")
        for row in rows
        if row.delta_vs_full_selector >= -0.005
    )
    return {
        "best_physical_term": best.ablation_name.removeprefix("without_"),
        "worst_physical_term": worst.ablation_name.removeprefix("without_"),
        "physical_terms_with_positive_ablation_effect": positive_terms,
        "physical_terms_rejected_as_noise": rejected,
    }


def ablation_interpretation(delta_vs_full_selector: float) -> str:
    if delta_vs_full_selector < -0.005:
        return "term_supports_selection"
    if delta_vs_full_selector > 0.005:
        return "term_adds_noise_in_current_proxy"
    return "term_effect_negligible"
