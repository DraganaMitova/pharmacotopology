from __future__ import annotations

from statistics import pstdev
from typing import Mapping, Sequence


CONTACT_LAW_FALSIFICATION_KIND = "contact_law_threshold_falsification_v1"


def threshold_std(thresholds: Sequence[float]) -> float:
    if len(thresholds) <= 1:
        return 0.0
    return round(pstdev(thresholds), 6)


def threshold_is_stable(thresholds: Sequence[float], *, max_std: float = 0.10) -> bool:
    return threshold_std(thresholds) < max_std


def current_scalar_law_rejected(
    *,
    best_global_f1: float,
    threshold_stable: bool,
    minimum_survival_f1: float = 0.22,
) -> bool:
    return (not threshold_stable) or best_global_f1 < minimum_survival_f1


def candidate_law_survives(
    *,
    candidate: Mapping[str, object],
    pair_only: Mapping[str, object],
    min_f1_gain: float = 0.02,
    max_threshold_std: float = 0.10,
) -> bool:
    candidate_f1 = float(candidate["loo_mean_test_f1"])
    pair_only_f1 = float(pair_only["loo_mean_test_f1"])
    candidate_false_rate = float(candidate["loo_mean_false_contact_rate"])
    pair_only_false_rate = float(pair_only["loo_mean_false_contact_rate"])
    candidate_long_recall = float(candidate["loo_mean_long_range_contact_recall"])
    pair_only_long_recall = float(pair_only["loo_mean_long_range_contact_recall"])
    candidate_std = float(candidate["loo_threshold_std"])
    return (
        candidate_f1 > pair_only_f1 + min_f1_gain
        and candidate_false_rate < pair_only_false_rate
        and candidate_std < max_threshold_std
        and candidate_long_recall > pair_only_long_recall
    )


def falsification_summary(
    *,
    current_scalar: Mapping[str, object],
    best_candidate: Mapping[str, object],
    pair_only: Mapping[str, object],
) -> dict[str, object]:
    current_rejected = current_scalar_law_rejected(
        best_global_f1=float(current_scalar["best_global_f1"]),
        threshold_stable=bool(current_scalar["threshold_stable"]),
    )
    candidate_survives = candidate_law_survives(
        candidate=best_candidate,
        pair_only=pair_only,
    )
    return {
        "falsification_kind": CONTACT_LAW_FALSIFICATION_KIND,
        "current_scalar_score_law_rejected": current_rejected,
        "candidate_contact_law_survives": candidate_survives,
        "law_generalizes": candidate_survives,
        "mechanism_discovery_claim_allowed": False,
        "folding_problem_solved": False,
    }
