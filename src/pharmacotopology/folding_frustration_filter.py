from __future__ import annotations

from dataclasses import asdict, dataclass

from pharmacotopology.folding_closure_geometry import rounded_unit
from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent


FRUSTRATION_FILTER_KIND = "sequence_only_closure_frustration_filter_v1"

FRUSTRATION_COST_LIMIT = 0.08
GEOMETRY_COST_LIMIT = 0.32
ENTROPY_TRAP_LIMIT = 0.86
LOW_REGISTRY_LIMIT = 0.10
FALSE_CONTACT_RISK_LIMIT = 0.78


@dataclass(frozen=True)
class FrustrationAssessment:
    row_id: str
    event_id: str
    false_contact_risk_proxy: float
    cluster_precision_proxy: float
    hydrophobic_exposure_mismatch: float
    long_span_early_pressure: float
    passed_filter: bool
    primary_rejection_reason: str
    native_truth_used_before_filtering: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def false_contact_risk_proxy(event: NucleusClosureEvent) -> float:
    return rounded_unit(
        0.24 * (1.0 - event.contact_cluster_gain)
        + 0.22 * (1.0 - event.registry_support)
        + 0.18 * event.loop_entropy_cost
        + 0.12 * event.isolation_penalty
        + 0.10 * event.geometry_violation_cost
        + 0.08 * event.frustration_cost
        + (0.06 if event.normalized_span >= 0.62 else 0.0)
    )


def cluster_precision_proxy(event: NucleusClosureEvent) -> float:
    return rounded_unit(
        0.44 * event.contact_cluster_gain
        + 0.28 * event.registry_support
        + 0.18 * event.closure_event_stability
        + 0.10 * event.hydrophobic_burial_gain
        - 0.18 * event.loop_entropy_cost
        - 0.10 * event.isolation_penalty
    )


def hydrophobic_exposure_mismatch(event: NucleusClosureEvent) -> float:
    return rounded_unit(
        max(0.0, 0.72 - event.hydrophobic_burial_gain)
        + 0.35 * event.loop_entropy_cost
        - 0.20 * event.contact_cluster_gain
    )


def long_span_early_pressure(event: NucleusClosureEvent) -> float:
    return rounded_unit(
        max(0.0, event.normalized_span - 0.58)
        + 0.30 * event.geometry_violation_cost
        + 0.18 * event.loop_entropy_cost
    )


def assess_event_frustration(event: NucleusClosureEvent) -> FrustrationAssessment:
    risk = false_contact_risk_proxy(event)
    precision_proxy = cluster_precision_proxy(event)
    mismatch = hydrophobic_exposure_mismatch(event)
    span_pressure = long_span_early_pressure(event)
    reason = "passed"
    if event.frustration_cost >= FRUSTRATION_COST_LIMIT:
        reason = "frustration_rejection"
    elif event.geometry_violation_cost >= GEOMETRY_COST_LIMIT:
        reason = "geometry_rejection"
    elif (
        event.loop_entropy_cost >= ENTROPY_TRAP_LIMIT
        and event.registry_support <= LOW_REGISTRY_LIMIT
    ):
        reason = "trap_rejection"
    elif risk >= FALSE_CONTACT_RISK_LIMIT and precision_proxy < 0.22:
        reason = "trap_rejection"
    return FrustrationAssessment(
        row_id=event.row_id,
        event_id=event.event_id,
        false_contact_risk_proxy=risk,
        cluster_precision_proxy=precision_proxy,
        hydrophobic_exposure_mismatch=mismatch,
        long_span_early_pressure=span_pressure,
        passed_filter=reason == "passed",
        primary_rejection_reason=reason,
    )

