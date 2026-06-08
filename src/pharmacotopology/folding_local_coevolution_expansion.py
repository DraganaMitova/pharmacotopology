from __future__ import annotations

"""Local co-evolution expansion around sparse DCA anchors.

Hypothesis under test:
    A true long-range contact should often live in a local co-evolutionary
    patch around sparse DCA anchors.  Instead of drawing straight geodesic
    lines through contact-map space, expand each anchor in a small +/- window
    and keep only pairs that have local evolutionary support.

Boundary:
* Raw MSA is not bundled in the locked benchmark ZIP.  The default path is
  therefore an explicit proxy over safe external coupling rows: local APC/raw
  support, nearby coupling patch support, anchor-conditioned decay, and a small
  sequence-chemistry term.
* No native contacts or coordinates are used before the contact map is frozen.
* Native contacts are attached only after selection for audit.
"""

from dataclasses import asdict, dataclass
from functools import lru_cache
from math import exp
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.folding_coarse_grain_md_geometry import CoarseGrainMDContactDecision, _matched_controls_for_report, _pair_hash, _rounded
from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_native_contact_eval import ContactMetricPacket, ContactPair, evaluate_contact_prediction, normalized_contact_pairs
from pharmacotopology.folding_real_coordinate_visual_benchmark import RealCoordinateVisualRow
from pharmacotopology.folding_sequence_physical_priors import secondary_structure_pair_score, predict_lightweight_secondary_structure

LOCAL_COEVOLUTION_KIND = "local_coevolution_anchor_expansion_v0"
LOCAL_COEVOLUTION_MODE = "external_dca_anchor_local_coevolution_proxy"
LOCAL_COEVOLUTION_DECISION_KIND = "local_coevolution_contact_decision_v0"
LOCAL_COEVOLUTION_RULE = (
    "top_safe_dca_anchors;window_plus_minus_residues;"
    "raw_msa_required_for_true_mi_but_absent_in_zip;"
    "proxy_local_coupling_from_nearby_external_apc_raw_patch_support;"
    "anchor_conditioned_expansion;degree_cap;native_audit_after_selection"
)
LOCAL_COEVOLUTION_CLAIM_RULE = (
    "claim_requires_precision_recall_long_range_ge_0_70_and_control_margins;"
    "proxy_mi_disallows_universal_physical_law_claim"
)

LONG_RANGE_THRESHOLD = 24
LOCAL_MI_LONG_RANGE_FLOOR = 30


@dataclass(frozen=True)
class LocalCoevolutionContact:
    row_id: str
    source_accession: str
    i: int
    j: int
    score: float
    channel: str
    parent_anchor_i: int
    parent_anchor_j: int
    offset_i: int
    offset_j: int
    nearby_external_support: float
    anchor_conditioned_support: float
    patch_support: float
    chemistry_support: float
    true_msa_mi_used: bool = False
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class LocalCoevolutionRowReport:
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    raw_msa_available_for_true_local_mi: bool
    local_mi_channel_is_proxy_not_new_msa_calculation: bool
    safe_anchor_count: int
    local_window: int
    threshold: float
    candidate_local_pair_count: int
    accepted_local_pair_count: int
    selected_contact_count: int
    selected_long_range_contact_count: int
    selected_contact_map_hash: str
    ion_bridge_candidate_count: int
    metric_after_native_audit: ContactMetricPacket
    matched_control_count: int
    best_control_f1_after_audit: float
    best_control_long_range_recall_after_audit: float
    f1_margin_vs_best_control: float
    long_range_recall_margin_vs_best_control: float
    row_local_coevolution_claim_allowed: bool
    row_universal_physical_law_claim_allowed: bool
    row_claim_rejection_reason: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    target_native_truth_attached_after_selection_for_evaluation: bool = True
    learned_geometry_prior_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = True

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["metric_after_native_audit"] = self.metric_after_native_audit.to_dict()
        return payload


@dataclass(frozen=True)
class LocalCoevolutionExpansionPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    raw_msa_available_for_true_local_mi: bool
    local_mi_channel_is_proxy_not_new_msa_calculation: bool
    top_safe_dca_anchor_count_requested: int
    local_window: int
    threshold: float
    mean_safe_anchor_count: float
    mean_accepted_local_pair_count: float
    mean_native_contact_precision_after_audit: float
    mean_native_contact_recall_after_audit: float
    mean_long_range_contact_recall_after_audit: float
    mean_contact_map_f1_after_audit: float
    min_row_precision_after_audit: float
    min_row_recall_after_audit: float
    min_row_long_range_recall_after_audit: float
    mean_f1_margin_vs_best_control: float
    mean_long_range_recall_margin_vs_best_control: float
    local_coevolution_claim_allowed: bool
    universal_physical_law_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    rows: tuple[LocalCoevolutionRowReport, ...]
    contacts: tuple[LocalCoevolutionContact, ...]
    decisions: tuple[CoarseGrainMDContactDecision, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "source_mode": self.source_mode,
            "row_count": self.row_count,
            "decision_rule": self.decision_rule,
            "claim_rule": self.claim_rule,
            "raw_msa_available_for_true_local_mi": self.raw_msa_available_for_true_local_mi,
            "local_mi_channel_is_proxy_not_new_msa_calculation": self.local_mi_channel_is_proxy_not_new_msa_calculation,
            "top_safe_dca_anchor_count_requested": self.top_safe_dca_anchor_count_requested,
            "local_window": self.local_window,
            "threshold": self.threshold,
            "mean_safe_anchor_count": self.mean_safe_anchor_count,
            "mean_accepted_local_pair_count": self.mean_accepted_local_pair_count,
            "mean_native_contact_precision_after_audit": self.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": self.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": self.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": self.mean_contact_map_f1_after_audit,
            "min_row_precision_after_audit": self.min_row_precision_after_audit,
            "min_row_recall_after_audit": self.min_row_recall_after_audit,
            "min_row_long_range_recall_after_audit": self.min_row_long_range_recall_after_audit,
            "mean_f1_margin_vs_best_control": self.mean_f1_margin_vs_best_control,
            "mean_long_range_recall_margin_vs_best_control": self.mean_long_range_recall_margin_vs_best_control,
            "local_coevolution_claim_allowed": self.local_coevolution_claim_allowed,
            "universal_physical_law_claim_allowed": self.universal_physical_law_claim_allowed,
            "folding_problem_solved": self.folding_problem_solved,
            "claim_rejection_reason": self.claim_rejection_reason,
            "rows": [r.to_dict() for r in self.rows],
            "contacts": [c.to_dict() for c in self.contacts],
            "decisions": [d.to_dict() for d in self.decisions],
        }


def _bounded01(value: float) -> float:
    return round(max(0.0, min(1.0, float(value))), 6)


def _mean(values: Sequence[float]) -> float:
    return _rounded(mean(values)) if values else 0.0


def _safe_anchors(row: RealCoordinateVisualRow, constraints: Sequence[CouplingConstraint], *, top_n: int) -> tuple[CouplingConstraint, ...]:
    safe = [
        c
        for c in constraints
        if (c.row_id == row.row_id or c.source_accession == row.source_accession)
        and not c.coordinate_truth_used_to_build_constraint
        and not c.native_truth_used_before_coupling_selection
        and not c.structure_model_used
        and c.i < c.j
        and c.sequence_separation >= LONG_RANGE_THRESHOLD
    ]
    safe.sort(key=lambda c: (-(c.apc_corrected_score or c.confidence), c.rank if c.rank else 999999, c.i, c.j))
    return tuple(safe[: max(1, top_n)])


def _constraint_support_map(row: RealCoordinateVisualRow, constraints: Sequence[CouplingConstraint]) -> dict[ContactPair, float]:
    out: dict[ContactPair, float] = {}
    for c in constraints:
        if (c.row_id != row.row_id and c.source_accession != row.source_accession) or c.i >= c.j:
            continue
        if c.coordinate_truth_used_to_build_constraint or c.native_truth_used_before_coupling_selection or c.structure_model_used:
            continue
        score = max(float(c.confidence), float(c.apc_corrected_score or 0.0), min(1.0, float(c.raw_score or 0.0)))
        pair = (c.i, c.j)
        out[pair] = max(out.get(pair, 0.0), _bounded01(score))
    return out


def _nearby_external_support(pair: ContactPair, support: Mapping[ContactPair, float], *, radius: int = 2) -> float:
    i, j = pair
    best = 0.0
    for di in range(-radius, radius + 1):
        for dj in range(-radius, radius + 1):
            q = (i + di, j + dj)
            if q[0] >= q[1]:
                continue
            decay = exp(-(abs(di) + abs(dj)) / 2.5)
            best = max(best, float(support.get(q, 0.0)) * decay)
    return _bounded01(best)


def _patch_support(pair: ContactPair, support: Mapping[ContactPair, float], *, radius: int = 3) -> float:
    i, j = pair
    total = 0.0
    count = 0
    for di in range(-radius, radius + 1):
        for dj in range(-radius, radius + 1):
            if di == 0 and dj == 0:
                continue
            q = (i + di, j + dj)
            if q[0] >= q[1]:
                continue
            if q in support:
                total += float(support[q]) * exp(-(abs(di) + abs(dj)) / 5.0)
                count += 1
    return _bounded01(total / max(1.0, min(12.0, float(count) + 2.0)))


@lru_cache(maxsize=64)
def _cached_ss(sequence: str) -> tuple[str, ...]:
    return predict_lightweight_secondary_structure(sequence)

def _chemistry_support(row: RealCoordinateVisualRow, pair: ContactPair) -> float:
    i, j = pair
    ss = _cached_ss(row.sequence)
    ss_score = secondary_structure_pair_score(pair, ss)
    a = row.sequence[i - 1]
    b = row.sequence[j - 1]
    hydrophobic = set("AILMFWVY")
    charged_pos = set("KRH")
    charged_neg = set("DE")
    pair_score = 0.0
    if a in hydrophobic and b in hydrophobic:
        pair_score = 0.78
    elif (a in charged_pos and b in charged_neg) or (a in charged_neg and b in charged_pos):
        pair_score = 0.72
    elif a in "STNQ" and b in "STNQ":
        pair_score = 0.46
    elif a == "G" or b == "G":
        pair_score = 0.22
    else:
        pair_score = 0.34
    return _bounded01(0.55 * ss_score + 0.45 * pair_score)


def _anchor_conditioned_support(anchor: CouplingConstraint, pair: ContactPair) -> float:
    di = abs(pair[0] - anchor.i)
    dj = abs(pair[1] - anchor.j)
    # This is the user's local co-evolution intuition translated to a bounded
    # proxy: nearby residues inherit some anchor support, but the support decays
    # sharply and never becomes evidence by itself unless the threshold is low.
    return _bounded01(float(anchor.confidence) * exp(-(di + dj) / 6.0))


def build_local_coevolution_contacts(
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    *,
    top_anchor_count: int = 50,
    window: int = 5,
    threshold: float = 0.5,
    degree_cap: int = 5,
) -> tuple[tuple[LocalCoevolutionContact, ...], tuple[ContactPair, ...], dict[ContactPair, float]]:
    anchors = _safe_anchors(row, constraints, top_n=top_anchor_count)
    support = _constraint_support_map(row, constraints)
    accepted: dict[ContactPair, LocalCoevolutionContact] = {}
    scores: dict[ContactPair, float] = {}

    def add_contact(record: LocalCoevolutionContact) -> None:
        pair = record.pair()
        prev = accepted.get(pair)
        if prev is None or record.score > prev.score:
            accepted[pair] = record
            scores[pair] = record.score

    for anchor in anchors:
        direct_pair = anchor.pair()
        direct_score = max(float(anchor.confidence), float(anchor.apc_corrected_score or 0.0))
        direct_record = LocalCoevolutionContact(
            row_id=row.row_id,
            source_accession=row.source_accession,
            i=direct_pair[0],
            j=direct_pair[1],
            score=_bounded01(direct_score),
            channel="direct_safe_dca_anchor",
            parent_anchor_i=anchor.i,
            parent_anchor_j=anchor.j,
            offset_i=0,
            offset_j=0,
            nearby_external_support=_bounded01(direct_score),
            anchor_conditioned_support=_bounded01(direct_score),
            patch_support=_patch_support(direct_pair, support),
            chemistry_support=_chemistry_support(row, direct_pair),
        )
        add_contact(direct_record)
        for di in range(-window, window + 1):
            for dj in range(-window, window + 1):
                if di == 0 and dj == 0:
                    continue
                i2 = anchor.i + di
                j2 = anchor.j + dj
                if i2 < 1 or j2 > row.sequence_length or i2 >= j2 or j2 - i2 < LOCAL_MI_LONG_RANGE_FLOOR:
                    continue
                pair = (i2, j2)
                nearby = _nearby_external_support(pair, support)
                conditioned = _anchor_conditioned_support(anchor, pair)
                patch = _patch_support(pair, support)
                chemistry = _chemistry_support(row, pair)
                # Proxy for compute_local_coupling(msa_filtered, i2, j2): the
                # candidate must have local external support or coherent patch
                # support; anchor inheritance alone is down-weighted.
                local_score = _bounded01(0.44 * nearby + 0.27 * patch + 0.20 * conditioned + 0.09 * chemistry)
                if local_score < threshold:
                    continue
                add_contact(LocalCoevolutionContact(
                    row_id=row.row_id,
                    source_accession=row.source_accession,
                    i=i2,
                    j=j2,
                    score=local_score,
                    channel="local_coevolution_proxy_window",
                    parent_anchor_i=anchor.i,
                    parent_anchor_j=anchor.j,
                    offset_i=di,
                    offset_j=dj,
                    nearby_external_support=nearby,
                    anchor_conditioned_support=conditioned,
                    patch_support=patch,
                    chemistry_support=chemistry,
                ))

    ranked = sorted(accepted.values(), key=lambda r: (-r.score, r.channel != "direct_safe_dca_anchor", r.i, r.j))
    degrees = [0] * row.sequence_length
    selected: list[ContactPair] = []
    selected_scores: dict[ContactPair, float] = {}
    for record in ranked:
        i, j = record.pair()
        if degrees[i - 1] >= degree_cap or degrees[j - 1] >= degree_cap:
            continue
        selected.append((i, j))
        selected_scores[(i, j)] = record.score
        degrees[i - 1] += 1
        degrees[j - 1] += 1
    return tuple(ranked), normalized_contact_pairs(selected), selected_scores


def _ion_bridge_candidate_count(row: RealCoordinateVisualRow, selected: Sequence[ContactPair]) -> int:
    pos = set("KRH")
    neg = set("DE")
    count = 0
    for i, j in selected:
        a = row.sequence[i - 1]
        b = row.sequence[j - 1]
        if j - i >= LOCAL_MI_LONG_RANGE_FLOOR and ((a in pos and b in neg) or (a in neg and b in pos)):
            count += 1
    return count


def _claim(metric: ContactMetricPacket, f1_margin: float, lr_margin: float) -> tuple[bool, bool, str]:
    claim = (
        metric.native_contact_precision >= 0.70
        and metric.native_contact_recall >= 0.70
        and metric.long_range_contact_recall >= 0.70
        and f1_margin > 0.03
        and lr_margin > 0.03
    )
    reason = "local_coevolution_expansion_survived_gate_not_universal_due_to_proxy_mi" if claim else "local_coevolution_claim_rejected_precision_recall_or_long_range_below_0_70"
    return claim, False, reason


def _decisions(row: RealCoordinateVisualRow, selected: Sequence[ContactPair], scores: Mapping[ContactPair, float]) -> tuple[CoarseGrainMDContactDecision, ...]:
    out: list[CoarseGrainMDContactDecision] = []
    for i, j in normalized_contact_pairs(selected):
        out.append(CoarseGrainMDContactDecision(
            kind=LOCAL_COEVOLUTION_DECISION_KIND,
            row_id=row.row_id,
            source_accession=row.source_accession,
            source_mode=LOCAL_COEVOLUTION_MODE,
            i=i,
            j=j,
            sequence_separation=j - i,
            final_distance_angstrom=0.0,
            geometry_contact_score=_bounded01(scores.get((i, j), 0.0)),
            selected=True,
            selected_from_final_structure=False,
            msa_dca_used_before_selection=True,
        ))
    return tuple(out)


def run_local_coevolution_expansion_row(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    top_anchor_count: int = 50,
    window: int = 5,
    threshold: float = 0.5,
    degree_cap: int = 5,
) -> tuple[LocalCoevolutionRowReport, tuple[LocalCoevolutionContact, ...], tuple[CoarseGrainMDContactDecision, ...]]:
    contacts, selected, scores = build_local_coevolution_contacts(
        row,
        constraints,
        top_anchor_count=top_anchor_count,
        window=window,
        threshold=threshold,
        degree_cap=degree_cap,
    )
    metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=selected)
    best_control_f1, best_control_lr, control_count = _matched_controls_for_report(row=row, selected_pairs=selected)
    f1_margin = _rounded(metric.contact_map_f1 - best_control_f1)
    lr_margin = _rounded(metric.long_range_contact_recall - best_control_lr)
    claim, universal, reason = _claim(metric, f1_margin, lr_margin)
    safe_anchor_count = len(_safe_anchors(row, constraints, top_n=top_anchor_count))
    accepted_local = sum(1 for c in contacts if c.channel == "local_coevolution_proxy_window")
    candidate_local = safe_anchor_count * ((2 * window + 1) ** 2 - 1)
    report = LocalCoevolutionRowReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        raw_msa_available_for_true_local_mi=False,
        local_mi_channel_is_proxy_not_new_msa_calculation=True,
        safe_anchor_count=safe_anchor_count,
        local_window=window,
        threshold=_rounded(threshold),
        candidate_local_pair_count=candidate_local,
        accepted_local_pair_count=accepted_local,
        selected_contact_count=len(selected),
        selected_long_range_contact_count=sum(1 for i, j in selected if j - i >= LONG_RANGE_THRESHOLD),
        selected_contact_map_hash=_pair_hash(selected),
        ion_bridge_candidate_count=_ion_bridge_candidate_count(row, selected),
        metric_after_native_audit=metric,
        matched_control_count=control_count,
        best_control_f1_after_audit=best_control_f1,
        best_control_long_range_recall_after_audit=best_control_lr,
        f1_margin_vs_best_control=f1_margin,
        long_range_recall_margin_vs_best_control=lr_margin,
        row_local_coevolution_claim_allowed=claim,
        row_universal_physical_law_claim_allowed=universal,
        row_claim_rejection_reason=reason,
    )
    return report, contacts, _decisions(row, selected, scores)


def run_local_coevolution_expansion_packet(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    constraints: Sequence[CouplingConstraint],
    evaluation_source_accessions: Sequence[str] | None = None,
    top_anchor_count: int = 50,
    window: int = 5,
    threshold: float = 0.5,
    degree_cap: int = 5,
) -> LocalCoevolutionExpansionPacket:
    evaluation_set = set(evaluation_source_accessions or ())
    eval_rows = tuple(row for row in rows if not evaluation_set or row.source_accession in evaluation_set)
    reports: list[LocalCoevolutionRowReport] = []
    all_contacts: list[LocalCoevolutionContact] = []
    all_decisions: list[CoarseGrainMDContactDecision] = []
    for row in eval_rows:
        report, contacts, decisions = run_local_coevolution_expansion_row(
            row=row,
            constraints=constraints,
            top_anchor_count=top_anchor_count,
            window=window,
            threshold=threshold,
            degree_cap=degree_cap,
        )
        reports.append(report)
        all_contacts.extend(contacts[:800])
        all_decisions.extend(decisions[:500])
    precisions = [r.metric_after_native_audit.native_contact_precision for r in reports]
    recalls = [r.metric_after_native_audit.native_contact_recall for r in reports]
    long_ranges = [r.metric_after_native_audit.long_range_contact_recall for r in reports]
    f1s = [r.metric_after_native_audit.contact_map_f1 for r in reports]
    f1_margins = [r.f1_margin_vs_best_control for r in reports]
    lr_margins = [r.long_range_recall_margin_vs_best_control for r in reports]
    claim = bool(reports) and all(r.row_local_coevolution_claim_allowed for r in reports) and _mean(precisions) >= 0.70 and _mean(recalls) >= 0.70 and _mean(long_ranges) >= 0.70
    reason = "local_coevolution_expansion_survived_gate_not_universal_physics" if claim else "local_coevolution_claim_rejected_for_rows:" + ",".join(r.source_accession for r in reports if not r.row_local_coevolution_claim_allowed)[:180]
    return LocalCoevolutionExpansionPacket(
        kind=LOCAL_COEVOLUTION_KIND,
        source_mode=LOCAL_COEVOLUTION_MODE,
        row_count=len(reports),
        decision_rule=LOCAL_COEVOLUTION_RULE,
        claim_rule=LOCAL_COEVOLUTION_CLAIM_RULE,
        raw_msa_available_for_true_local_mi=False,
        local_mi_channel_is_proxy_not_new_msa_calculation=True,
        top_safe_dca_anchor_count_requested=top_anchor_count,
        local_window=window,
        threshold=_rounded(threshold),
        mean_safe_anchor_count=_mean([r.safe_anchor_count for r in reports]),
        mean_accepted_local_pair_count=_mean([r.accepted_local_pair_count for r in reports]),
        mean_native_contact_precision_after_audit=_mean(precisions),
        mean_native_contact_recall_after_audit=_mean(recalls),
        mean_long_range_contact_recall_after_audit=_mean(long_ranges),
        mean_contact_map_f1_after_audit=_mean(f1s),
        min_row_precision_after_audit=_rounded(min(precisions)) if precisions else 0.0,
        min_row_recall_after_audit=_rounded(min(recalls)) if recalls else 0.0,
        min_row_long_range_recall_after_audit=_rounded(min(long_ranges)) if long_ranges else 0.0,
        mean_f1_margin_vs_best_control=_mean(f1_margins),
        mean_long_range_recall_margin_vs_best_control=_mean(lr_margins),
        local_coevolution_claim_allowed=claim,
        universal_physical_law_claim_allowed=False,
        folding_problem_solved=claim,
        claim_rejection_reason=reason,
        rows=tuple(reports),
        contacts=tuple(all_contacts),
        decisions=tuple(all_decisions),
    )
