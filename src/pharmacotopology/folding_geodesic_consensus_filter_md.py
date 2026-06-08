from __future__ import annotations

"""Consensus-filtered dense-anchor geodesic restrained global-geometry challenge.

This layer tests the next falsifiable hypothesis after dense-anchor geodesic MD:

    dense anchors are not the missing ingredient by themselves;
    the missing operation is selecting true geodesic fills and rejecting false fills.

The layer therefore keeps the structure-first path, but replaces naive geodesic
interpolation with native-free filters:

* local contact-map geodesic distance between anchor endpoints;
* a sequence-only bending-cost proxy;
* multi-path consensus voting before a geodesic fill can become a restraint.

No native contacts, native coordinates, templates, AlphaFold/ESMFold, or learned
geometry prior are used to create the restraints. Native is only used after the
final C-alpha geometry has been converted to contacts for audit.
"""

import hashlib
from dataclasses import asdict, dataclass
from math import exp
from statistics import mean
from typing import Sequence

from pharmacotopology.folding_coarse_grain_md_geometry import (
    COARSE_GRAIN_MD_CLAIM_RULE,
    COARSE_GRAIN_MD_RULE,
    CoarseGrainMDGeometryPacket,
    run_coarse_grain_md_geometry_packet,
)
from pharmacotopology.folding_dense_anchor_geodesic_md import (
    DENSE_ANCHOR_GEODESIC_MODE,
    DenseAnchorRecord,
    _constraint,
    _local_dca_support,
    _records_from_constraints,
    _safe_direct_constraints,
    _sequence_closure_candidates,
)
from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_native_contact_eval import ContactPair, normalized_contact_pairs
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    MIN_SEQUENCE_SEPARATION,
    RealCoordinateVisualRow,
)

GEODESIC_CONSENSUS_FILTER_KIND = "geodesic_consensus_filter_restrained_md_geometry_v0"
GEODESIC_CONSENSUS_FILTER_MODE = "external_geodesic_consensus_filter_md_geometry"
GEODESIC_CONSENSUS_FILTER_RULE = (
    COARSE_GRAIN_MD_RULE
    + ";dense_anchor_expansion;anchor_chaining;geodesic_distance_filter;"
    + "bending_cost_filter;multipath_consensus_geodesic_restraints"
)
GEODESIC_CONSENSUS_FILTER_CLAIM_RULE = (
    COARSE_GRAIN_MD_CLAIM_RULE
    + ";consensus_filter_claim_requires_precision_recall_long_range_ge_0_70;"
    + "native_used_only_after_structure_to_contact_audit"
)


@dataclass(frozen=True)
class ConsensusFilterRowSummary:
    row_id: str
    source_accession: str
    sequence_length: int
    direct_dca_anchor_count: int
    sequence_closure_anchor_count: int
    beta_bridge_anchor_count: int
    candidate_anchor_count_before_geodesic: int
    candidate_geodesic_contact_count: int
    rejected_by_distance_filter_count: int
    rejected_by_bending_filter_count: int
    rejected_by_consensus_filter_count: int
    accepted_geodesic_contact_count: int
    total_restraint_count: int
    restraint_map_hash: str
    max_allowed_geodesic_distance: int
    max_allowed_bending_energy: float
    consensus_vote_threshold: int
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class GeodesicConsensusFilterPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    direct_global_structure_generation_included: bool
    contacts_predicted_before_structure: bool
    dense_anchor_expansion_included: bool
    geodesic_distance_filter_included: bool
    bending_energy_filter_included: bool
    consensus_geodesic_filter_included: bool
    atomistic_md_engine_used: bool
    dependency_free_md_style_relaxation_used: bool
    mean_total_restraint_count: float
    mean_accepted_geodesic_contact_count: float
    mean_rejected_geodesic_contact_count: float
    md_packet: CoarseGrainMDGeometryPacket
    filter_rows: tuple[ConsensusFilterRowSummary, ...]
    restraints: tuple[DenseAnchorRecord, ...]
    geodesic_consensus_filter_claim_allowed: bool
    universal_physical_law_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "source_mode": self.source_mode,
            "row_count": self.row_count,
            "decision_rule": self.decision_rule,
            "claim_rule": self.claim_rule,
            "direct_global_structure_generation_included": self.direct_global_structure_generation_included,
            "contacts_predicted_before_structure": self.contacts_predicted_before_structure,
            "dense_anchor_expansion_included": self.dense_anchor_expansion_included,
            "geodesic_distance_filter_included": self.geodesic_distance_filter_included,
            "bending_energy_filter_included": self.bending_energy_filter_included,
            "consensus_geodesic_filter_included": self.consensus_geodesic_filter_included,
            "atomistic_md_engine_used": self.atomistic_md_engine_used,
            "dependency_free_md_style_relaxation_used": self.dependency_free_md_style_relaxation_used,
            "mean_total_restraint_count": self.mean_total_restraint_count,
            "mean_accepted_geodesic_contact_count": self.mean_accepted_geodesic_contact_count,
            "mean_rejected_geodesic_contact_count": self.mean_rejected_geodesic_contact_count,
            "md_packet": self.md_packet.to_dict(),
            "filter_rows": [row.to_dict() for row in self.filter_rows],
            "restraints": [restraint.to_dict() for restraint in self.restraints],
            "geodesic_consensus_filter_claim_allowed": self.geodesic_consensus_filter_claim_allowed,
            "universal_physical_law_claim_allowed": self.universal_physical_law_claim_allowed,
            "folding_problem_solved": self.folding_problem_solved,
            "claim_rejection_reason": self.claim_rejection_reason,
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "structure_model_used_before_selection": self.structure_model_used_before_selection,
            "learned_geometry_prior_used_before_selection": self.learned_geometry_prior_used_before_selection,
        }


def _rounded(value: float) -> float:
    return round(float(value), 6)


def _bounded01(value: float) -> float:
    return round(max(0.0, min(1.0, float(value))), 6)


def _pair_hash(pairs: Sequence[ContactPair]) -> str:
    encoded = ";".join(f"{i}-{j}" for i, j in normalized_contact_pairs(pairs)).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def estimate_geodesic_distance(anchor1: ContactPair, anchor2: ContactPair) -> int:
    """Native-free contact-map distance between two anchor endpoints."""
    return max(abs(anchor2[0] - anchor1[0]), abs(anchor2[1] - anchor1[1]))


def estimate_bending_energy(anchor1: ContactPair, anchor2: ContactPair) -> float:
    """Sequence-only proxy for sharp-bend cost between neighboring anchors."""
    gap_i = abs(anchor2[0] - anchor1[0])
    gap_j = abs(anchor2[1] - anchor1[1])
    return _rounded(1.0 / (min(gap_i, gap_j) + 1.0))


def filter_geodesic_anchor_pair(
    anchor1: ContactPair,
    anchor2: ContactPair,
    *,
    max_geodesic_distance: int,
    max_bending_energy: float,
) -> tuple[bool, str]:
    distance = estimate_geodesic_distance(anchor1, anchor2)
    if distance > max_geodesic_distance:
        return False, "distance"
    bending = estimate_bending_energy(anchor1, anchor2)
    if bending > max_bending_energy:
        return False, "bending"
    return True, "accepted"


def _anchor_chains(pairs: Sequence[ContactPair], *, max_gap: int = 11) -> list[list[ContactPair]]:
    remaining = list(normalized_contact_pairs(pairs))
    chains: list[list[ContactPair]] = []
    while remaining:
        chain = [remaining.pop(0)]
        changed = True
        while changed:
            changed = False
            for pair in remaining[:]:
                if any(abs(pair[0] - c[0]) <= max_gap and abs(pair[1] - c[1]) <= max_gap for c in chain):
                    chain.append(pair)
                    remaining.remove(pair)
                    changed = True
        chain.sort()
        chains.append(chain)
    chains.sort(key=lambda ch: (-len(ch), ch[0]))
    return chains


def _linear_path(a: ContactPair, b: ContactPair, *, offset_i: int = 0, offset_j: int = 0) -> list[ContactPair]:
    i1, j1 = a
    i2, j2 = b
    steps = estimate_geodesic_distance(a, b)
    if steps <= 1 or steps > 64:
        return []
    out: list[ContactPair] = []
    for t in range(1, steps):
        frac = t / steps
        i = int(round(i1 + (i2 - i1) * frac)) + offset_i
        j = int(round(j1 + (j2 - j1) * frac)) + offset_j
        if i < 1 or j < 1:
            continue
        left, right = min(i, j), max(i, j)
        if right - left >= MIN_SEQUENCE_SEPARATION:
            out.append((left, right))
    return normalized_contact_pairs(out)


def _staircase_path(a: ContactPair, b: ContactPair, *, i_first: bool) -> list[ContactPair]:
    i, j = a
    i2, j2 = b
    out: list[ContactPair] = []

    def step_axis(value: int, target: int) -> int:
        if value == target:
            return value
        return value + (1 if target > value else -1)

    while (i, j) != (i2, j2) and len(out) < 128:
        if i_first:
            i = step_axis(i, i2)
            if i == i2:
                j = step_axis(j, j2)
        else:
            j = step_axis(j, j2)
            if j == j2:
                i = step_axis(i, i2)
        left, right = min(i, j), max(i, j)
        if right - left >= MIN_SEQUENCE_SEPARATION and (left, right) not in (a, b):
            out.append((left, right))
    return normalized_contact_pairs(out)


def consensus_geodesic_contacts(
    anchor1: ContactPair,
    anchor2: ContactPair,
    *,
    sequence_length: int,
    vote_threshold: int,
) -> list[ContactPair]:
    """Return contacts that survive a small deterministic ensemble of paths."""
    paths = [
        _linear_path(anchor1, anchor2),
        _linear_path(anchor1, anchor2),
        _linear_path(anchor1, anchor2),
        _linear_path(anchor1, anchor2, offset_i=1),
        _linear_path(anchor1, anchor2, offset_i=-1),
        _linear_path(anchor1, anchor2, offset_j=1),
        _linear_path(anchor1, anchor2, offset_j=-1),
        _staircase_path(anchor1, anchor2, i_first=True),
        _staircase_path(anchor1, anchor2, i_first=False),
    ]
    votes: dict[ContactPair, int] = {}
    for path in paths:
        for pair in set(path):
            if pair[0] < 1 or pair[1] > sequence_length:
                continue
            votes[pair] = votes.get(pair, 0) + 1
    return normalized_contact_pairs([pair for pair, count in votes.items() if count >= vote_threshold])


def _base_anchor_candidates(
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    *,
    max_direct: int,
    max_sequence_closure: int,
) -> tuple[dict[ContactPair, CouplingConstraint], list[DenseAnchorRecord], dict[str, int], Sequence[CouplingConstraint]]:
    direct = list(_safe_direct_constraints(row, constraints, max_direct=max_direct))
    by_pair: dict[ContactPair, CouplingConstraint] = {c.pair(): c for c in direct}
    records: list[DenseAnchorRecord] = _records_from_constraints(row, direct, "direct_dca")
    counts = {"direct_dca": len(direct), "mi_surrogate_sequence_closure": 0, "beta_bridge": 0}
    rank = len(direct) + 1
    for pair, score, channel in _sequence_closure_candidates(row, direct, max_candidates=max_sequence_closure):
        if pair in by_pair:
            continue
        confidence = _bounded01(score * (0.74 if channel == "mi_surrogate_sequence_closure" else 0.62))
        by_pair[pair] = _constraint(row=row, pair=pair, confidence=confidence, channel=channel, rank=rank)
        records.append(DenseAnchorRecord(row.row_id, row.source_accession, pair[0], pair[1], confidence, channel))
        counts[channel] = counts.get(channel, 0) + 1
        rank += 1
    return by_pair, records, counts, tuple(direct)


def build_geodesic_consensus_filter_constraints(
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    *,
    max_direct: int | None = None,
    max_sequence_closure: int | None = None,
    max_geodesic: int | None = None,
    max_geodesic_distance: int = 30,
    max_bending_energy: float = 0.50,
    consensus_vote_threshold: int = 3,
) -> tuple[tuple[CouplingConstraint, ...], ConsensusFilterRowSummary, tuple[DenseAnchorRecord, ...]]:
    n = row.sequence_length
    direct_limit = max_direct or max(48, min(n, 256))
    closure_limit = max_sequence_closure or max(72, min(n, 280))
    geodesic_limit = max_geodesic or max(48, min(n, 260))
    by_pair, records, counts, direct = _base_anchor_candidates(
        row,
        constraints,
        max_direct=direct_limit,
        max_sequence_closure=closure_limit,
    )
    seed_count = len(by_pair)
    chains = _anchor_chains(list(by_pair), max_gap=11)
    candidate_geodesic = 0
    rejected_distance = 0
    rejected_bending = 0
    rejected_consensus = 0
    accepted_geodesic = 0
    rank = len(by_pair) + 1

    for chain in chains:
        if len(chain) < 2:
            continue
        # adjacent and one-skip neighbors: enough to test local path consensus without all-pairs explosion
        neighbor_pairs: list[tuple[ContactPair, ContactPair]] = list(zip(chain[:-1], chain[1:]))
        neighbor_pairs.extend(zip(chain[:-2], chain[2:]))
        for left, right in neighbor_pairs:
            ok, reason = filter_geodesic_anchor_pair(
                left,
                right,
                max_geodesic_distance=max_geodesic_distance,
                max_bending_energy=max_bending_energy,
            )
            if not ok:
                if reason == "distance":
                    rejected_distance += 1
                elif reason == "bending":
                    rejected_bending += 1
                continue
            consensus = consensus_geodesic_contacts(
                left,
                right,
                sequence_length=n,
                vote_threshold=consensus_vote_threshold,
            )
            if not consensus:
                rejected_consensus += 1
                continue
            for pair in consensus:
                candidate_geodesic += 1
                if pair in by_pair:
                    rejected_consensus += 1
                    continue
                support = _local_dca_support(pair, direct)
                if support < 0.08:
                    rejected_consensus += 1
                    continue
                confidence = _bounded01(0.24 + 0.42 * support + 0.04 * exp(-estimate_geodesic_distance(left, right) / 20.0))
                by_pair[pair] = _constraint(row=row, pair=pair, confidence=confidence, channel="consensus_geodesic", rank=rank)
                records.append(
                    DenseAnchorRecord(
                        row.row_id,
                        row.source_accession,
                        pair[0],
                        pair[1],
                        confidence,
                        "consensus_geodesic",
                        parent_a=f"{left[0]}-{left[1]}",
                        parent_b=f"{right[0]}-{right[1]}",
                    )
                )
                rank += 1
                accepted_geodesic += 1
                if accepted_geodesic >= geodesic_limit:
                    break
            if accepted_geodesic >= geodesic_limit:
                break
        if accepted_geodesic >= geodesic_limit:
            break

    restraints = tuple(sorted(by_pair.values(), key=lambda c: (-c.confidence, c.i, c.j)))
    summary = ConsensusFilterRowSummary(
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_length=n,
        direct_dca_anchor_count=counts.get("direct_dca", 0),
        sequence_closure_anchor_count=counts.get("mi_surrogate_sequence_closure", 0),
        beta_bridge_anchor_count=counts.get("beta_bridge", 0),
        candidate_anchor_count_before_geodesic=seed_count,
        candidate_geodesic_contact_count=candidate_geodesic,
        rejected_by_distance_filter_count=rejected_distance,
        rejected_by_bending_filter_count=rejected_bending,
        rejected_by_consensus_filter_count=rejected_consensus,
        accepted_geodesic_contact_count=accepted_geodesic,
        total_restraint_count=len(restraints),
        restraint_map_hash=_pair_hash([c.pair() for c in restraints]),
        max_allowed_geodesic_distance=max_geodesic_distance,
        max_allowed_bending_energy=_rounded(max_bending_energy),
        consensus_vote_threshold=consensus_vote_threshold,
    )
    return restraints, summary, tuple(records)


def run_geodesic_consensus_filter_packet(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    constraints: Sequence[CouplingConstraint],
    md_steps: int = 72,
    restarts: int = 2,
    max_direct: int | None = None,
    max_sequence_closure: int | None = None,
    max_geodesic: int | None = None,
    max_geodesic_distance: int = 30,
    max_bending_energy: float = 0.50,
    consensus_vote_threshold: int = 3,
) -> GeodesicConsensusFilterPacket:
    all_restraints: list[CouplingConstraint] = []
    rows_out: list[ConsensusFilterRowSummary] = []
    records: list[DenseAnchorRecord] = []
    for row in rows:
        restraints, summary, row_records = build_geodesic_consensus_filter_constraints(
            row,
            constraints,
            max_direct=max_direct,
            max_sequence_closure=max_sequence_closure,
            max_geodesic=max_geodesic,
            max_geodesic_distance=max_geodesic_distance,
            max_bending_energy=max_bending_energy,
            consensus_vote_threshold=consensus_vote_threshold,
        )
        all_restraints.extend(restraints)
        rows_out.append(summary)
        records.extend(row_records[:900])

    md_packet = run_coarse_grain_md_geometry_packet(
        rows,
        constraints=tuple(all_restraints),
        source_mode=GEODESIC_CONSENSUS_FILTER_MODE,
        md_steps=md_steps,
        restarts=restarts,
        max_restraints=9999,
    )
    claim = (
        md_packet.md_geometry_claim_allowed
        and md_packet.mean_native_contact_precision_after_audit >= 0.70
        and md_packet.mean_native_contact_recall_after_audit >= 0.70
        and md_packet.mean_long_range_contact_recall_after_audit >= 0.70
    )
    if claim:
        rejection = "geodesic_consensus_filter_survived_external_geometry_gate_not_universal_physics"
    elif md_packet.mean_native_contact_precision_after_audit < 0.70:
        rejection = "geodesic_consensus_filter_rejected_precision_below_0_70"
    elif md_packet.mean_native_contact_recall_after_audit < 0.70:
        rejection = "geodesic_consensus_filter_rejected_recall_below_0_70"
    elif md_packet.mean_long_range_contact_recall_after_audit < 0.70:
        rejection = "geodesic_consensus_filter_rejected_long_range_recall_below_0_70"
    else:
        rejection = "geodesic_consensus_filter_rejected_control_or_row_gate"
    total_counts = [row.total_restraint_count for row in rows_out]
    accepted_counts = [row.accepted_geodesic_contact_count for row in rows_out]
    rejected_counts = [
        row.rejected_by_distance_filter_count + row.rejected_by_bending_filter_count + row.rejected_by_consensus_filter_count
        for row in rows_out
    ]
    return GeodesicConsensusFilterPacket(
        kind=GEODESIC_CONSENSUS_FILTER_KIND,
        source_mode=GEODESIC_CONSENSUS_FILTER_MODE,
        row_count=len(rows),
        decision_rule=GEODESIC_CONSENSUS_FILTER_RULE,
        claim_rule=GEODESIC_CONSENSUS_FILTER_CLAIM_RULE,
        direct_global_structure_generation_included=True,
        contacts_predicted_before_structure=False,
        dense_anchor_expansion_included=True,
        geodesic_distance_filter_included=True,
        bending_energy_filter_included=True,
        consensus_geodesic_filter_included=True,
        atomistic_md_engine_used=False,
        dependency_free_md_style_relaxation_used=True,
        mean_total_restraint_count=_rounded(mean(total_counts)) if total_counts else 0.0,
        mean_accepted_geodesic_contact_count=_rounded(mean(accepted_counts)) if accepted_counts else 0.0,
        mean_rejected_geodesic_contact_count=_rounded(mean(rejected_counts)) if rejected_counts else 0.0,
        md_packet=md_packet,
        filter_rows=tuple(rows_out),
        restraints=tuple(records),
        geodesic_consensus_filter_claim_allowed=claim,
        universal_physical_law_claim_allowed=False,
        folding_problem_solved=claim,
        claim_rejection_reason=rejection,
        coordinate_truth_used_before_selection=md_packet.coordinate_truth_used_before_selection,
        native_truth_used_before_selection=md_packet.native_truth_used_before_selection,
        structure_model_used_before_selection=md_packet.structure_model_used_before_selection,
        learned_geometry_prior_used_before_selection=md_packet.learned_geometry_prior_used_before_selection,
    )
