from __future__ import annotations

"""Dense-anchor geodesic restrained global-geometry challenge.

This layer tests the user's proposed next move:

    sparse external DCA anchors
    -> dense native-free anchor expansion
    -> anchor chaining / geodesic interpolation in contact-map space
    -> structure-first coarse C-alpha global relaxation
    -> contacts extracted from the final geometry
    -> native audit only after selection.

It deliberately does not use native contacts, native coordinates, templates,
ESMFold/AlphaFold, or any learned geometry prior to build anchors.  Raw MSA is
not bundled in the benchmark, so the "MI" channel is represented as a registered
surrogate support channel derived from safe external coupling raw/APC scores plus
sequence-only chemistry agreement; the certificate records this explicitly rather
than pretending a new MSA calculation happened.
"""

import hashlib
from dataclasses import asdict, dataclass
from math import exp
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.folding_coarse_grain_md_geometry import (
    COARSE_GRAIN_MD_CLAIM_RULE,
    COARSE_GRAIN_MD_RULE,
    MULTISTART_DCA_MD_MODE,
    CoarseGrainMDGeometryPacket,
    run_coarse_grain_md_geometry_packet,
)
from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_native_contact_eval import ContactPair, normalized_contact_pairs
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    MIN_SEQUENCE_SEPARATION,
    RealCoordinateVisualRow,
)
from pharmacotopology.folding_sequence_physical_priors import (
    build_sequence_physical_prior_scores,
    predict_lightweight_secondary_structure,
    secondary_structure_pair_score,
)

DENSE_ANCHOR_GEODESIC_MD_KIND = "dense_anchor_geodesic_restrained_md_geometry_v0"
DENSE_ANCHOR_GEODESIC_MODE = "external_dense_anchor_geodesic_md_geometry"
DENSE_ANCHOR_GEODESIC_RULE = (
    COARSE_GRAIN_MD_RULE
    + ";dense_anchor_expansion_from_safe_dca_raw_apc_sequence_closure_ss;"
    + "anchor_chaining;geodesic_interpolation;weak_geodesic_restraints"
)
DENSE_ANCHOR_GEODESIC_CLAIM_RULE = (
    COARSE_GRAIN_MD_CLAIM_RULE
    + ";dense_anchor_claim_requires_precision_recall_long_range_ge_0_70;"
    + "dense_anchors_may_support_external_dca_geometry_claim_but_not_universal_physics"
)


@dataclass(frozen=True)
class DenseAnchorRecord:
    row_id: str
    source_accession: str
    i: int
    j: int
    confidence: float
    channel: str
    parent_a: str = ""
    parent_b: str = ""
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class DenseAnchorRowSummary:
    row_id: str
    source_accession: str
    sequence_length: int
    direct_dca_anchor_count: int
    mi_surrogate_anchor_count: int
    sequence_closure_anchor_count: int
    beta_bridge_anchor_count: int
    geodesic_anchor_count: int
    total_dense_anchor_count: int
    dense_anchor_map_hash: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class DenseAnchorGeodesicPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    raw_msa_available_for_true_mi: bool
    mi_channel_is_surrogate_not_new_msa_calculation: bool
    direct_global_structure_generation_included: bool
    contacts_predicted_before_structure: bool
    dense_anchor_expansion_included: bool
    anchor_chaining_included: bool
    geodesic_interpolation_included: bool
    atomistic_md_engine_used: bool
    dependency_free_md_style_relaxation_used: bool
    mean_dense_anchor_count: float
    min_dense_anchor_count: int
    md_packet: CoarseGrainMDGeometryPacket
    anchor_rows: tuple[DenseAnchorRowSummary, ...]
    anchors: tuple[DenseAnchorRecord, ...]
    dense_anchor_geodesic_claim_allowed: bool
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
            "raw_msa_available_for_true_mi": self.raw_msa_available_for_true_mi,
            "mi_channel_is_surrogate_not_new_msa_calculation": self.mi_channel_is_surrogate_not_new_msa_calculation,
            "direct_global_structure_generation_included": self.direct_global_structure_generation_included,
            "contacts_predicted_before_structure": self.contacts_predicted_before_structure,
            "dense_anchor_expansion_included": self.dense_anchor_expansion_included,
            "anchor_chaining_included": self.anchor_chaining_included,
            "geodesic_interpolation_included": self.geodesic_interpolation_included,
            "atomistic_md_engine_used": self.atomistic_md_engine_used,
            "dependency_free_md_style_relaxation_used": self.dependency_free_md_style_relaxation_used,
            "mean_dense_anchor_count": self.mean_dense_anchor_count,
            "min_dense_anchor_count": self.min_dense_anchor_count,
            "md_packet": self.md_packet.to_dict(),
            "anchor_rows": [row.to_dict() for row in self.anchor_rows],
            "anchors": [anchor.to_dict() for anchor in self.anchors],
            "dense_anchor_geodesic_claim_allowed": self.dense_anchor_geodesic_claim_allowed,
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


def _safe_direct_constraints(
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    *,
    max_direct: int,
) -> tuple[CouplingConstraint, ...]:
    safe = [
        c
        for c in constraints
        if (c.row_id == row.row_id or c.source_accession == row.source_accession)
        and not c.coordinate_truth_used_to_build_constraint
        and not c.native_truth_used_before_coupling_selection
        and not c.structure_model_used
        and c.i < c.j
        and c.sequence_separation >= MIN_SEQUENCE_SEPARATION
    ]
    safe.sort(key=lambda c: (-float(c.confidence), c.rank if c.rank else 999999, c.i, c.j))
    return tuple(safe[: max(1, max_direct)])


def _local_dca_support(pair: ContactPair, direct: Sequence[CouplingConstraint]) -> float:
    if not direct:
        return 0.0
    left, right = pair
    best = 0.0
    for c in direct:
        dist = abs(left - c.i) + abs(right - c.j)
        support = float(c.confidence) * exp(-dist / 8.0)
        if support > best:
            best = support
    return _bounded01(best)


def _constraint(
    *,
    row: RealCoordinateVisualRow,
    pair: ContactPair,
    confidence: float,
    channel: str,
    rank: int,
) -> CouplingConstraint:
    i, j = pair
    sep = j - i
    return CouplingConstraint(
        row_id=row.row_id,
        source_accession=row.source_accession,
        constraint_id=f"densegeo_{row.row_id}_{i}_{j}_{channel}_{rank}",
        i=i,
        j=j,
        sequence_separation=sep,
        normalized_separation=_bounded01(sep / max(1, row.sequence_length)),
        confidence=_bounded01(confidence),
        constraint_class=f"dense_anchor_{channel}",
        source_kind="native_free_dense_anchor_geodesic_v0",
        coordinate_truth_used_to_build_constraint=False,
        native_truth_used_before_coupling_selection=False,
        structure_model_used=False,
        raw_sequence_exposed=False,
        raw_score=_rounded(confidence),
        apc_corrected_score=_rounded(confidence),
        rank=rank,
        rank_fraction=_bounded01(rank / max(1, row.sequence_length * 2)),
        msa_source_kind="external_safe_dca_plus_sequence_only_surrogate",
    )


def _records_from_constraints(
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    channel: str,
) -> list[DenseAnchorRecord]:
    return [
        DenseAnchorRecord(
            row_id=row.row_id,
            source_accession=row.source_accession,
            i=c.i,
            j=c.j,
            confidence=_bounded01(c.confidence),
            channel=channel,
        )
        for c in constraints
    ]


def _sequence_closure_candidates(
    row: RealCoordinateVisualRow,
    direct: Sequence[CouplingConstraint],
    *,
    max_candidates: int,
) -> list[tuple[ContactPair, float, str]]:
    n = row.sequence_length
    candidate_pairs = [
        (i, j)
        for i in range(1, n + 1)
        for j in range(i + MIN_SEQUENCE_SEPARATION, n + 1)
        if j - i >= 8
    ]
    priors = build_sequence_physical_prior_scores(row=row, candidate_pairs=candidate_pairs)
    ss = predict_lightweight_secondary_structure(row.sequence)
    scored: list[tuple[float, ContactPair, str]] = []
    for pair in candidate_pairs:
        support = _local_dca_support(pair, direct)
        if support < 0.16:
            continue
        prior = priors.get(pair)
        if prior is None:
            continue
        ss_score = secondary_structure_pair_score(pair, ss)
        sep = pair[1] - pair[0]
        mi_surrogate = _bounded01(0.55 * support + 0.25 * prior.physical_prior_score + 0.20 * ss_score)
        if mi_surrogate >= 0.42:
            scored.append((mi_surrogate, pair, "mi_surrogate_sequence_closure"))
        elif ss[pair[0] - 1] == "E" and ss[pair[1] - 1] == "E" and sep >= 24 and support >= 0.10:
            scored.append((_bounded01(0.40 * support + 0.60 * ss_score), pair, "beta_bridge"))
    scored.sort(key=lambda item: (-item[0], item[1][0], item[1][1]))
    out: list[tuple[ContactPair, float, str]] = []
    seen: set[ContactPair] = set()
    for score, pair, channel in scored:
        if pair in seen:
            continue
        seen.add(pair)
        out.append((pair, score, channel))
        if len(out) >= max_candidates:
            break
    return out


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


def _geodesic_interpolation(a: ContactPair, b: ContactPair) -> list[ContactPair]:
    i1, j1 = a
    i2, j2 = b
    steps = max(abs(i2 - i1), abs(j2 - j1))
    if steps <= 1 or steps > 32:
        return []
    out: list[ContactPair] = []
    for t in range(1, steps):
        i = int(round(i1 + (i2 - i1) * t / steps))
        j = int(round(j1 + (j2 - j1) * t / steps))
        if j - i >= MIN_SEQUENCE_SEPARATION:
            out.append((min(i, j), max(i, j)))
    return normalized_contact_pairs(out)


def build_dense_anchor_geodesic_constraints(
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    *,
    max_direct: int | None = None,
    max_sequence_closure: int | None = None,
    max_geodesic: int | None = None,
) -> tuple[tuple[CouplingConstraint, ...], DenseAnchorRowSummary, tuple[DenseAnchorRecord, ...]]:
    n = row.sequence_length
    direct_limit = max_direct or max(48, min(n, 256))
    closure_limit = max_sequence_closure or max(72, min(n, 280))
    geodesic_limit = max_geodesic or max(96, min(n * 2, 520))
    direct = list(_safe_direct_constraints(row, constraints, max_direct=direct_limit))
    dense_by_pair: dict[ContactPair, CouplingConstraint] = {c.pair(): c for c in direct}
    records: list[DenseAnchorRecord] = _records_from_constraints(row, direct, "direct_dca")
    channel_counts = {"direct_dca": len(direct), "mi_surrogate_sequence_closure": 0, "beta_bridge": 0, "geodesic": 0}

    rank = len(direct) + 1
    for pair, score, channel in _sequence_closure_candidates(row, direct, max_candidates=closure_limit):
        if pair not in dense_by_pair:
            dense_by_pair[pair] = _constraint(row=row, pair=pair, confidence=score * 0.72, channel=channel, rank=rank)
            records.append(DenseAnchorRecord(row.row_id, row.source_accession, pair[0], pair[1], _bounded01(score * 0.72), channel))
            channel_counts[channel] = channel_counts.get(channel, 0) + 1
            rank += 1

    seed_pairs = list(dense_by_pair)
    chains = _anchor_chains(seed_pairs)
    geodesic_added = 0
    for chain in chains:
        if len(chain) < 2:
            continue
        for left, right in zip(chain[:-1], chain[1:]):
            for pair in _geodesic_interpolation(left, right):
                if pair in dense_by_pair:
                    continue
                if pair[0] < 1 or pair[1] > n:
                    continue
                support = max(_local_dca_support(pair, direct), 0.25)
                confidence = _bounded01(0.30 + 0.35 * support)
                dense_by_pair[pair] = _constraint(row=row, pair=pair, confidence=confidence, channel="geodesic", rank=rank)
                records.append(
                    DenseAnchorRecord(
                        row.row_id,
                        row.source_accession,
                        pair[0],
                        pair[1],
                        confidence,
                        "geodesic",
                        parent_a=f"{left[0]}-{left[1]}",
                        parent_b=f"{right[0]}-{right[1]}",
                    )
                )
                rank += 1
                geodesic_added += 1
                if geodesic_added >= geodesic_limit:
                    break
            if geodesic_added >= geodesic_limit:
                break
        if geodesic_added >= geodesic_limit:
            break
    channel_counts["geodesic"] = geodesic_added

    dense = tuple(sorted(dense_by_pair.values(), key=lambda c: (-c.confidence, c.i, c.j)))
    summary = DenseAnchorRowSummary(
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_length=n,
        direct_dca_anchor_count=channel_counts.get("direct_dca", 0),
        mi_surrogate_anchor_count=channel_counts.get("mi_surrogate_sequence_closure", 0),
        sequence_closure_anchor_count=channel_counts.get("mi_surrogate_sequence_closure", 0),
        beta_bridge_anchor_count=channel_counts.get("beta_bridge", 0),
        geodesic_anchor_count=channel_counts.get("geodesic", 0),
        total_dense_anchor_count=len(dense),
        dense_anchor_map_hash=_pair_hash([c.pair() for c in dense]),
    )
    return dense, summary, tuple(records)


def run_dense_anchor_geodesic_packet(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    constraints: Sequence[CouplingConstraint],
    md_steps: int = 96,
    restarts: int = 2,
    max_direct: int | None = None,
    max_sequence_closure: int | None = None,
    max_geodesic: int | None = None,
) -> DenseAnchorGeodesicPacket:
    dense_constraints: list[CouplingConstraint] = []
    row_summaries: list[DenseAnchorRowSummary] = []
    anchor_records: list[DenseAnchorRecord] = []
    for row in rows:
        dense, summary, records = build_dense_anchor_geodesic_constraints(
            row,
            constraints,
            max_direct=max_direct,
            max_sequence_closure=max_sequence_closure,
            max_geodesic=max_geodesic,
        )
        dense_constraints.extend(dense)
        row_summaries.append(summary)
        anchor_records.extend(records[:900])

    md_packet = run_coarse_grain_md_geometry_packet(
        rows,
        constraints=tuple(dense_constraints),
        source_mode=DENSE_ANCHOR_GEODESIC_MODE,
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
        rejection = "dense_anchor_geodesic_md_survived_external_geometry_gate_not_universal_physics"
    elif md_packet.mean_native_contact_precision_after_audit < 0.70:
        rejection = "dense_anchor_geodesic_claim_rejected_precision_below_0_70"
    elif md_packet.mean_native_contact_recall_after_audit < 0.70:
        rejection = "dense_anchor_geodesic_claim_rejected_recall_below_0_70"
    elif md_packet.mean_long_range_contact_recall_after_audit < 0.70:
        rejection = "dense_anchor_geodesic_claim_rejected_long_range_recall_below_0_70"
    else:
        rejection = "dense_anchor_geodesic_claim_rejected_control_or_row_gate"
    counts = [row.total_dense_anchor_count for row in row_summaries]
    return DenseAnchorGeodesicPacket(
        kind=DENSE_ANCHOR_GEODESIC_MD_KIND,
        source_mode=DENSE_ANCHOR_GEODESIC_MODE,
        row_count=len(rows),
        decision_rule=DENSE_ANCHOR_GEODESIC_RULE,
        claim_rule=DENSE_ANCHOR_GEODESIC_CLAIM_RULE,
        raw_msa_available_for_true_mi=False,
        mi_channel_is_surrogate_not_new_msa_calculation=True,
        direct_global_structure_generation_included=True,
        contacts_predicted_before_structure=False,
        dense_anchor_expansion_included=True,
        anchor_chaining_included=True,
        geodesic_interpolation_included=True,
        atomistic_md_engine_used=False,
        dependency_free_md_style_relaxation_used=True,
        mean_dense_anchor_count=_rounded(mean(counts)) if counts else 0.0,
        min_dense_anchor_count=min(counts) if counts else 0,
        md_packet=md_packet,
        anchor_rows=tuple(row_summaries),
        anchors=tuple(anchor_records),
        dense_anchor_geodesic_claim_allowed=claim,
        universal_physical_law_claim_allowed=False,
        folding_problem_solved=claim,
        claim_rejection_reason=rejection,
        coordinate_truth_used_before_selection=md_packet.coordinate_truth_used_before_selection,
        native_truth_used_before_selection=md_packet.native_truth_used_before_selection,
        structure_model_used_before_selection=md_packet.structure_model_used_before_selection,
        learned_geometry_prior_used_before_selection=md_packet.learned_geometry_prior_used_before_selection,
    )
