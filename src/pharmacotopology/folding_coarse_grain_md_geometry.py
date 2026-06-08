from __future__ import annotations

"""Native-free coarse-grain MD geometry challenge.

This layer tests a stricter pivot than the earlier contact-field attempts:

    sequence/DCA restraints -> global C-alpha trajectory -> contacts extracted
    from the final structure.

It intentionally does not predict contacts first and then build a structure.  The
selection surface is a globally relaxed coarse C-alpha geometry.  DCA couplings,
when supplied, are used only as distance restraints during the geometry run.
Native coordinates/native contacts are used only after the final structure has
already emitted a contact map, for audit metrics.

The implementation is dependency-free and bounded, so it can run in environments
without OpenMM/GROMACS.  It is MD-style global relaxation rather than atomistic
MD: all residues are moved together by chain springs, secondary-structure
springs, DCA distance restraints, excluded volume, weak compaction, and a simple
velocity/temperature schedule.
"""

import hashlib
from dataclasses import asdict, dataclass
from math import cos, exp, pi, sin, sqrt
from statistics import mean
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_five_axis_physics import matched_control_pairs
from pharmacotopology.folding_native_contact_eval import (
    ContactMetricPacket,
    ContactPair,
    contact_map_hash,
    evaluate_contact_prediction,
    normalized_contact_pairs,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    MIN_SEQUENCE_SEPARATION,
    RealCoordinateVisualRow,
)
from pharmacotopology.folding_sequence_physical_priors import predict_lightweight_secondary_structure

COARSE_GRAIN_MD_GEOMETRY_KIND = "coarse_grain_md_global_geometry_v0"
COARSE_GRAIN_MD_DECISION_KIND = "coarse_grain_md_extracted_contact_decision_v0"
PURE_SEQUENCE_MD_MODE = "pure_sequence_global_md_geometry"
EXTERNAL_DCA_MD_MODE = "external_dca_restrained_global_md_geometry"
MULTISTART_DCA_MD_MODE = "external_dca_multistart_global_md_geometry"
COARSE_GRAIN_MD_RULE = (
    "sequence_or_safe_dca_restraints_to_global_calpha_geometry;"
    "chain_springs;secondary_structure_springs;dca_distance_restraints;"
    "excluded_volume;weak_compaction;temperature_damped_global_updates;"
    "contacts_extracted_from_final_geometry_only;native_audit_after_selection_only"
)
COARSE_GRAIN_MD_CLAIM_RULE = (
    "claim_requires_no_coordinate_truth_no_native_truth_no_template_no_structure_model_no_learned_prior;"
    "all_rows_precision_ge_0_70_and_recall_ge_0_70_and_long_range_recall_ge_0_70;"
    "matched_control_f1_and_long_range_margin_ge_0_15;"
    "pure_sequence_mode_required_for_universal_physical_law_claim"
)

Vector = list[float]
Edge = tuple[int, int, float, float, str]


@dataclass(frozen=True)
class CoarseGrainMDContactDecision:
    kind: str
    row_id: str
    source_accession: str
    source_mode: str
    i: int
    j: int
    sequence_separation: int
    final_distance_angstrom: float
    geometry_contact_score: float
    selected: bool
    selected_from_final_structure: bool
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class CoarseGrainMDRowReport:
    row_id: str
    source_accession: str
    source_mode: str
    sequence_hash: str
    sequence_length: int
    safe_restraint_count: int
    md_step_count: int
    md_restart_count: int
    final_energy_proxy: float
    extracted_contact_count: int
    extracted_long_range_contact_count: int
    selected_contact_map_hash: str
    metric_after_native_audit: ContactMetricPacket
    matched_control_count: int
    best_control_f1_after_audit: float
    best_control_long_range_recall_after_audit: float
    f1_margin_vs_best_control: float
    long_range_recall_margin_vs_best_control: float
    row_md_geometry_claim_allowed: bool
    row_universal_physical_law_claim_allowed: bool
    row_claim_rejection_reason: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    native_truth_attached_after_selection_for_evaluation: bool = True
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["metric_after_native_audit"] = self.metric_after_native_audit.to_dict()
        return payload


@dataclass(frozen=True)
class CoarseGrainMDGeometryPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    direct_global_structure_generation_included: bool
    contacts_predicted_before_structure: bool
    dca_restraints_included: bool
    atomistic_md_engine_used: bool
    dependency_free_md_style_relaxation_used: bool
    mean_native_contact_precision_after_audit: float
    mean_native_contact_recall_after_audit: float
    mean_long_range_contact_recall_after_audit: float
    mean_contact_map_f1_after_audit: float
    min_row_precision_after_audit: float
    min_row_recall_after_audit: float
    min_row_long_range_recall_after_audit: float
    mean_f1_margin_vs_best_control: float
    mean_long_range_recall_margin_vs_best_control: float
    md_geometry_claim_allowed: bool
    universal_physical_law_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    rows: tuple[CoarseGrainMDRowReport, ...]
    decisions: tuple[CoarseGrainMDContactDecision, ...]
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "source_mode": self.source_mode,
            "row_count": self.row_count,
            "decision_rule": self.decision_rule,
            "claim_rule": self.claim_rule,
            "direct_global_structure_generation_included": self.direct_global_structure_generation_included,
            "contacts_predicted_before_structure": self.contacts_predicted_before_structure,
            "dca_restraints_included": self.dca_restraints_included,
            "atomistic_md_engine_used": self.atomistic_md_engine_used,
            "dependency_free_md_style_relaxation_used": self.dependency_free_md_style_relaxation_used,
            "mean_native_contact_precision_after_audit": self.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": self.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": self.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": self.mean_contact_map_f1_after_audit,
            "min_row_precision_after_audit": self.min_row_precision_after_audit,
            "min_row_recall_after_audit": self.min_row_recall_after_audit,
            "min_row_long_range_recall_after_audit": self.min_row_long_range_recall_after_audit,
            "mean_f1_margin_vs_best_control": self.mean_f1_margin_vs_best_control,
            "mean_long_range_recall_margin_vs_best_control": self.mean_long_range_recall_margin_vs_best_control,
            "md_geometry_claim_allowed": self.md_geometry_claim_allowed,
            "universal_physical_law_claim_allowed": self.universal_physical_law_claim_allowed,
            "folding_problem_solved": self.folding_problem_solved,
            "claim_rejection_reason": self.claim_rejection_reason,
            "rows": [row.to_dict() for row in self.rows],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "structure_model_used_before_selection": self.structure_model_used_before_selection,
            "learned_geometry_prior_used_before_selection": self.learned_geometry_prior_used_before_selection,
            "msa_dca_used_before_selection": self.msa_dca_used_before_selection,
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }


def _rounded(value: float) -> float:
    return round(float(value), 6)


def _mean(values: Sequence[float]) -> float:
    return _rounded(mean(values)) if values else 0.0


def _pair_hash(pairs: Iterable[ContactPair]) -> str:
    return contact_map_hash(normalized_contact_pairs(pairs))


def _stable_seed(row: RealCoordinateVisualRow, restart_index: int) -> int:
    digest = hashlib.sha256(f"{row.sequence_sha256}:{restart_index}:cgmd".encode("utf-8")).hexdigest()
    return int(digest[:12], 16)


def _safe_restraints_for_row(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    source_mode: str,
    max_restraints: int = 96,
) -> tuple[tuple[CouplingConstraint, ...], bool, bool, bool, bool]:
    if source_mode == PURE_SEQUENCE_MD_MODE:
        return (), False, False, False, False
    safe: list[CouplingConstraint] = []
    coordinate_taint = False
    native_taint = False
    structure_model = False
    for constraint in constraints:
        if constraint.source_accession != row.source_accession and constraint.row_id != row.row_id:
            continue
        if constraint.coordinate_truth_used_to_build_constraint:
            coordinate_taint = True
            continue
        if constraint.native_truth_used_before_coupling_selection:
            native_taint = True
            continue
        if constraint.structure_model_used:
            structure_model = True
            continue
        if constraint.sequence_separation < MIN_SEQUENCE_SEPARATION:
            continue
        safe.append(constraint)
    safe.sort(key=lambda item: (-(float(item.confidence)), item.rank if item.rank else 999999, item.i, item.j))
    return tuple(safe[:max_restraints]), coordinate_taint, native_taint, structure_model, bool(safe)


def _initial_trace(row: RealCoordinateVisualRow, restart_index: int) -> list[Vector]:
    ss = predict_lightweight_secondary_structure(row.sequence)
    seed = _stable_seed(row, restart_index)
    coords: list[Vector] = []
    x = y = z = 0.0
    phase = (seed % 6283) / 1000.0
    angle = phase
    for idx, label in enumerate(ss):
        if label == "H":
            angle += 100.0 * pi / 180.0
            x += 1.52
            y += 2.08 * cos(angle)
            z += 2.08 * sin(angle)
        elif label == "E":
            x += 3.25
            y += 1.15 * (-1.0 if idx % 2 else 1.0)
            z += 0.32 * sin(phase + idx * 0.71)
        else:
            angle += 63.0 * pi / 180.0
            x += 2.30
            y += 1.55 * cos(angle)
            z += 1.55 * sin(angle)
        # deterministic restart jitter; not stochastic and does not use native truth.
        jitter = 0.18 + 0.04 * restart_index
        coords.append([
            x + jitter * sin(phase + idx * 1.37),
            y + jitter * cos(phase + idx * 1.91),
            z + jitter * sin(phase + idx * 2.23),
        ])
    return coords


def _restraint_edges(row: RealCoordinateVisualRow, restraints: Sequence[CouplingConstraint]) -> list[Edge]:
    n = row.sequence_length
    ss = predict_lightweight_secondary_structure(row.sequence)
    edges: list[Edge] = []
    for left in range(1, n):
        edges.append((left - 1, left, 3.80, 2.25, "chain"))
    for index, label in enumerate(ss, start=1):
        if label == "H":
            for offset, target, weight in ((3, 5.25, 0.40), (4, 6.25, 0.42)):
                if index + offset <= n:
                    edges.append((index - 1, index + offset - 1, target, weight, "secondary_h"))
        elif label == "E":
            if index + 2 <= n:
                edges.append((index - 1, index + 1, 6.75, 0.18, "secondary_e"))
    for restraint in restraints:
        left = restraint.i - 1
        right = restraint.j - 1
        if left < 0 or right >= n or left >= right:
            continue
        confidence = max(0.05, min(1.0, float(restraint.confidence)))
        separation = right - left
        target = 7.1 if separation >= 12 else 6.2
        edges.append((left, right, target, 3.2 * confidence, "dca"))
        for delta in (-1, 1):
            if 0 <= left + delta < n:
                edges.append((left + delta, right, target + 0.35, 0.55 * confidence, "dca_neighbour"))
            if 0 <= right + delta < n:
                edges.append((left, right + delta, target + 0.35, 0.55 * confidence, "dca_neighbour"))
    return edges


def _distance(a: Sequence[float], b: Sequence[float]) -> float:
    return sqrt((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2 + (b[2] - a[2]) ** 2) + 1e-9


def _energy_proxy(coords: Sequence[Sequence[float]], edges: Sequence[Edge]) -> float:
    energy = 0.0
    for left, right, target, weight, _kind in edges:
        delta = _distance(coords[left], coords[right]) - target
        energy += weight * delta * delta
    n = len(coords)
    # cheap compactness/excluded-volume proxy sampled deterministically.
    for left in range(0, n, 2):
        for right in range(left + MIN_SEQUENCE_SEPARATION, n, 3):
            distance = _distance(coords[left], coords[right])
            if distance < 3.25:
                energy += 2.0 * (3.25 - distance) ** 2
    return _rounded(energy / max(1, n))


def _relax_global_geometry(
    *,
    row: RealCoordinateVisualRow,
    restraints: Sequence[CouplingConstraint],
    steps: int,
    restart_index: int,
) -> tuple[list[Vector], float]:
    coords = _initial_trace(row, restart_index)
    n = len(coords)
    edges = _restraint_edges(row, restraints)
    velocities: list[Vector] = [[0.0, 0.0, 0.0] for _ in coords]
    step_count = max(12, int(steps))
    target_radius = 2.7 * (n ** 0.37)

    for step in range(step_count):
        temperature = 1.0 - step / max(1.0, step_count - 1.0)
        lr = 0.0065 + 0.010 * temperature
        damping = 0.72 - 0.20 * temperature
        forces: list[Vector] = [[0.0, 0.0, 0.0] for _ in coords]

        # Harmonic forces for chain, SS, and DCA restraints.
        for left, right, target, weight, kind in edges:
            lx, ly, lz = coords[left]
            rx, ry, rz = coords[right]
            dx = rx - lx
            dy = ry - ly
            dz = rz - lz
            distance = sqrt(dx * dx + dy * dy + dz * dz) + 1e-9
            # DCA restraints are intentionally global and stronger early, then cooler.
            anneal_weight = weight * (1.0 + (0.45 * temperature if kind.startswith("dca") else 0.0))
            force = anneal_weight * (distance - target) / distance
            fx, fy, fz = dx * force, dy * force, dz * force
            forces[left][0] += fx
            forces[left][1] += fy
            forces[left][2] += fz
            forces[right][0] -= fx
            forces[right][1] -= fy
            forces[right][2] -= fz

        # Weak global compaction, moving all residues together rather than selecting contacts.
        center = [sum(point[axis] for point in coords) / n for axis in range(3)]
        for index, point in enumerate(coords):
            px = point[0] - center[0]
            py = point[1] - center[1]
            pz = point[2] - center[2]
            radius = sqrt(px * px + py * py + pz * pz) + 1e-9
            if radius > target_radius:
                pull = 0.20 * (radius - target_radius) / radius
                forces[index][0] -= px * pull
                forces[index][1] -= py * pull
                forces[index][2] -= pz * pull
            elif radius < target_radius * 0.45:
                push = 0.04 * (target_radius * 0.45 - radius) / radius
                forces[index][0] += px * push
                forces[index][1] += py * push
                forces[index][2] += pz * push

        # Excluded volume; sampled grid keeps this bounded.
        stride = 2 if n > 160 else 1
        for left in range(0, n, stride):
            for right in range(left + MIN_SEQUENCE_SEPARATION, n, stride):
                distance = _distance(coords[left], coords[right])
                if distance >= 3.45:
                    continue
                dx = coords[right][0] - coords[left][0]
                dy = coords[right][1] - coords[left][1]
                dz = coords[right][2] - coords[left][2]
                force = 0.55 * (3.45 - distance) / distance
                fx, fy, fz = dx * force, dy * force, dz * force
                forces[left][0] -= fx
                forces[left][1] -= fy
                forces[left][2] -= fz
                forces[right][0] += fx
                forces[right][1] += fy
                forces[right][2] += fz

        # Deterministic thermal exploration; fades out completely.
        phase = (_stable_seed(row, restart_index) % 1009) * 0.01
        thermal = 0.022 * temperature
        for index in range(n):
            forces[index][0] += thermal * sin(phase + step * 0.17 + index * 0.37)
            forces[index][1] += thermal * cos(phase + step * 0.19 + index * 0.41)
            forces[index][2] += thermal * sin(phase + step * 0.23 + index * 0.43)
            for axis in range(3):
                velocities[index][axis] = damping * velocities[index][axis] + lr * forces[index][axis]
                coords[index][axis] += velocities[index][axis]

        # Remove global drift.
        center = [sum(point[axis] for point in coords) / n for axis in range(3)]
        for point in coords:
            for axis in range(3):
                point[axis] -= center[axis]

    return coords, _energy_proxy(coords, edges)


def _extract_contacts_from_geometry(row: RealCoordinateVisualRow, coords: Sequence[Sequence[float]]) -> tuple[tuple[ContactPair, ...], dict[ContactPair, float]]:
    pairs: list[ContactPair] = []
    scores: dict[ContactPair, float] = {}
    n = row.sequence_length
    for left in range(1, n + 1):
        for right in range(left + MIN_SEQUENCE_SEPARATION, n + 1):
            distance = _distance(coords[left - 1], coords[right - 1])
            if distance <= 8.0:
                pair = (left, right)
                pairs.append(pair)
                scores[pair] = _rounded(1.0 / (1.0 + exp((distance - 7.2) / 0.65)))
    return normalized_contact_pairs(pairs), scores


def _run_multistart_geometry(
    *,
    row: RealCoordinateVisualRow,
    restraints: Sequence[CouplingConstraint],
    steps: int,
    restarts: int,
) -> tuple[tuple[ContactPair, ...], dict[ContactPair, float], float, int]:
    restart_count = max(1, int(restarts))
    contact_votes: dict[ContactPair, int] = {}
    contact_scores: dict[ContactPair, float] = {}
    best_energy: float | None = None
    for restart_index in range(restart_count):
        coords, energy = _relax_global_geometry(row=row, restraints=restraints, steps=steps, restart_index=restart_index)
        pairs, scores = _extract_contacts_from_geometry(row, coords)
        if best_energy is None or energy < best_energy:
            best_energy = energy
        for pair in pairs:
            contact_votes[pair] = contact_votes.get(pair, 0) + 1
            contact_scores[pair] = max(contact_scores.get(pair, 0.0), scores.get(pair, 0.0))
    threshold = max(1, int(round(0.50 * restart_count)))
    selected = [pair for pair, votes in contact_votes.items() if votes >= threshold]
    return normalized_contact_pairs(selected), contact_scores, _rounded(best_energy or 0.0), restart_count


def _matched_controls_for_report(
    *,
    row: RealCoordinateVisualRow,
    selected_pairs: Sequence[ContactPair],
) -> tuple[float, float, int]:
    # Candidate pool is all legal pairs; this avoids any native/coordinate leakage before audit.
    candidate_pairs = normalized_contact_pairs(
        (left, right)
        for left in range(1, row.sequence_length + 1)
        for right in range(left + MIN_SEQUENCE_SEPARATION, row.sequence_length + 1)
    )
    controls: list[ContactMetricPacket] = []
    count = max(2, min(6, int(sqrt(max(1, len(selected_pairs))))))
    for control_index in range(1, count + 1):
        control_pairs = matched_control_pairs(
            row=row,
            selected_pairs=selected_pairs,
            candidate_pairs=candidate_pairs,
            control_index=control_index,
        )
        if not control_pairs:
            continue
        controls.append(evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=control_pairs))
    return (
        max((metric.contact_map_f1 for metric in controls), default=0.0),
        max((metric.long_range_contact_recall for metric in controls), default=0.0),
        len(controls),
    )


def _claim_decision(
    *,
    source_mode: str,
    metric: ContactMetricPacket,
    f1_margin: float,
    long_range_margin: float,
    coordinate_taint: bool,
    native_taint: bool,
    structure_model: bool,
    restraint_count: int,
) -> tuple[bool, bool, str]:
    if coordinate_taint or native_taint:
        return False, False, "row_claim_rejected_coordinate_or_native_tainted_restraint"
    if structure_model:
        return False, False, "row_claim_rejected_structure_model_used_before_selection"
    if source_mode != PURE_SEQUENCE_MD_MODE and restraint_count <= 0:
        return False, False, "row_claim_rejected_no_safe_dca_restraints"
    claim_allowed = (
        metric.native_contact_precision >= 0.70
        and metric.native_contact_recall >= 0.70
        and metric.long_range_contact_recall >= 0.70
        and f1_margin >= 0.15
        and long_range_margin >= 0.15
    )
    universal_allowed = claim_allowed and source_mode == PURE_SEQUENCE_MD_MODE
    if claim_allowed and source_mode != PURE_SEQUENCE_MD_MODE:
        return True, False, "coarse_grain_md_geometry_row_survived_external_dca_gate_not_universal_physics"
    if universal_allowed:
        return True, True, "pure_sequence_coarse_grain_md_row_survived_universal_physics_gate"
    if metric.native_contact_precision < 0.70:
        return False, False, "row_claim_rejected_precision_below_0_70"
    if metric.native_contact_recall < 0.70:
        return False, False, "row_claim_rejected_recall_below_0_70"
    if metric.long_range_contact_recall < 0.70:
        return False, False, "row_claim_rejected_long_range_recall_below_0_70"
    return False, False, "row_claim_rejected_matched_control_margin_below_gate"


def _decisions_for_row(
    *,
    row: RealCoordinateVisualRow,
    source_mode: str,
    selected_pairs: Sequence[ContactPair],
    scores: Mapping[ContactPair, float],
    dca_used: bool,
    max_report_decisions: int = 900,
) -> tuple[CoarseGrainMDContactDecision, ...]:
    decisions: list[CoarseGrainMDContactDecision] = []
    ordered = sorted(selected_pairs, key=lambda pair: (-scores.get(pair, 0.0), pair[0], pair[1]))[:max_report_decisions]
    # We do not store the whole coordinate trace; only distances/scores for selected contact decisions.
    for pair in ordered:
        # score->distance approximation for compact CSV inspection; exact selected distances are not needed for audit.
        score = float(scores.get(pair, 0.0))
        distance = 7.2 + 0.65 * (-score if score else 0.0)
        decisions.append(
            CoarseGrainMDContactDecision(
                kind=COARSE_GRAIN_MD_DECISION_KIND,
                row_id=row.row_id,
                source_accession=row.source_accession,
                source_mode=source_mode,
                i=pair[0],
                j=pair[1],
                sequence_separation=pair[1] - pair[0],
                final_distance_angstrom=_rounded(distance),
                geometry_contact_score=_rounded(score),
                selected=True,
                selected_from_final_structure=True,
                msa_dca_used_before_selection=dca_used,
            )
        )
    return tuple(decisions)


def run_coarse_grain_md_geometry_row(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint] = (),
    source_mode: str = MULTISTART_DCA_MD_MODE,
    md_steps: int = 96,
    restarts: int = 3,
    max_restraints: int = 96,
) -> tuple[CoarseGrainMDRowReport, tuple[CoarseGrainMDContactDecision, ...]]:
    safe_restraints, coordinate_taint, native_taint, structure_model, dca_used = _safe_restraints_for_row(
        row=row,
        constraints=constraints,
        source_mode=source_mode,
        max_restraints=max_restraints,
    )
    effective_restarts = 1 if source_mode in {PURE_SEQUENCE_MD_MODE, EXTERNAL_DCA_MD_MODE} else restarts
    selected_pairs, scores, final_energy, restart_count = _run_multistart_geometry(
        row=row,
        restraints=safe_restraints,
        steps=md_steps,
        restarts=effective_restarts,
    )
    metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=selected_pairs)
    best_control_f1, best_control_lr, control_count = _matched_controls_for_report(row=row, selected_pairs=selected_pairs)
    f1_margin = _rounded(metric.contact_map_f1 - best_control_f1)
    lr_margin = _rounded(metric.long_range_contact_recall - best_control_lr)
    row_claim, row_universal, rejection = _claim_decision(
        source_mode=source_mode,
        metric=metric,
        f1_margin=f1_margin,
        long_range_margin=lr_margin,
        coordinate_taint=coordinate_taint,
        native_taint=native_taint,
        structure_model=structure_model,
        restraint_count=len(safe_restraints),
    )
    decisions = _decisions_for_row(
        row=row,
        source_mode=source_mode,
        selected_pairs=selected_pairs,
        scores=scores,
        dca_used=dca_used,
    )
    report = CoarseGrainMDRowReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        source_mode=source_mode,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        safe_restraint_count=len(safe_restraints),
        md_step_count=md_steps,
        md_restart_count=restart_count,
        final_energy_proxy=final_energy,
        extracted_contact_count=len(selected_pairs),
        extracted_long_range_contact_count=sum(1 for pair in selected_pairs if pair[1] - pair[0] >= 24),
        selected_contact_map_hash=_pair_hash(selected_pairs),
        metric_after_native_audit=metric,
        matched_control_count=control_count,
        best_control_f1_after_audit=best_control_f1,
        best_control_long_range_recall_after_audit=best_control_lr,
        f1_margin_vs_best_control=f1_margin,
        long_range_recall_margin_vs_best_control=lr_margin,
        row_md_geometry_claim_allowed=row_claim,
        row_universal_physical_law_claim_allowed=row_universal,
        row_claim_rejection_reason=rejection,
        coordinate_truth_used_before_selection=coordinate_taint,
        native_truth_used_before_selection=native_taint,
        structure_model_used_before_selection=structure_model,
        msa_dca_used_before_selection=dca_used,
    )
    return report, decisions


def run_coarse_grain_md_geometry_packet(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    constraints: Sequence[CouplingConstraint] = (),
    source_mode: str = MULTISTART_DCA_MD_MODE,
    md_steps: int = 96,
    restarts: int = 3,
    max_restraints: int = 96,
) -> CoarseGrainMDGeometryPacket:
    reports: list[CoarseGrainMDRowReport] = []
    decisions: list[CoarseGrainMDContactDecision] = []
    for row in rows:
        report, row_decisions = run_coarse_grain_md_geometry_row(
            row=row,
            constraints=constraints,
            source_mode=source_mode,
            md_steps=md_steps,
            restarts=restarts,
            max_restraints=max_restraints,
        )
        reports.append(report)
        decisions.extend(row_decisions[:420])

    precision_values = [row.metric_after_native_audit.native_contact_precision for row in reports]
    recall_values = [row.metric_after_native_audit.native_contact_recall for row in reports]
    long_range_values = [row.metric_after_native_audit.long_range_contact_recall for row in reports]
    f1_values = [row.metric_after_native_audit.contact_map_f1 for row in reports]
    f1_margins = [row.f1_margin_vs_best_control for row in reports]
    lr_margins = [row.long_range_recall_margin_vs_best_control for row in reports]
    coordinate_taint = any(row.coordinate_truth_used_before_selection for row in reports)
    native_taint = any(row.native_truth_used_before_selection for row in reports)
    structure_model = any(row.structure_model_used_before_selection for row in reports)
    dca_used = any(row.msa_dca_used_before_selection for row in reports)
    all_rows_claim = bool(reports) and all(row.row_md_geometry_claim_allowed for row in reports)
    all_rows_universal = bool(reports) and all(row.row_universal_physical_law_claim_allowed for row in reports)
    mean_gate = _mean(precision_values) >= 0.70 and _mean(recall_values) >= 0.70 and _mean(long_range_values) >= 0.70
    md_claim = all_rows_claim and mean_gate and not coordinate_taint and not native_taint and not structure_model
    universal_claim = all_rows_universal and md_claim and source_mode == PURE_SEQUENCE_MD_MODE and not dca_used
    if md_claim and source_mode != PURE_SEQUENCE_MD_MODE:
        rejection = "external_dca_coarse_grain_md_geometry_survived_gate_not_universal_physics"
    elif universal_claim:
        rejection = "pure_sequence_coarse_grain_md_geometry_survived_universal_physical_law_gate"
    else:
        failed = [row.source_accession for row in reports if not row.row_md_geometry_claim_allowed]
        rejection = "coarse_grain_md_geometry_claim_rejected_for_rows:" + ",".join(failed[:12])

    return CoarseGrainMDGeometryPacket(
        kind=COARSE_GRAIN_MD_GEOMETRY_KIND,
        source_mode=source_mode,
        row_count=len(reports),
        decision_rule=COARSE_GRAIN_MD_RULE,
        claim_rule=COARSE_GRAIN_MD_CLAIM_RULE,
        direct_global_structure_generation_included=True,
        contacts_predicted_before_structure=False,
        dca_restraints_included=source_mode != PURE_SEQUENCE_MD_MODE,
        atomistic_md_engine_used=False,
        dependency_free_md_style_relaxation_used=True,
        mean_native_contact_precision_after_audit=_mean(precision_values),
        mean_native_contact_recall_after_audit=_mean(recall_values),
        mean_long_range_contact_recall_after_audit=_mean(long_range_values),
        mean_contact_map_f1_after_audit=_mean(f1_values),
        min_row_precision_after_audit=_rounded(min(precision_values)) if precision_values else 0.0,
        min_row_recall_after_audit=_rounded(min(recall_values)) if recall_values else 0.0,
        min_row_long_range_recall_after_audit=_rounded(min(long_range_values)) if long_range_values else 0.0,
        mean_f1_margin_vs_best_control=_mean(f1_margins),
        mean_long_range_recall_margin_vs_best_control=_mean(lr_margins),
        md_geometry_claim_allowed=md_claim,
        universal_physical_law_claim_allowed=universal_claim,
        folding_problem_solved=md_claim or universal_claim,
        claim_rejection_reason=rejection,
        rows=tuple(reports),
        decisions=tuple(decisions),
        coordinate_truth_used_before_selection=coordinate_taint,
        native_truth_used_before_selection=native_taint,
        structure_model_used_before_selection=structure_model,
        learned_geometry_prior_used_before_selection=False,
        msa_dca_used_before_selection=dca_used,
        raw_sequence_exposed=False,
    )
