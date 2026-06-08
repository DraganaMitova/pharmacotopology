from __future__ import annotations

"""Native-free iterative distance-geometry diffusion for contact maps.

This module is the strongest no-AlphaFold/no-template fallback in the project. It
turns weak DCA anchors plus sequence-only physical priors into a coarse 3D
C-alpha geometry, then diffuses the contact map through repeated geometry
reconstruction. Native coordinates are never read inside the predictor; metrics
are attached only by the caller for audit.
"""

from collections import defaultdict
from dataclasses import asdict, dataclass
from math import exp, sqrt
from statistics import mean
from typing import Mapping, Sequence

import numpy as np

from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_native_contact_eval import ContactPair, contact_map_hash, normalized_contact_pairs
from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent
from pharmacotopology.folding_real_coordinate_visual_benchmark import CONTACT_CUTOFF_ANGSTROM, MIN_SEQUENCE_SEPARATION, RealCoordinateVisualRow
from pharmacotopology.folding_sequence_physical_priors import (
    build_sequence_physical_prior_scores,
    predict_lightweight_secondary_structure,
    summarize_physical_priors,
)


ITERATIVE_DISTANCE_GEOMETRY_DIFFUSION_KIND = "native_free_iterative_distance_geometry_diffusion_v0"
ITERATIVE_DISTANCE_GEOMETRY_DIFFUSION_RULE = (
    "sequence_chain_metric+external_dca_anchor_field+sequence_physics+event_region_prior;"
    "iterative_graph_shortest_path_mds;native_truth_after_audit_only;no_alphafold_no_templates"
)


@dataclass(frozen=True)
class DiffusionContactScore:
    row_id: str
    source_accession: str
    i: int
    j: int
    sequence_separation: int
    final_score: float
    coupling_score: float
    anchor_field_score: float
    event_region_score: float
    physical_prior_score: float
    secondary_structure_pair_score: float
    geometry_score: float
    embedded_distance_angstrom: float | None
    selected: bool
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class DiffusionIteration:
    iteration_index: int
    candidate_pair_count: int
    seed_anchor_pair_count: int
    attraction_pair_count: int
    embedded_contact_pair_count: int
    selected_pair_count: int
    selected_contact_map_hash: str
    mean_selected_score: float
    mean_selected_embedded_distance_angstrom: float
    mean_selected_physical_prior_score: float
    mean_selected_geometry_score: float
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class DiffusionPacket:
    kind: str
    row_id: str
    source_accession: str
    sequence_length: int
    candidate_pair_count: int
    seed_anchor_pair_count: int
    final_pair_count: int
    final_long_range_pair_count: int
    final_contact_map_hash: str
    iteration_count: int
    decision_rule: str
    selected_pairs: tuple[ContactPair, ...]
    iterations: tuple[DiffusionIteration, ...]
    scores: tuple[DiffusionContactScore, ...]
    mean_final_contact_energy_score: float
    mean_final_secondary_structure_score: float
    mean_final_degree_consistency_score: float
    mean_final_physical_prior_score: float
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    structure_template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_report_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "row_id": self.row_id,
            "source_accession": self.source_accession,
            "sequence_length": self.sequence_length,
            "candidate_pair_count": self.candidate_pair_count,
            "seed_anchor_pair_count": self.seed_anchor_pair_count,
            "final_pair_count": self.final_pair_count,
            "final_long_range_pair_count": self.final_long_range_pair_count,
            "final_contact_map_hash": self.final_contact_map_hash,
            "iteration_count": self.iteration_count,
            "decision_rule": self.decision_rule,
            "mean_final_contact_energy_score": self.mean_final_contact_energy_score,
            "mean_final_secondary_structure_score": self.mean_final_secondary_structure_score,
            "mean_final_degree_consistency_score": self.mean_final_degree_consistency_score,
            "mean_final_physical_prior_score": self.mean_final_physical_prior_score,
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "alphafold_used_before_selection": self.alphafold_used_before_selection,
            "structure_template_used_before_selection": self.structure_template_used_before_selection,
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }



def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, float(value))), 6)



def _score(value: float) -> float:
    return round(float(value), 6)



def _mean(values: Sequence[float]) -> float:
    return _score(mean(values)) if values else 0.0



def _all_candidate_pairs(row: RealCoordinateVisualRow, *, minimum_sequence_separation: int) -> tuple[ContactPair, ...]:
    return tuple(
        (i, j)
        for i in range(1, row.sequence_length + 1)
        for j in range(i + minimum_sequence_separation, row.sequence_length + 1)
    )



def _safe_coupling_scores(row: RealCoordinateVisualRow, constraints: Sequence[CouplingConstraint]) -> dict[ContactPair, float]:
    scores: dict[ContactPair, float] = {}
    for constraint in constraints:
        if constraint.row_id != row.row_id and constraint.source_accession != row.source_accession:
            continue
        if constraint.coordinate_truth_used_to_build_constraint:
            continue
        if constraint.native_truth_used_before_coupling_selection:
            continue
        if constraint.structure_model_used:
            continue
        pair = constraint.pair()
        if pair[0] < 1 or pair[1] > row.sequence_length or pair[1] - pair[0] < MIN_SEQUENCE_SEPARATION:
            continue
        scores[pair] = max(scores.get(pair, 0.0), float(constraint.confidence))
    return scores



def _event_scores(row: RealCoordinateVisualRow, events: Sequence[NucleusClosureEvent]) -> dict[ContactPair, float]:
    output: dict[ContactPair, float] = {}
    for event in events:
        if event.row_id != row.row_id and event.source_accession != row.source_accession:
            continue
        if event.native_truth_used_before_event_generation:
            continue
        for raw_pair in event.candidate_region_pairs():
            if raw_pair[0] < 1 or raw_pair[1] > row.sequence_length or raw_pair[1] - raw_pair[0] < MIN_SEQUENCE_SEPARATION:
                continue
            output[raw_pair] = max(output.get(raw_pair, 0.0), float(event.nucleus_score))
    return output



def _rank_percentiles(raw: Mapping[ContactPair, float]) -> dict[ContactPair, float]:
    positive = sorted({float(value) for value in raw.values() if float(value) > 0.0})
    if not positive:
        return {pair: 0.0 for pair in raw}
    denom = max(1, len(positive))
    ranks = {value: (index + 1) / denom for index, value in enumerate(positive)}
    return {pair: _rounded(ranks.get(float(value), 0.0)) if float(value) > 0.0 else 0.0 for pair, value in raw.items()}



def _anchor_field(candidate_pairs: Sequence[ContactPair], coupling_scores: Mapping[ContactPair, float], *, max_anchors: int = 96) -> dict[ContactPair, float]:
    anchors = tuple(sorted(coupling_scores.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))[:max_anchors])
    if not anchors:
        return {pair: 0.0 for pair in candidate_pairs}
    output: dict[ContactPair, float] = {}
    for pair in candidate_pairs:
        i, j = pair
        best = 0.0
        for (ai, aj), confidence in anchors:
            # Support diffuses along the contact-map grid and also along the
            # anti-diagonal registry, which helps beta/hinge-like alternatives.
            direct = abs(i - ai) + abs(j - aj)
            swapped = abs(i - aj) + abs(j - ai)
            grid = min(direct, swapped)
            registry = abs((i - j) - (ai - aj))
            value = float(confidence) * max(exp(-grid / 6.5), 0.55 * exp(-registry / 7.5))
            if value > best:
                best = value
        output[pair] = _rounded(best)
    return output



def _initial_local_geometry_prior(row: RealCoordinateVisualRow, pair: ContactPair, ss: Sequence[str]) -> float:
    i, j = pair
    sep = j - i
    left = ss[i - 1]
    right = ss[j - 1]
    if sep == 3 and left == "H" and right == "H":
        return 0.84
    if sep == 4 and left == "H" and right == "H":
        return 1.0
    if sep == 5 and left == "H" and right == "H":
        return 0.74
    if sep <= 5:
        return 0.45
    if left == "E" and right == "E" and sep >= 8:
        return 0.55
    return 0.0



def _embedded_distances(coords: np.ndarray, candidate_pairs: Sequence[ContactPair]) -> dict[ContactPair, float]:
    output: dict[ContactPair, float] = {}
    for i, j in candidate_pairs:
        delta = coords[i - 1] - coords[j - 1]
        output[(i, j)] = _score(float(np.linalg.norm(delta)))
    return output



def _geometry_scores(distances: Mapping[ContactPair, float]) -> dict[ContactPair, float]:
    return {pair: _rounded(1.0 / (1.0 + exp((float(distance) - CONTACT_CUTOFF_ANGSTROM) / 1.35))) for pair, distance in distances.items()}



def _combined_scores(
    *,
    row: RealCoordinateVisualRow,
    candidate_pairs: Sequence[ContactPair],
    coupling_scores: Mapping[ContactPair, float],
    anchor_field_scores: Mapping[ContactPair, float],
    event_region_scores: Mapping[ContactPair, float],
    physical_scores: Mapping[ContactPair, float],
    local_geometry_scores: Mapping[ContactPair, float],
    geometry_scores: Mapping[ContactPair, float] | None,
    iteration_index: int,
) -> dict[ContactPair, float]:
    raw_layers = {
        "coupling": {pair: float(coupling_scores.get(pair, 0.0)) for pair in candidate_pairs},
        "anchor_field": {pair: float(anchor_field_scores.get(pair, 0.0)) for pair in candidate_pairs},
        "event_region": {pair: float(event_region_scores.get(pair, 0.0)) for pair in candidate_pairs},
        "physical": {pair: float(physical_scores.get(pair, 0.0)) for pair in candidate_pairs},
        "local_geometry": {pair: float(local_geometry_scores.get(pair, 0.0)) for pair in candidate_pairs},
        "embedded_geometry": {pair: float((geometry_scores or {}).get(pair, 0.0)) for pair in candidate_pairs},
    }
    ranks = {name: _rank_percentiles(values) for name, values in raw_layers.items()}
    geometry_weight = min(0.30, 0.10 + 0.045 * max(0, iteration_index - 1))
    weights = {
        "coupling": 0.26,
        "anchor_field": 0.18,
        "event_region": 0.10,
        "physical": 0.20,
        "local_geometry": 0.10,
        "embedded_geometry": geometry_weight,
    }
    non_geom_total = sum(value for name, value in weights.items() if name != "embedded_geometry")
    scale = (1.0 - geometry_weight) / max(non_geom_total, 1e-9)
    weights = {name: (value if name == "embedded_geometry" else value * scale) for name, value in weights.items()}
    return {
        pair: _score(sum(weights[name] * ranks[name][pair] for name in weights))
        for pair in candidate_pairs
    }



def _degree_limited_selection(
    *,
    scores: Mapping[ContactPair, float],
    budget: int,
    max_degree: int,
    required_pairs: Sequence[ContactPair] = (),
) -> tuple[ContactPair, ...]:
    degrees: dict[int, int] = defaultdict(int)
    selected: list[ContactPair] = []
    selected_set: set[ContactPair] = set()
    for pair in required_pairs:
        if pair not in scores:
            continue
        i, j = pair
        if degrees[i] >= max_degree or degrees[j] >= max_degree:
            continue
        selected.append(pair)
        selected_set.add(pair)
        degrees[i] += 1
        degrees[j] += 1
        if len(selected) >= budget:
            return tuple(normalized_contact_pairs(selected))
    for pair, _value in sorted(scores.items(), key=lambda item: (-item[1], item[0][0], item[0][1])):
        if pair in selected_set:
            continue
        i, j = pair
        if degrees[i] >= max_degree or degrees[j] >= max_degree:
            continue
        selected.append(pair)
        selected_set.add(pair)
        degrees[i] += 1
        degrees[j] += 1
        if len(selected) >= budget:
            break
    return tuple(normalized_contact_pairs(selected))



def _distance_matrix_from_contacts(
    *,
    row: RealCoordinateVisualRow,
    contacts: Sequence[ContactPair],
    contact_scores: Mapping[ContactPair, float],
    secondary_structure: Sequence[str],
) -> np.ndarray:
    n = row.sequence_length
    distances = np.full((n, n), np.inf, dtype=float)
    np.fill_diagonal(distances, 0.0)
    # Hard C-alpha chain geometry and local backbone upper bounds. These are
    # generic polymer facts, not target-structure evidence.
    for idx in range(n - 1):
        distances[idx, idx + 1] = distances[idx + 1, idx] = 3.80
    for idx in range(n - 2):
        distances[idx, idx + 2] = distances[idx + 2, idx] = 5.60
    for idx in range(n - 3):
        target = 5.30 if secondary_structure[idx] == "H" and secondary_structure[idx + 3] == "H" else 7.40
        distances[idx, idx + 3] = distances[idx + 3, idx] = target
    for idx in range(n - 4):
        target = 6.10 if secondary_structure[idx] == "H" and secondary_structure[idx + 4] == "H" else 9.40
        distances[idx, idx + 4] = distances[idx + 4, idx] = min(distances[idx, idx + 4], target)

    for pair in contacts:
        i, j = pair
        score_value = max(0.0, min(1.0, float(contact_scores.get(pair, 0.5))))
        target = 5.70 + 2.45 * (1.0 - score_value)
        left = i - 1
        right = j - 1
        if target < distances[left, right]:
            distances[left, right] = distances[right, left] = target

    # Floyd-Warshall all-pairs upper-bound closure. n is <= a few hundred in the
    # locked benchmark, so this bounded deterministic O(n^3) pass is fine.
    for k in range(n):
        through = distances[:, [k]] + distances[[k], :]
        distances = np.minimum(distances, through)
    return distances



def _classical_mds(distance_matrix: np.ndarray, *, dimensions: int = 3) -> np.ndarray:
    n = distance_matrix.shape[0]
    finite = distance_matrix[np.isfinite(distance_matrix)]
    fallback = float(np.max(finite)) if finite.size else 3.8 * n
    d = np.where(np.isfinite(distance_matrix), distance_matrix, fallback)
    d = (d + d.T) / 2.0
    np.fill_diagonal(d, 0.0)
    j = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * j @ (d ** 2) @ j
    eigvals, eigvecs = np.linalg.eigh(gram)
    order = np.argsort(eigvals)[::-1]
    coords = np.zeros((n, dimensions), dtype=float)
    kept = 0
    for eig_index in order:
        value = float(eigvals[eig_index])
        if value <= 1e-9:
            continue
        coords[:, kept] = eigvecs[:, eig_index] * sqrt(value)
        kept += 1
        if kept == dimensions:
            break
    if kept < dimensions:
        # Deterministic tiny helix-like fallback for degenerate matrices.
        for idx in range(n):
            coords[idx, kept:] = 0.01 * np.array([
                np.cos(idx * 1.618 + axis) for axis in range(dimensions - kept)
            ])
    coords -= coords.mean(axis=0, keepdims=True)
    return coords



def _contacts_from_coords(
    coords: np.ndarray,
    *,
    minimum_sequence_separation: int,
    contact_cutoff_angstrom: float,
) -> tuple[ContactPair, ...]:
    n = coords.shape[0]
    pairs: list[ContactPair] = []
    for i in range(1, n + 1):
        left = coords[i - 1]
        for j in range(i + minimum_sequence_separation, n + 1):
            if float(np.linalg.norm(left - coords[j - 1])) <= contact_cutoff_angstrom:
                pairs.append((i, j))
    return normalized_contact_pairs(pairs)



def run_iterative_distance_geometry_diffusion(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    events: Sequence[NucleusClosureEvent] = (),
    iteration_count: int | None = None,
    attraction_budget_fraction: float = 1.72,
    final_budget_fraction: float = 2.22,
    max_degree: int = 9,
    minimum_sequence_separation: int = MIN_SEQUENCE_SEPARATION,
) -> DiffusionPacket:
    candidate_pairs = _all_candidate_pairs(row, minimum_sequence_separation=minimum_sequence_separation)
    coupling_scores = _safe_coupling_scores(row, constraints)
    seed_pairs = tuple(sorted(coupling_scores, key=lambda pair: (-coupling_scores[pair], pair[0], pair[1])))
    event_region_scores = _event_scores(row, events)
    secondary_structure = predict_lightweight_secondary_structure(row.sequence)
    priors = build_sequence_physical_prior_scores(row=row, candidate_pairs=candidate_pairs, current_pairs=())
    physical_scores = {pair: priors[pair].physical_prior_score for pair in candidate_pairs if pair in priors}
    ss_pair_scores = {pair: priors[pair].secondary_structure_score for pair in candidate_pairs if pair in priors}
    local_geometry_scores = {pair: _initial_local_geometry_prior(row, pair, secondary_structure) for pair in candidate_pairs}
    anchor_scores = _anchor_field(candidate_pairs, coupling_scores)

    # Default is deliberately bounded: this is a rescue path, not a long MD job.
    # Callers can opt in to more iterations with --iteration-count.
    iteration_budget = iteration_count if iteration_count is not None else max(3, min(5, int(sqrt(row.sequence_length)) // 3 + 2))
    attraction_budget = max(len(seed_pairs), int(round(row.sequence_length * attraction_budget_fraction)))
    final_budget = max(attraction_budget, int(round(row.sequence_length * final_budget_fraction)))

    coords: np.ndarray | None = None
    embedded_distance_by_pair: dict[ContactPair, float] = {}
    geometry_scores: dict[ContactPair, float] = {}
    selected_pairs: tuple[ContactPair, ...] = tuple(seed_pairs[:attraction_budget])
    iterations: list[DiffusionIteration] = []
    seen_hashes: set[str] = set()

    for iteration_index in range(1, iteration_budget + 1):
        scores = _combined_scores(
            row=row,
            candidate_pairs=candidate_pairs,
            coupling_scores=coupling_scores,
            anchor_field_scores=anchor_scores,
            event_region_scores=event_region_scores,
            physical_scores=physical_scores,
            local_geometry_scores=local_geometry_scores,
            geometry_scores=geometry_scores,
            iteration_index=iteration_index,
        )
        required = seed_pairs[: min(len(seed_pairs), max(8, len(seed_pairs) // 2))]
        attraction_pairs = _degree_limited_selection(
            scores=scores,
            budget=attraction_budget,
            max_degree=max_degree,
            required_pairs=required,
        )
        distances = _distance_matrix_from_contacts(
            row=row,
            contacts=attraction_pairs,
            contact_scores=scores,
            secondary_structure=secondary_structure,
        )
        coords = _classical_mds(distances)
        embedded_contact_pairs = _contacts_from_coords(
            coords,
            minimum_sequence_separation=minimum_sequence_separation,
            contact_cutoff_angstrom=CONTACT_CUTOFF_ANGSTROM,
        )
        embedded_distance_by_pair = _embedded_distances(coords, candidate_pairs)
        geometry_scores = _geometry_scores(embedded_distance_by_pair)
        # Geometry adds evidence, but the final selection still goes through
        # sequence/DCA priors and degree limits so that an over-compressed MDS
        # solution cannot flood the contact map.
        final_scores = _combined_scores(
            row=row,
            candidate_pairs=candidate_pairs,
            coupling_scores=coupling_scores,
            anchor_field_scores=anchor_scores,
            event_region_scores=event_region_scores,
            physical_scores=physical_scores,
            local_geometry_scores=local_geometry_scores,
            geometry_scores=geometry_scores,
            iteration_index=iteration_index + 1,
        )
        embedded_required = tuple(
            pair for pair in sorted(embedded_contact_pairs, key=lambda pair: (-final_scores.get(pair, 0.0), pair[0], pair[1]))[: max(0, final_budget // 5)]
        )
        selected_pairs = _degree_limited_selection(
            scores=final_scores,
            budget=final_budget,
            max_degree=max_degree,
            required_pairs=tuple(required) + embedded_required,
        )
        selected_hash = contact_map_hash(selected_pairs)
        selected_distances = [embedded_distance_by_pair.get(pair, 0.0) for pair in selected_pairs]
        selected_geometry = [geometry_scores.get(pair, 0.0) for pair in selected_pairs]
        selected_physical = [physical_scores.get(pair, 0.0) for pair in selected_pairs]
        selected_scores = [final_scores.get(pair, 0.0) for pair in selected_pairs]
        iterations.append(
            DiffusionIteration(
                iteration_index=iteration_index,
                candidate_pair_count=len(candidate_pairs),
                seed_anchor_pair_count=len(seed_pairs),
                attraction_pair_count=len(attraction_pairs),
                embedded_contact_pair_count=len(embedded_contact_pairs),
                selected_pair_count=len(selected_pairs),
                selected_contact_map_hash=selected_hash,
                mean_selected_score=_mean(selected_scores),
                mean_selected_embedded_distance_angstrom=_mean(selected_distances),
                mean_selected_physical_prior_score=_mean(selected_physical),
                mean_selected_geometry_score=_mean(selected_geometry),
            )
        )
        if selected_hash in seen_hashes:
            break
        seen_hashes.add(selected_hash)

    if coords is None:
        # No iteration happened only if a caller requested zero iterations. Still
        # return a deterministic native-free packet.
        scores = _combined_scores(
            row=row,
            candidate_pairs=candidate_pairs,
            coupling_scores=coupling_scores,
            anchor_field_scores=anchor_scores,
            event_region_scores=event_region_scores,
            physical_scores=physical_scores,
            local_geometry_scores=local_geometry_scores,
            geometry_scores=None,
            iteration_index=1,
        )
        selected_pairs = _degree_limited_selection(scores=scores, budget=final_budget, max_degree=max_degree, required_pairs=seed_pairs)
    else:
        scores = _combined_scores(
            row=row,
            candidate_pairs=candidate_pairs,
            coupling_scores=coupling_scores,
            anchor_field_scores=anchor_scores,
            event_region_scores=event_region_scores,
            physical_scores=physical_scores,
            local_geometry_scores=local_geometry_scores,
            geometry_scores=geometry_scores,
            iteration_index=iteration_budget + 1,
        )

    selected_pairs = normalized_contact_pairs(selected_pairs)
    selected_set = set(selected_pairs)
    prior_summary = summarize_physical_priors(priors, selected_pairs)
    score_rows: list[DiffusionContactScore] = []
    for pair, value in sorted(scores.items(), key=lambda item: (-item[1], item[0][0], item[0][1])):
        if pair in selected_set or len(score_rows) < max(500, final_budget + 64):
            score_rows.append(
                DiffusionContactScore(
                    row_id=row.row_id,
                    source_accession=row.source_accession,
                    i=pair[0],
                    j=pair[1],
                    sequence_separation=pair[1] - pair[0],
                    final_score=_score(value),
                    coupling_score=_rounded(coupling_scores.get(pair, 0.0)),
                    anchor_field_score=_rounded(anchor_scores.get(pair, 0.0)),
                    event_region_score=_rounded(event_region_scores.get(pair, 0.0)),
                    physical_prior_score=_rounded(physical_scores.get(pair, 0.0)),
                    secondary_structure_pair_score=_rounded(ss_pair_scores.get(pair, 0.0)),
                    geometry_score=_rounded(geometry_scores.get(pair, 0.0)),
                    embedded_distance_angstrom=(
                        _score(embedded_distance_by_pair[pair]) if pair in embedded_distance_by_pair else None
                    ),
                    selected=pair in selected_set,
                )
            )
    return DiffusionPacket(
        kind=ITERATIVE_DISTANCE_GEOMETRY_DIFFUSION_KIND,
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_length=row.sequence_length,
        candidate_pair_count=len(candidate_pairs),
        seed_anchor_pair_count=len(seed_pairs),
        final_pair_count=len(selected_pairs),
        final_long_range_pair_count=sum(1 for pair in selected_pairs if pair[1] - pair[0] >= 24),
        final_contact_map_hash=contact_map_hash(selected_pairs),
        iteration_count=len(iterations),
        decision_rule=ITERATIVE_DISTANCE_GEOMETRY_DIFFUSION_RULE,
        selected_pairs=selected_pairs,
        iterations=tuple(iterations),
        scores=tuple(score_rows),
        mean_final_contact_energy_score=prior_summary["mean_contact_energy_score"],
        mean_final_secondary_structure_score=prior_summary["mean_secondary_structure_score"],
        mean_final_degree_consistency_score=prior_summary["mean_degree_consistency_score"],
        mean_final_physical_prior_score=prior_summary["mean_physical_prior_score"],
    )
