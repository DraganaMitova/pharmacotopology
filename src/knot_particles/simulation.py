from __future__ import annotations

from dataclasses import asdict, dataclass
from math import atan2, cos, isfinite, pi, sin
from typing import Mapping, Optional, Sequence, Tuple


Position = Tuple[int, int]

RING_OFFSETS: tuple[Position, ...] = (
    (0, -1),
    (1, -1),
    (1, 0),
    (1, 1),
    (0, 1),
    (-1, 1),
    (-1, 0),
    (-1, -1),
)

DEFAULT_LINK_OFFSETS: tuple[Position, ...] = RING_OFFSETS


@dataclass(frozen=True)
class KnotFieldNode:
    state: float
    phase: float
    charge: int
    memory: float


@dataclass(frozen=True)
class KnotSeed:
    center: Position
    charge: int = 1
    winding: int = 1
    phase_offset: float = 0.0
    label: str = "seed_loop"


@dataclass(frozen=True)
class KnotParticleConfig:
    width: int = 48
    height: int = 32
    steps: int = 48
    activation_threshold: float = 0.32
    loop_threshold: float = 0.74
    decay_rate: float = 0.24
    spread_gain: float = 0.03
    resonance_gain: float = 0.02
    difference_decay: float = 0.22
    memory_gain: float = 0.18
    memory_decay: float = 0.06
    loop_gain: float = 0.36
    attention_gain: float = 0.96
    attention_drain: float = 0.78
    phase_step: float = pi / 5
    attention_vector: Position = (1, 0)
    decay_emission_mode: str = "periodic"
    decay_emission_gain: float = 0.0
    decay_emission_drain: float = 0.0
    decay_emission_start: int = 24
    decay_emission_period: int = 24
    decay_emission_distance: int = 5
    decay_emission_threshold: float = 0.92
    decay_emission_resource_cost: float = 0.0
    stress_emission_threshold: float = 1.0
    stress_emission_accumulation_gain: float = 0.08
    stress_emission_decay: float = 0.02
    stress_emission_relief: float = 0.85
    stress_emission_crowding_gain: float = 0.0
    resource_budget: float = 0.0
    resource_overbranch_drain: float = 0.0
    bound_structure_radius: float = 5.0
    bound_structure_discount: float = 0.0
    confinement_enabled: bool = False
    confinement_min_neighbors: int = 1
    fragment_grace_period: int = 2
    fragment_unbound_support: float = 1.0
    fragment_unbound_drain: float = 0.0
    min_track_persistence: int = 10
    min_motion_distance: float = 6.0
    link_offsets: tuple[Position, ...] = DEFAULT_LINK_OFFSETS


@dataclass(frozen=True)
class KnotLoopObservation:
    time: int
    object_id: int
    center: Position
    charge: int
    mass: float
    memory: float
    phase_winding: float
    spin: float
    stability_score: float
    family: str
    cells: int


@dataclass(frozen=True)
class KnotFrameSummary:
    time: int
    active_nodes: int
    total_mass: float
    total_charge: float
    loop_count: int
    collision_count: int
    annihilation_count: int
    decay_emission_count: int
    fragment_count: int
    bound_fragment_count: int
    unbound_fragment_count: int


@dataclass(frozen=True)
class KnotTrackSummary:
    object_id: int
    charge: int
    family: str
    first_time: int
    last_time: int
    persistence: int
    path_length: float
    displacement: float
    max_stability_score: float
    mean_mass: float
    fragment_like: bool
    bound_frames: int
    unbound_frames: int
    confined: bool
    stable: bool
    moving: bool


@dataclass(frozen=True)
class KnotParticleMetrics:
    benchmark_kind: str
    simulation_only: bool
    toy_universe_only: bool
    known_particle_claim_created: bool
    first_goal_met: bool
    steps: int
    frames_recorded: int
    stable_loop_count: int
    moving_loop_count: int
    loop_family_count: int
    collisions_detected: int
    annihilations_detected: int
    decay_emissions_detected: int
    fragment_track_count: int
    confined_fragment_track_count: int
    max_bound_fragment_count: int
    initial_charge: float
    final_charge: float
    charge_like_conservation_error: float
    max_track_persistence: int


@dataclass(frozen=True)
class KnotParticleRun:
    config: KnotParticleConfig
    seeds: tuple[KnotSeed, ...]
    frames: tuple[KnotFrameSummary, ...]
    observations: tuple[KnotLoopObservation, ...]
    tracks: tuple[KnotTrackSummary, ...]
    metrics: KnotParticleMetrics

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class _LoopCandidate:
    center: Position
    charge: int
    mass: float
    memory: float
    phase_winding: float
    spin: float
    stability_score: float
    family: str
    cells: int


@dataclass
class _TrackState:
    object_id: int
    charge: int
    family: str
    first_time: int
    last_time: int
    centers: list[Position]
    masses: list[float]
    max_stability_score: float
    fragment_like: bool
    bound_frames: int = 0
    unbound_frames: int = 0


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def _rounded(value: float) -> float:
    return round(value, 6)


def _sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def _wrap_position(position: Position, config: KnotParticleConfig) -> Position:
    x, y = position
    return x % config.width, y % config.height


def _add_position(left: Position, right: Position, config: KnotParticleConfig) -> Position:
    return _wrap_position((left[0] + right[0], left[1] + right[1]), config)


def _subtract_position(
    left: Position,
    right: Position,
    config: KnotParticleConfig,
) -> Position:
    return _wrap_position((left[0] - right[0], left[1] - right[1]), config)


def _phase_wrap(value: float) -> float:
    return value % (2.0 * pi)


def _phase_delta(current: float, previous: float) -> float:
    return ((current - previous + pi) % (2.0 * pi)) - pi


def _phase_from_vector(x_value: float, y_value: float, fallback: float) -> float:
    if abs(x_value) + abs(y_value) < 1e-9:
        return _phase_wrap(fallback)
    return _phase_wrap(atan2(y_value, x_value))


def _torus_distance(
    left: Position,
    right: Position,
    config: KnotParticleConfig,
) -> float:
    dx = abs(left[0] - right[0])
    dy = abs(left[1] - right[1])
    dx = min(dx, config.width - dx)
    dy = min(dy, config.height - dy)
    return (dx * dx + dy * dy) ** 0.5


def _path_step_distance(
    left: Position,
    right: Position,
    config: KnotParticleConfig,
) -> float:
    return _torus_distance(left, right, config)


def default_loop_seed(config: Optional[KnotParticleConfig] = None) -> KnotSeed:
    active_config = config or KnotParticleConfig()
    return KnotSeed(
        center=(active_config.width // 4, active_config.height // 2),
        charge=1,
        winding=1,
        phase_offset=0.0,
        label="positive_winding_loop_seed",
    )


def initialize_knot_field(
    seeds: Sequence[KnotSeed],
    config: KnotParticleConfig,
) -> dict[Position, KnotFieldNode]:
    if config.width < 8 or config.height < 8:
        raise ValueError("KNOT particle simulation requires at least an 8x8 field.")

    nodes: dict[Position, KnotFieldNode] = {}
    for seed in seeds:
        if seed.charge not in (-1, 1):
            raise ValueError("KnotSeed.charge must be -1 or 1.")
        if seed.winding == 0:
            raise ValueError("KnotSeed.winding must not be zero.")

        center = _wrap_position(seed.center, config)
        nodes[center] = KnotFieldNode(
            state=0.18,
            phase=_phase_wrap(seed.phase_offset),
            charge=0,
            memory=1.0,
        )

        for index, offset in enumerate(RING_OFFSETS):
            phase = _phase_wrap(
                seed.phase_offset + seed.winding * index * (2.0 * pi / len(RING_OFFSETS))
            )
            nodes[_add_position(center, offset, config)] = KnotFieldNode(
                state=1.0,
                phase=phase,
                charge=seed.charge,
                memory=0.86,
            )

    return nodes


def _candidate_centers(
    nodes: Mapping[Position, KnotFieldNode],
    config: KnotParticleConfig,
) -> tuple[Position, ...]:
    centers: set[Position] = set()
    for position, node in nodes.items():
        if node.state < config.activation_threshold:
            continue
        for offset in RING_OFFSETS:
            centers.add(_subtract_position(position, offset, config))
    return tuple(sorted(centers))


def _loop_candidate_at(
    center: Position,
    nodes: Mapping[Position, KnotFieldNode],
    config: KnotParticleConfig,
) -> Optional[_LoopCandidate]:
    ring_nodes: list[KnotFieldNode] = []
    for offset in RING_OFFSETS:
        ring_nodes.append(
            nodes.get(
                _add_position(center, offset, config),
                KnotFieldNode(state=0.0, phase=0.0, charge=0, memory=0.0),
            )
        )

    active_nodes = [
        node for node in ring_nodes if node.state >= config.activation_threshold
    ]
    cells = len(active_nodes)
    if not active_nodes:
        return None

    mass = sum(node.state for node in ring_nodes)
    memory = sum(node.memory for node in ring_nodes) / len(RING_OFFSETS)
    charge_pressure = sum(node.state * node.charge for node in ring_nodes)
    charge = _sign(charge_pressure)
    if charge == 0:
        return None

    occupancy_score = cells / len(RING_OFFSETS)
    charge_score = abs(charge_pressure) / max(mass, 1e-9)

    winding_total = 0.0
    for previous, current in zip(ring_nodes, ring_nodes[1:] + ring_nodes[:1]):
        winding_total += _phase_delta(current.phase, previous.phase)
    phase_winding = winding_total / (2.0 * pi)
    nearest_winding = round(phase_winding)
    if nearest_winding == 0:
        return None
    winding_score = max(0.0, 1.0 - abs(phase_winding - nearest_winding))

    stability_score = _clamp(
        0.42 * occupancy_score
        + 0.28 * charge_score
        + 0.22 * winding_score
        + 0.08 * memory
    )
    if stability_score < config.loop_threshold:
        return None

    spin = phase_winding * charge
    family = _family_label(charge=charge, winding=nearest_winding)
    return _LoopCandidate(
        center=center,
        charge=charge,
        mass=_rounded(mass),
        memory=_rounded(memory),
        phase_winding=_rounded(phase_winding),
        spin=_rounded(spin),
        stability_score=_rounded(stability_score),
        family=family,
        cells=cells,
    )


def _family_label(*, charge: int, winding: int) -> str:
    charge_label = "positive" if charge > 0 else "negative"
    winding_label = "clockwise" if winding >= 0 else "counterclockwise"
    return f"{charge_label}_{winding_label}_closure"


def detect_loop_candidates(
    nodes: Mapping[Position, KnotFieldNode],
    config: KnotParticleConfig,
) -> tuple[_LoopCandidate, ...]:
    candidates = [
        candidate
        for center in _candidate_centers(nodes, config)
        if (candidate := _loop_candidate_at(center, nodes, config)) is not None
    ]
    candidates.sort(
        key=lambda candidate: (
            candidate.stability_score,
            candidate.mass,
            -candidate.center[0],
            -candidate.center[1],
        ),
        reverse=True,
    )

    selected: list[_LoopCandidate] = []
    for candidate in candidates:
        if any(
            _torus_distance(candidate.center, existing.center, config) <= 1.0
            for existing in selected
        ):
            continue
        selected.append(candidate)
    return tuple(selected)


def _empty_node() -> KnotFieldNode:
    return KnotFieldNode(state=0.0, phase=0.0, charge=0, memory=0.0)


def _add_boost(
    boosts: dict[Position, dict[int, list[float]]],
    position: Position,
    charge: int,
    state: float,
    phase: float,
    memory: float,
) -> None:
    if charge not in (-1, 1):
        return
    charge_boost = boosts.setdefault(position, {}).setdefault(
        charge,
        [0.0, 0.0, 0.0, 0.0],
    )
    charge_boost[0] += state
    charge_boost[1] += state * cos(phase)
    charge_boost[2] += state * sin(phase)
    charge_boost[3] += memory


def _unit_direction(vector: Position) -> Position:
    x_value, y_value = vector
    if x_value == 0 and y_value == 0:
        return (1, 0)
    if abs(x_value) >= abs(y_value):
        return (_sign(float(x_value)), 0)
    return (0, _sign(float(y_value)))


def _emission_offsets(config: KnotParticleConfig) -> tuple[Position, Position]:
    direction = _unit_direction(config.attention_vector)
    perpendicular = (-direction[1], direction[0])
    distance = max(2, config.decay_emission_distance)
    forward = (direction[0] * distance, direction[1] * distance)
    lateral = (perpendicular[0] * distance, perpendicular[1] * distance)
    return (
        (forward[0] + lateral[0], forward[1] + lateral[1]),
        (forward[0] - lateral[0], forward[1] - lateral[1]),
    )


def _resource_support(
    closures: Sequence[_LoopCandidate],
    config: KnotParticleConfig,
) -> dict[Position, float]:
    if config.resource_budget <= 0.0 or not closures:
        return {closure.center: 1.0 for closure in closures}

    discount = _clamp(config.bound_structure_discount)
    weights: dict[Position, float] = {}
    for closure in closures:
        bound_neighbor = any(
            other.center != closure.center
            and other.charge == -closure.charge
            and _torus_distance(closure.center, other.center, config)
            <= config.bound_structure_radius
            for other in closures
        )
        weights[closure.center] = 1.0 - discount if bound_neighbor else 1.0

    demand = sum(weights.values())
    if demand <= config.resource_budget:
        return {closure.center: 1.0 for closure in closures}

    base_support = _clamp(config.resource_budget / max(demand, 1e-9))
    return {
        center: _clamp(base_support / max(weight, 1e-9))
        for center, weight in weights.items()
    }


def _bound_observation_ids(
    observations: Sequence[KnotLoopObservation],
    config: KnotParticleConfig,
) -> frozenset[int]:
    bound_ids: set[int] = set()
    for observation in observations:
        neighbor_count = sum(
            1
            for other in observations
            if other.object_id != observation.object_id
            and other.charge == -observation.charge
            and _torus_distance(observation.center, other.center, config)
            <= config.bound_structure_radius
        )
        if neighbor_count >= config.confinement_min_neighbors:
            bound_ids.add(observation.object_id)
    return frozenset(bound_ids)


def _update_fragment_binding(
    observations: Sequence[KnotLoopObservation],
    tracks: Mapping[int, _TrackState],
    config: KnotParticleConfig,
) -> tuple[int, int, int]:
    if not config.confinement_enabled:
        return (0, 0, 0)

    bound_ids = _bound_observation_ids(observations, config)
    fragment_count = 0
    bound_count = 0
    unbound_count = 0
    for observation in observations:
        track = tracks.get(observation.object_id)
        if track is None or not track.fragment_like:
            continue
        fragment_count += 1
        if observation.object_id in bound_ids:
            track.bound_frames += 1
            bound_count += 1
        else:
            track.unbound_frames += 1
            unbound_count += 1
    return fragment_count, bound_count, unbound_count


def _fragment_confinement_forces(
    observations: Sequence[KnotLoopObservation],
    tracks: Mapping[int, _TrackState],
    config: KnotParticleConfig,
) -> tuple[dict[Position, float], dict[Position, float]]:
    support: dict[Position, float] = {}
    drains: dict[Position, float] = {}
    if not config.confinement_enabled:
        return support, drains

    bound_ids = _bound_observation_ids(observations, config)
    for observation in observations:
        track = tracks.get(observation.object_id)
        if track is None or not track.fragment_like:
            continue
        age = observation.time - track.first_time
        if observation.object_id in bound_ids:
            support[observation.center] = 1.0
            continue
        if age <= config.fragment_grace_period:
            support[observation.center] = 1.0
            continue
        support[observation.center] = _clamp(config.fragment_unbound_support)
        drains[observation.center] = config.fragment_unbound_drain
    return support, drains


def _should_emit_decay(config: KnotParticleConfig, time: int) -> bool:
    if config.decay_emission_gain <= 0.0:
        return False
    if config.decay_emission_mode != "periodic":
        return False
    if config.decay_emission_period <= 0:
        return False
    if time < config.decay_emission_start:
        return False
    return (time - config.decay_emission_start) % config.decay_emission_period == 0


def _emit_decay_closures(
    closures: Sequence[_LoopCandidate],
    nodes: Mapping[Position, KnotFieldNode],
    config: KnotParticleConfig,
    time: int,
    boosts: dict[Position, dict[int, list[float]]],
    drains: dict[Position, float],
    resource_support: Mapping[Position, float],
    emission_centers: Optional[frozenset[Position]] = None,
) -> int:
    if emission_centers is None and not _should_emit_decay(config, time):
        return 0

    emission_count = 0
    daughter_offsets = _emission_offsets(config)
    for closure in closures:
        if emission_centers is not None and closure.center not in emission_centers:
            continue
        if closure.stability_score < config.decay_emission_threshold:
            continue
        support_scale = resource_support.get(closure.center, 1.0)
        winding = round(closure.phase_winding) or 1
        daughter_specs = (
            (daughter_offsets[0], closure.charge, winding),
            (daughter_offsets[1], -closure.charge, -winding),
        )
        for center_offset, daughter_charge, daughter_winding in daughter_specs:
            daughter_center = _add_position(closure.center, center_offset, config)
            phase_offset = time * config.phase_step
            for index, ring_offset in enumerate(RING_OFFSETS):
                target_position = _add_position(daughter_center, ring_offset, config)
                phase = _phase_wrap(
                    phase_offset
                    + daughter_winding * index * (2.0 * pi / len(RING_OFFSETS))
                )
                _add_boost(
                    boosts,
                    target_position,
                    daughter_charge,
                    config.decay_emission_gain
                    * closure.stability_score
                    * support_scale,
                    phase,
                    config.memory_gain * 0.5,
                )
            emission_count += 1

        for ring_offset in RING_OFFSETS:
            source_position = _add_position(closure.center, ring_offset, config)
            source_node = nodes.get(source_position, _empty_node())
            drains[source_position] = drains.get(source_position, 0.0) + (
                (
                    config.decay_emission_drain
                    + config.decay_emission_resource_cost
                    + config.resource_overbranch_drain * (1.0 - support_scale)
                )
                * closure.stability_score
                * max(source_node.state, config.activation_threshold)
            )

    return emission_count


def _closure_forces(
    closures: Sequence[_LoopCandidate],
    nodes: Mapping[Position, KnotFieldNode],
    config: KnotParticleConfig,
    time: int,
    emission_centers: Optional[frozenset[Position]] = None,
    confinement_support: Optional[Mapping[Position, float]] = None,
    confinement_drains: Optional[Mapping[Position, float]] = None,
) -> tuple[dict[Position, dict[int, list[float]]], dict[Position, float], int]:
    boosts: dict[Position, dict[int, list[float]]] = {}
    drains: dict[Position, float] = {}
    moving_attention = config.attention_vector != (0, 0)
    support_by_center = _resource_support(closures, config)
    confinement_support = confinement_support or {}
    confinement_drains = confinement_drains or {}

    for closure in closures:
        support_scale = (
            support_by_center.get(closure.center, 1.0)
            * confinement_support.get(closure.center, 1.0)
        )
        target_center = (
            _add_position(closure.center, config.attention_vector, config)
            if moving_attention
            else closure.center
        )
        for offset in RING_OFFSETS:
            source_position = _add_position(closure.center, offset, config)
            target_position = _add_position(target_center, offset, config)
            source_node = nodes.get(source_position, _empty_node())
            projected_phase = _phase_wrap(source_node.phase + config.phase_step)
            drains[source_position] = drains.get(source_position, 0.0) + (
                confinement_drains.get(closure.center, 0.0)
                * closure.stability_score
            )
            if moving_attention:
                drains[source_position] = drains.get(source_position, 0.0) + (
                    config.attention_drain * closure.stability_score
                    + config.resource_overbranch_drain
                    * (1.0 - support_scale)
                    * closure.stability_score
                )
                _add_boost(
                    boosts,
                    target_position,
                    closure.charge,
                    config.attention_gain
                    * closure.stability_score
                    * support_scale,
                    projected_phase,
                    config.memory_gain,
                )
            else:
                drains[source_position] = drains.get(source_position, 0.0) + (
                    config.resource_overbranch_drain
                    * (1.0 - support_scale)
                    * closure.stability_score
                )
                _add_boost(
                    boosts,
                    source_position,
                    closure.charge,
                    config.loop_gain * closure.stability_score * support_scale,
                    projected_phase,
                    config.memory_gain,
                )

    emission_count = _emit_decay_closures(
        closures,
        nodes,
        config,
        time,
        boosts,
        drains,
        support_by_center,
        emission_centers,
    )
    return boosts, drains, emission_count


def advance_knot_field(
    nodes: Mapping[Position, KnotFieldNode],
    config: KnotParticleConfig,
    time: int = 0,
    emission_centers: Optional[frozenset[Position]] = None,
    confinement_support: Optional[Mapping[Position, float]] = None,
    confinement_drains: Optional[Mapping[Position, float]] = None,
) -> tuple[dict[Position, KnotFieldNode], dict[str, int]]:
    closures = detect_loop_candidates(nodes, config)
    boosts, drains, emission_count = _closure_forces(
        closures,
        nodes,
        config,
        time,
        emission_centers,
        confinement_support,
        confinement_drains,
    )

    candidate_positions = set(nodes)
    candidate_positions.update(boosts)
    for position in tuple(nodes):
        for offset in config.link_offsets:
            candidate_positions.add(_add_position(position, offset, config))

    next_nodes: dict[Position, KnotFieldNode] = {}
    collisions = 0
    annihilations = 0
    for position in sorted(candidate_positions):
        old = nodes.get(position, _empty_node())
        neighbors = [
            nodes.get(_add_position(position, offset, config), _empty_node())
            for offset in config.link_offsets
        ]
        neighbor_mass = sum(node.state for node in neighbors)
        neighbor_charge = sum(node.state * node.charge for node in neighbors)
        neighbor_phase_x = sum(node.state * cos(node.phase) for node in neighbors)
        neighbor_phase_y = sum(node.state * sin(node.phase) for node in neighbors)
        same_phase = sum(
            node.state * max(0.0, cos(_phase_delta(node.phase, old.phase)))
            for node in neighbors
        )
        opposite_pressure = sum(
            node.state
            for node in neighbors
            if old.charge != 0 and node.charge == -old.charge
        )

        charge_boosts = boosts.get(position, {})
        positive = charge_boosts.get(1, [0.0, 0.0, 0.0, 0.0])
        negative = charge_boosts.get(-1, [0.0, 0.0, 0.0, 0.0])
        conflict = min(positive[0], negative[0])
        if conflict > 0.0:
            collisions += 1
        if conflict >= config.activation_threshold:
            annihilations += 1

        boost_state = max(0.0, positive[0] + negative[0] - 2.0 * conflict)
        boost_charge_pressure = positive[0] - negative[0]
        boost_phase_x = positive[1] + negative[1]
        boost_phase_y = positive[2] + negative[2]
        boost_memory = positive[3] + negative[3]

        state = (
            old.state * (1.0 - config.decay_rate)
            + config.spread_gain * (neighbor_mass / max(1, len(neighbors)))
            + config.resonance_gain * (same_phase / max(1, len(neighbors)))
            + boost_state
            - config.difference_decay * opposite_pressure
            - drains.get(position, 0.0)
            - conflict
        )
        state = _clamp(state)
        memory = _clamp(
            old.memory * (1.0 - config.memory_decay)
            + boost_memory
            + config.memory_gain * max(0.0, state - config.activation_threshold)
        )
        charge_pressure = (
            old.state * old.charge
            + 0.35 * neighbor_charge
            + boost_charge_pressure
        )
        charge = _sign(charge_pressure)
        if state < config.activation_threshold * 0.25:
            charge = 0
        if boost_state > 0.0:
            phase = _phase_from_vector(
                boost_phase_x,
                boost_phase_y,
                old.phase + config.phase_step,
            )
        else:
            phase = _phase_from_vector(
                old.state * cos(old.phase) + config.spread_gain * neighbor_phase_x,
                old.state * sin(old.phase) + config.spread_gain * neighbor_phase_y,
                old.phase + config.phase_step,
            )
        if state >= config.activation_threshold * 0.25 or memory > 0.05:
            next_nodes[position] = KnotFieldNode(
                state=_rounded(state),
                phase=_rounded(phase),
                charge=charge,
                memory=_rounded(memory),
            )

    return next_nodes, {
        "collision_count": collisions,
        "annihilation_count": annihilations,
        "decay_emission_count": emission_count,
    }


def _frame_summary(
    time: int,
    nodes: Mapping[Position, KnotFieldNode],
    config: KnotParticleConfig,
    loop_count: int,
    interactions: Mapping[str, int],
) -> KnotFrameSummary:
    active_nodes = sum(
        1 for node in nodes.values() if node.state >= config.activation_threshold
    )
    total_mass = sum(node.state for node in nodes.values())
    total_charge = sum(node.state * node.charge for node in nodes.values())
    return KnotFrameSummary(
        time=time,
        active_nodes=active_nodes,
        total_mass=_rounded(total_mass),
        total_charge=_rounded(total_charge),
        loop_count=loop_count,
        collision_count=int(interactions.get("collision_count", 0)),
        annihilation_count=int(interactions.get("annihilation_count", 0)),
        decay_emission_count=int(interactions.get("decay_emission_count", 0)),
        fragment_count=int(interactions.get("fragment_count", 0)),
        bound_fragment_count=int(interactions.get("bound_fragment_count", 0)),
        unbound_fragment_count=int(interactions.get("unbound_fragment_count", 0)),
    )


def _assign_observations(
    time: int,
    candidates: Sequence[_LoopCandidate],
    tracks: dict[int, _TrackState],
    config: KnotParticleConfig,
    next_object_id: int,
) -> tuple[list[KnotLoopObservation], int]:
    observations: list[KnotLoopObservation] = []
    unmatched_track_ids = set(tracks)

    for candidate in candidates:
        best_id: Optional[int] = None
        best_distance = float("inf")
        for track_id in unmatched_track_ids:
            track = tracks[track_id]
            if track.charge != candidate.charge:
                continue
            distance = _torus_distance(candidate.center, track.centers[-1], config)
            if distance <= 2.25 and distance < best_distance:
                best_id = track_id
                best_distance = distance

        if best_id is None:
            best_id = next_object_id
            next_object_id += 1
            tracks[best_id] = _TrackState(
                object_id=best_id,
                charge=candidate.charge,
                family=candidate.family,
                first_time=time,
                last_time=time,
                centers=[candidate.center],
                masses=[candidate.mass],
                max_stability_score=candidate.stability_score,
                fragment_like=time > 0,
            )
        else:
            unmatched_track_ids.remove(best_id)
            track = tracks[best_id]
            track.last_time = time
            track.centers.append(candidate.center)
            track.masses.append(candidate.mass)
            track.max_stability_score = max(
                track.max_stability_score,
                candidate.stability_score,
            )

        observations.append(
            KnotLoopObservation(
                time=time,
                object_id=best_id,
                center=candidate.center,
                charge=candidate.charge,
                mass=candidate.mass,
                memory=candidate.memory,
                phase_winding=candidate.phase_winding,
                spin=candidate.spin,
                stability_score=candidate.stability_score,
                family=candidate.family,
                cells=candidate.cells,
            )
        )

    return observations, next_object_id


def _stress_increment(
    observation: KnotLoopObservation,
    loop_count: int,
    config: KnotParticleConfig,
) -> float:
    mass_pressure = _clamp(observation.mass / 8.0)
    memory_pressure = _clamp(observation.memory)
    winding_pressure = _clamp(abs(observation.phase_winding))
    crowding_pressure = max(0, loop_count - 1) * config.stress_emission_crowding_gain
    return config.stress_emission_accumulation_gain * (
        0.42 * observation.stability_score
        + 0.24 * mass_pressure
        + 0.18 * memory_pressure
        + 0.12 * winding_pressure
        + 0.04 * crowding_pressure
    )


def _stress_emission_centers(
    observations: Sequence[KnotLoopObservation],
    *,
    loop_count: int,
    config: KnotParticleConfig,
    time: int,
    stress_state: dict[int, float],
) -> frozenset[Position]:
    if config.decay_emission_mode != "stress":
        return frozenset()
    if config.decay_emission_gain <= 0.0:
        return frozenset()

    active_ids = {observation.object_id for observation in observations}
    for object_id in tuple(stress_state):
        if object_id not in active_ids:
            stress_state[object_id] = _rounded(
                stress_state[object_id] * (1.0 - config.stress_emission_decay)
            )

    emission_centers: set[Position] = set()
    for observation in observations:
        stress = stress_state.get(observation.object_id, 0.0)
        stress = stress * (1.0 - config.stress_emission_decay)
        stress += _stress_increment(observation, loop_count, config)
        if (
            time >= config.decay_emission_start
            and stress >= config.stress_emission_threshold
            and observation.stability_score >= config.decay_emission_threshold
        ):
            emission_centers.add(observation.center)
            stress = max(
                0.0,
                stress
                - config.stress_emission_threshold
                * config.stress_emission_relief,
            )
        stress_state[observation.object_id] = _rounded(stress)

    return frozenset(emission_centers)


def _track_summaries(
    tracks: Mapping[int, _TrackState],
    config: KnotParticleConfig,
) -> tuple[KnotTrackSummary, ...]:
    summaries: list[KnotTrackSummary] = []
    for track in tracks.values():
        path_length = 0.0
        for previous, current in zip(track.centers, track.centers[1:]):
            path_length += _path_step_distance(previous, current, config)
        displacement = (
            _torus_distance(track.centers[0], track.centers[-1], config)
            if track.centers
            else 0.0
        )
        persistence = len(track.centers)
        stable = persistence >= config.min_track_persistence
        moving = stable and path_length >= config.min_motion_distance
        confined = (
            track.fragment_like
            and track.bound_frames > 0
            and track.bound_frames >= track.unbound_frames
            and stable
        )
        summaries.append(
            KnotTrackSummary(
                object_id=track.object_id,
                charge=track.charge,
                family=track.family,
                first_time=track.first_time,
                last_time=track.last_time,
                persistence=persistence,
                path_length=_rounded(path_length),
                displacement=_rounded(displacement),
                max_stability_score=_rounded(track.max_stability_score),
                mean_mass=_rounded(sum(track.masses) / max(1, len(track.masses))),
                fragment_like=track.fragment_like,
                bound_frames=track.bound_frames,
                unbound_frames=track.unbound_frames,
                confined=confined,
                stable=stable,
                moving=moving,
            )
        )

    summaries.sort(key=lambda item: (item.object_id, item.first_time))
    return tuple(summaries)


def _metrics(
    *,
    config: KnotParticleConfig,
    frames: Sequence[KnotFrameSummary],
    tracks: Sequence[KnotTrackSummary],
) -> KnotParticleMetrics:
    initial_charge = frames[0].total_charge if frames else 0.0
    final_charge = frames[-1].total_charge if frames else 0.0
    families = {track.family for track in tracks}
    stable_loop_count = sum(1 for track in tracks if track.stable)
    moving_loop_count = sum(1 for track in tracks if track.moving)
    fragment_track_count = sum(1 for track in tracks if track.fragment_like)
    confined_fragment_track_count = sum(1 for track in tracks if track.confined)
    return KnotParticleMetrics(
        benchmark_kind="knot_particle_like_closure_toy_universe",
        simulation_only=True,
        toy_universe_only=True,
        known_particle_claim_created=False,
        first_goal_met=moving_loop_count > 0,
        steps=config.steps,
        frames_recorded=len(frames),
        stable_loop_count=stable_loop_count,
        moving_loop_count=moving_loop_count,
        loop_family_count=len(families),
        collisions_detected=sum(frame.collision_count for frame in frames),
        annihilations_detected=sum(frame.annihilation_count for frame in frames),
        decay_emissions_detected=sum(
            frame.decay_emission_count for frame in frames
        ),
        fragment_track_count=fragment_track_count,
        confined_fragment_track_count=confined_fragment_track_count,
        max_bound_fragment_count=max(
            (frame.bound_fragment_count for frame in frames),
            default=0,
        ),
        initial_charge=_rounded(initial_charge),
        final_charge=_rounded(final_charge),
        charge_like_conservation_error=_rounded(abs(final_charge - initial_charge)),
        max_track_persistence=max(
            (track.persistence for track in tracks),
            default=0,
        ),
    )


def run_knot_particle_simulation(
    config: Optional[KnotParticleConfig] = None,
    seeds: Optional[Sequence[KnotSeed]] = None,
) -> KnotParticleRun:
    active_config = config or KnotParticleConfig()
    active_seeds = tuple(seeds or (default_loop_seed(active_config),))
    nodes = initialize_knot_field(active_seeds, active_config)

    frames: list[KnotFrameSummary] = []
    observations: list[KnotLoopObservation] = []
    tracks: dict[int, _TrackState] = {}
    stress_state: dict[int, float] = {}
    next_object_id = 1
    interactions = {
        "collision_count": 0,
        "annihilation_count": 0,
        "decay_emission_count": 0,
        "fragment_count": 0,
        "bound_fragment_count": 0,
        "unbound_fragment_count": 0,
    }

    for time in range(active_config.steps + 1):
        candidates = detect_loop_candidates(nodes, active_config)
        assigned, next_object_id = _assign_observations(
            time,
            candidates,
            tracks,
            active_config,
            next_object_id,
        )
        observations.extend(assigned)
        fragment_count, bound_fragment_count, unbound_fragment_count = (
            _update_fragment_binding(assigned, tracks, active_config)
        )
        confinement_support, confinement_drains = _fragment_confinement_forces(
            assigned,
            tracks,
            active_config,
        )
        stress_emission_centers: Optional[frozenset[Position]] = None
        if active_config.decay_emission_mode == "stress":
            stress_emission_centers = _stress_emission_centers(
                assigned,
                loop_count=len(candidates),
                config=active_config,
                time=time,
                stress_state=stress_state,
            )
        frames.append(
            _frame_summary(
                time,
                nodes,
                active_config,
                loop_count=len(candidates),
                interactions={
                    **interactions,
                    "fragment_count": fragment_count,
                    "bound_fragment_count": bound_fragment_count,
                    "unbound_fragment_count": unbound_fragment_count,
                },
            )
        )
        if time == active_config.steps:
            break
        nodes, interactions = advance_knot_field(
            nodes,
            active_config,
            time=time,
            emission_centers=stress_emission_centers,
            confinement_support=confinement_support,
            confinement_drains=confinement_drains,
        )

    track_summaries = _track_summaries(tracks, active_config)
    return KnotParticleRun(
        config=active_config,
        seeds=active_seeds,
        frames=tuple(frames),
        observations=tuple(observations),
        tracks=track_summaries,
        metrics=_metrics(
            config=active_config,
            frames=frames,
            tracks=track_summaries,
        ),
    )


def run_seed_to_stable_moving_loop_experiment() -> KnotParticleRun:
    return run_knot_particle_simulation(
        KnotParticleConfig(
            width=48,
            height=32,
            steps=48,
            attention_vector=(1, 0),
            min_track_persistence=10,
            min_motion_distance=6.0,
        )
    )


def validate_config(config: KnotParticleConfig) -> None:
    numeric_values = {
        "activation_threshold": config.activation_threshold,
        "loop_threshold": config.loop_threshold,
        "decay_rate": config.decay_rate,
        "spread_gain": config.spread_gain,
        "resonance_gain": config.resonance_gain,
        "difference_decay": config.difference_decay,
        "memory_gain": config.memory_gain,
        "memory_decay": config.memory_decay,
        "loop_gain": config.loop_gain,
        "attention_gain": config.attention_gain,
        "attention_drain": config.attention_drain,
        "phase_step": config.phase_step,
        "decay_emission_gain": config.decay_emission_gain,
        "decay_emission_drain": config.decay_emission_drain,
        "decay_emission_threshold": config.decay_emission_threshold,
        "decay_emission_resource_cost": config.decay_emission_resource_cost,
        "stress_emission_threshold": config.stress_emission_threshold,
        "stress_emission_accumulation_gain": (
            config.stress_emission_accumulation_gain
        ),
        "stress_emission_decay": config.stress_emission_decay,
        "stress_emission_relief": config.stress_emission_relief,
        "stress_emission_crowding_gain": config.stress_emission_crowding_gain,
        "resource_budget": config.resource_budget,
        "resource_overbranch_drain": config.resource_overbranch_drain,
        "bound_structure_radius": config.bound_structure_radius,
        "bound_structure_discount": config.bound_structure_discount,
        "fragment_unbound_support": config.fragment_unbound_support,
        "fragment_unbound_drain": config.fragment_unbound_drain,
    }
    for name, value in numeric_values.items():
        if not isfinite(value):
            raise ValueError(f"{name} must be finite.")
    if config.width < 8 or config.height < 8:
        raise ValueError("width and height must be at least 8.")
    if config.steps < 1:
        raise ValueError("steps must be positive.")
    if config.decay_emission_mode not in {"periodic", "stress"}:
        raise ValueError("decay_emission_mode must be 'periodic' or 'stress'.")
    if config.decay_emission_distance < 2:
        raise ValueError("decay_emission_distance must be at least 2.")
    if config.resource_budget < 0.0:
        raise ValueError("resource_budget must not be negative.")
    if config.confinement_min_neighbors < 1:
        raise ValueError("confinement_min_neighbors must be at least 1.")
    if config.fragment_grace_period < 0:
        raise ValueError("fragment_grace_period must not be negative.")


def run_validated_knot_particle_simulation(
    config: Optional[KnotParticleConfig] = None,
    seeds: Optional[Sequence[KnotSeed]] = None,
) -> KnotParticleRun:
    active_config = config or KnotParticleConfig()
    validate_config(active_config)
    return run_knot_particle_simulation(active_config, seeds)
