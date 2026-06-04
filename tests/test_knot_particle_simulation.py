from pathlib import Path
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from knot_particles.simulation import (
    KnotParticleConfig,
    default_loop_seed,
    detect_loop_candidates,
    initialize_knot_field,
    run_knot_particle_simulation,
    run_seed_to_stable_moving_loop_experiment,
)


def test_knot_particle_seed_becomes_stable_moving_loop() -> None:
    run = run_seed_to_stable_moving_loop_experiment()

    assert run.metrics.benchmark_kind == "knot_particle_like_closure_toy_universe"
    assert run.metrics.simulation_only is True
    assert run.metrics.toy_universe_only is True
    assert run.metrics.known_particle_claim_created is False
    assert run.metrics.first_goal_met is True
    assert run.metrics.stable_loop_count >= 1
    assert run.metrics.moving_loop_count >= 1
    assert run.metrics.max_track_persistence >= 10

    moving_tracks = [track for track in run.tracks if track.moving]
    assert moving_tracks
    assert moving_tracks[0].path_length >= 6.0


def test_knot_particle_node_schema_matches_mark_phase_charge_memory_links() -> None:
    config = KnotParticleConfig(width=24, height=16, steps=4)
    nodes = initialize_knot_field((default_loop_seed(config),), config)
    candidates = detect_loop_candidates(nodes, config)

    assert candidates
    assert config.link_offsets
    for node in nodes.values():
        assert 0.0 <= node.state <= 1.0
        assert 0.0 <= node.phase <= 2.0 * 3.141593
        assert node.charge in (-1, 0, 1)
        assert 0.0 <= node.memory <= 1.0


def test_stationary_attention_stabilizes_without_motion() -> None:
    run = run_knot_particle_simulation(
        KnotParticleConfig(
            width=32,
            height=24,
            steps=16,
            attention_vector=(0, 0),
            min_track_persistence=8,
            min_motion_distance=3.0,
        )
    )

    assert run.metrics.first_goal_met is False
    assert run.metrics.stable_loop_count >= 1
    assert run.metrics.moving_loop_count == 0
    assert max(track.displacement for track in run.tracks) == 0.0


def test_one_seed_decay_emission_creates_daughter_loop_family() -> None:
    run = run_knot_particle_simulation(
        KnotParticleConfig(
            width=72,
            height=40,
            steps=72,
            attention_gain=0.82,
            decay_rate=0.4,
            decay_emission_gain=0.5,
            decay_emission_drain=0.5,
            decay_emission_start=24,
            decay_emission_period=36,
            decay_emission_distance=7,
            decay_emission_threshold=0.94,
            min_track_persistence=8,
            min_motion_distance=4.0,
        )
    )

    assert run.metrics.decay_emissions_detected >= 2
    assert run.metrics.loop_family_count >= 2
    assert run.metrics.stable_loop_count >= 3
    assert run.metrics.known_particle_claim_created is False


def test_stress_emission_with_resource_budget_limits_runaway_branching() -> None:
    run = run_knot_particle_simulation(
        KnotParticleConfig(
            width=56,
            height=36,
            steps=100,
            attention_gain=0.82,
            decay_rate=0.4,
            decay_emission_mode="stress",
            decay_emission_gain=0.5,
            decay_emission_drain=0.35,
            decay_emission_resource_cost=0.35,
            decay_emission_start=0,
            decay_emission_distance=7,
            decay_emission_threshold=0.94,
            stress_emission_threshold=1.55,
            stress_emission_accumulation_gain=0.075,
            stress_emission_decay=0.015,
            stress_emission_relief=0.9,
            stress_emission_crowding_gain=0.2,
            resource_budget=4.0,
            resource_overbranch_drain=0.45,
            bound_structure_radius=6.0,
            bound_structure_discount=0.35,
            min_track_persistence=8,
            min_motion_distance=4.0,
        )
    )

    assert run.metrics.decay_emissions_detected > 0
    assert run.metrics.loop_family_count >= 2
    assert run.frames[-1].loop_count < run.metrics.stable_loop_count
    assert run.frames[-1].loop_count <= 4
    assert run.metrics.known_particle_claim_created is False


def test_confinement_prunes_isolated_fragments_but_preserves_bound_fragments() -> None:
    run = run_knot_particle_simulation(
        KnotParticleConfig(
            width=64,
            height=44,
            steps=100,
            attention_gain=0.82,
            attention_drain=0.78,
            decay_rate=0.4,
            decay_emission_mode="stress",
            decay_emission_gain=0.5,
            decay_emission_drain=0.35,
            decay_emission_resource_cost=0.35,
            decay_emission_start=0,
            decay_emission_distance=7,
            decay_emission_threshold=0.94,
            stress_emission_threshold=1.55,
            stress_emission_accumulation_gain=0.075,
            stress_emission_decay=0.015,
            stress_emission_relief=0.9,
            stress_emission_crowding_gain=0.2,
            resource_budget=4.0,
            resource_overbranch_drain=0.45,
            bound_structure_radius=10.0,
            bound_structure_discount=0.35,
            confinement_enabled=True,
            fragment_grace_period=3,
            fragment_unbound_support=0.32,
            fragment_unbound_drain=0.35,
            min_track_persistence=8,
            min_motion_distance=4.0,
        )
    )
    stable_unbound_fragments = [
        track
        for track in run.tracks
        if track.fragment_like and track.stable and not track.confined
    ]

    assert run.metrics.fragment_track_count > 0
    assert run.metrics.confined_fragment_track_count > 0
    assert run.metrics.max_bound_fragment_count > 0
    assert stable_unbound_fragments == []
    assert run.metrics.known_particle_claim_created is False


def test_knot_particle_runner_writes_bounded_artifacts(tmp_path: Path) -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "run_knot_particle_simulation.py"),
            "--run-dir",
            str(tmp_path),
            "--width",
            "32",
            "--height",
            "24",
            "--steps",
            "16",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(completed.stdout)
    report = json.loads(
        (tmp_path / "knot_particle_toy_universe_report.json").read_text(
            encoding="utf-8"
        )
    )

    assert payload["known_particle_claim_created"] is False
    assert payload["report_path"].endswith("knot_particle_toy_universe_report.json")
    assert report["metrics"]["simulation_only"] is True
    assert report["metrics"]["known_particle_claim_created"] is False
    assert (tmp_path / "knot_particle_frames.csv").exists()
    assert (tmp_path / "knot_particle_loop_observations.csv").exists()
    assert (tmp_path / "knot_particle_tracks.csv").exists()
