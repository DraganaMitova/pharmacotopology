from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.artifact_io import write_csv_rows  # noqa: E402
from pharmacotopology.knot_particle_simulation import (  # noqa: E402
    KnotParticleConfig,
    run_validated_knot_particle_simulation,
)


DEFAULT_RUN_DIR = Path(
    "first_contact_clean_pharmacotopology_layer_run/knot_particle_toy_universe"
)


def _center_text(center: tuple[int, int]) -> str:
    return f"{center[0]},{center[1]}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the KNOT particle-like closure toy-universe simulation. "
            "This is a bounded emergence experiment, not a real-particle claim."
        )
    )
    parser.add_argument(
        "--run-dir",
        default=str(DEFAULT_RUN_DIR),
        help="Directory for KNOT particle toy-universe artifacts.",
    )
    parser.add_argument("--width", type=int, default=48)
    parser.add_argument("--height", type=int, default=32)
    parser.add_argument("--steps", type=int, default=48)
    parser.add_argument(
        "--attention-dx",
        type=int,
        default=1,
        help="Horizontal attention/selection pressure for coherent closures.",
    )
    parser.add_argument(
        "--attention-dy",
        type=int,
        default=0,
        help="Vertical attention/selection pressure for coherent closures.",
    )
    parser.add_argument("--attention-gain", type=float, default=0.96)
    parser.add_argument("--attention-drain", type=float, default=0.78)
    parser.add_argument("--decay-rate", type=float, default=0.24)
    parser.add_argument(
        "--decay-emission-gain",
        type=float,
        default=0.0,
        help="Opt-in closure-shedding gain. Keep 0.0 for no daughter emission.",
    )
    parser.add_argument(
        "--decay-emission-mode",
        choices=("periodic", "stress"),
        default="periodic",
        help="Use periodic clocked emission or tracked internal stress emission.",
    )
    parser.add_argument("--decay-emission-drain", type=float, default=0.0)
    parser.add_argument("--decay-emission-start", type=int, default=24)
    parser.add_argument("--decay-emission-period", type=int, default=24)
    parser.add_argument("--decay-emission-distance", type=int, default=5)
    parser.add_argument("--decay-emission-threshold", type=float, default=0.92)
    parser.add_argument("--decay-emission-resource-cost", type=float, default=0.0)
    parser.add_argument("--stress-emission-threshold", type=float, default=1.0)
    parser.add_argument(
        "--stress-emission-accumulation-gain",
        type=float,
        default=0.08,
    )
    parser.add_argument("--stress-emission-decay", type=float, default=0.02)
    parser.add_argument("--stress-emission-relief", type=float, default=0.85)
    parser.add_argument("--stress-emission-crowding-gain", type=float, default=0.0)
    parser.add_argument("--resource-budget", type=float, default=0.0)
    parser.add_argument("--resource-overbranch-drain", type=float, default=0.0)
    parser.add_argument("--bound-structure-radius", type=float, default=5.0)
    parser.add_argument("--bound-structure-discount", type=float, default=0.0)
    parser.add_argument(
        "--confinement-enabled",
        action="store_true",
        help="Penalize emitted fragment-like tracks unless they are locally bound.",
    )
    parser.add_argument("--confinement-min-neighbors", type=int, default=1)
    parser.add_argument("--fragment-grace-period", type=int, default=2)
    parser.add_argument("--fragment-unbound-support", type=float, default=1.0)
    parser.add_argument("--fragment-unbound-drain", type=float, default=0.0)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    config = KnotParticleConfig(
        width=args.width,
        height=args.height,
        steps=args.steps,
        attention_vector=(args.attention_dx, args.attention_dy),
        attention_gain=args.attention_gain,
        attention_drain=args.attention_drain,
        decay_rate=args.decay_rate,
        decay_emission_mode=args.decay_emission_mode,
        decay_emission_gain=args.decay_emission_gain,
        decay_emission_drain=args.decay_emission_drain,
        decay_emission_start=args.decay_emission_start,
        decay_emission_period=args.decay_emission_period,
        decay_emission_distance=args.decay_emission_distance,
        decay_emission_threshold=args.decay_emission_threshold,
        decay_emission_resource_cost=args.decay_emission_resource_cost,
        stress_emission_threshold=args.stress_emission_threshold,
        stress_emission_accumulation_gain=(
            args.stress_emission_accumulation_gain
        ),
        stress_emission_decay=args.stress_emission_decay,
        stress_emission_relief=args.stress_emission_relief,
        stress_emission_crowding_gain=args.stress_emission_crowding_gain,
        resource_budget=args.resource_budget,
        resource_overbranch_drain=args.resource_overbranch_drain,
        bound_structure_radius=args.bound_structure_radius,
        bound_structure_discount=args.bound_structure_discount,
        confinement_enabled=args.confinement_enabled,
        confinement_min_neighbors=args.confinement_min_neighbors,
        fragment_grace_period=args.fragment_grace_period,
        fragment_unbound_support=args.fragment_unbound_support,
        fragment_unbound_drain=args.fragment_unbound_drain,
    )
    run = run_validated_knot_particle_simulation(config)

    report_path = run_dir / "knot_particle_toy_universe_report.json"
    frames_path = run_dir / "knot_particle_frames.csv"
    observations_path = run_dir / "knot_particle_loop_observations.csv"
    tracks_path = run_dir / "knot_particle_tracks.csv"

    report_path.write_text(
        json.dumps(run.to_dict(), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_csv_rows([asdict(frame) for frame in run.frames], frames_path)
    write_csv_rows(
        [
            {
                **asdict(observation),
                "center": _center_text(observation.center),
            }
            for observation in run.observations
        ],
        observations_path,
    )
    write_csv_rows([asdict(track) for track in run.tracks], tracks_path)

    print(
        json.dumps(
            {
                "report_path": str(report_path),
                "frames_path": str(frames_path),
                "observations_path": str(observations_path),
                "tracks_path": str(tracks_path),
                "first_goal_met": run.metrics.first_goal_met,
                "stable_loop_count": run.metrics.stable_loop_count,
                "moving_loop_count": run.metrics.moving_loop_count,
                "loop_family_count": run.metrics.loop_family_count,
                "decay_emissions_detected": run.metrics.decay_emissions_detected,
                "fragment_track_count": run.metrics.fragment_track_count,
                "confined_fragment_track_count": (
                    run.metrics.confined_fragment_track_count
                ),
                "max_bound_fragment_count": run.metrics.max_bound_fragment_count,
                "charge_like_conservation_error": (
                    run.metrics.charge_like_conservation_error
                ),
                "known_particle_claim_created": (
                    run.metrics.known_particle_claim_created
                ),
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
