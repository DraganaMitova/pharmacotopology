#!/usr/bin/env python3
from __future__ import annotations

"""Run a dependency-aware coarse C-alpha closure attempt from DCA anchors.

Primary mode:
  - uses OpenMM if available with harmonic distance restraints for DCA anchors.

Fallback mode (explicitly enabled):
  - reuses existing dependency-free bounded C-alpha geometry runner seeded with the
    same anchor set, so the command still produces artifacts in environments
    without OpenMM.
"""

import argparse
import csv
import json
import math
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_native_contact_eval import (
    evaluate_contact_prediction,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    load_real_coordinate_visual_rows,
    RealCoordinateVisualRow,
)
from pharmacotopology.folding_coarse_grain_md_geometry import (
    CoarseGrainMDGeometryPacket,
    EXTERNAL_DCA_MD_MODE,
    run_coarse_grain_md_geometry_packet,
)


DEFAULT_BENCHMARK_FILE = (
    REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
)
DEFAULT_EXTERNAL_COUPLING_FILE = (
    REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
    / "query_centered_pfam00406_external_couplings_v0.locked.json"
)
DEFAULT_OUTPUT_DIR = (
    REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "openmm_dca_closure_md_v0"
)

CA_DISTANCE_TARGET_ANGSTROM = 3.8
CA_DISTANCE_METERS = CA_DISTANCE_TARGET_ANGSTROM / 10.0
CONTACT_CUTOFF_ANGSTROM = 8.0


@dataclass(frozen=True)
class AnchorRecord:
    i: int
    j: int
    confidence: float
    distance_angstrom: float
    force_constant_kj_per_mol_nm2: float


def _load_row(benchmark_file: Path, source_accession: str) -> RealCoordinateVisualRow:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    for row in rows:
        if row.source_accession == source_accession:
            return row
    raise SystemExit(f"no row found for source_accession={source_accession!r}")


def _extract_top_anchors(
    row: RealCoordinateVisualRow,
    external_coupling_path: Path,
    top_anchor_count: int,
    min_anchor_confidence: float,
    minimum_sequence_separation: int,
) -> tuple[AnchorRecord, ...]:
    dataset = load_coupling_dataset(external_coupling_path)
    anchors = []
    for constraint in dataset.constraints:
        if constraint.source_accession != row.source_accession:
            continue
        if constraint.row_id != row.row_id:
            continue
        if constraint.coordinate_truth_used_to_build_constraint:
            continue
        if constraint.native_truth_used_before_coupling_selection:
            continue
        if constraint.structure_model_used:
            continue
        if constraint.i < 1 or constraint.j > row.sequence_length:
            continue
        if (constraint.j - constraint.i) < minimum_sequence_separation:
            continue
        confidence = float(constraint.confidence)
        if confidence < min_anchor_confidence:
            continue
        distance = 6.0
        anchors.append(
            AnchorRecord(
                i=constraint.i,
                j=constraint.j,
                confidence=confidence,
                distance_angstrom=distance,
                force_constant_kj_per_mol_nm2=100.0,
            )
        )
    anchors.sort(key=lambda item: (-item.confidence, item.i, item.j))
    return tuple(anchors[:top_anchor_count])


def _write_csv(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _build_linear_positions(length: int) -> list[tuple[float, float, float]]:
    return [
        (i * CA_DISTANCE_METERS, 0.0, 0.0)
        for i in range(max(length, 1))
    ]


def _build_random_coil_positions(length: int, seed: int = 11) -> list[tuple[float, float, float]]:
    random.seed(seed)
    positions = []
    x = y = z = 0.0
    for _ in range(max(length, 1)):
        x += CA_DISTANCE_METERS * (0.7 + 0.6 * random.random())
        y += 0.2 * (random.random() - 0.5)
        z += 0.2 * (random.random() - 0.5)
        positions.append((x, y, z))
    return positions


def _parse_last_ca_coords(pdb_path: Path) -> list[tuple[int, tuple[float, float, float]]]:
    lines = pdb_path.read_text(encoding="utf-8").splitlines()
    model_active = False
    last_coords: list[tuple[int, tuple[float, float, float]]] = []
    model_coords: list[tuple[int, tuple[float, float, float]]] = []
    pending_model = False
    for line in lines:
        if line.startswith("MODEL"):
            pending_model = True
            model_coords = []
            continue
        if line.startswith("ENDMDL"):
            last_coords = model_coords
            pending_model = False
            continue
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        chain = line[21].strip()
        if chain and chain != "A":
            continue
        try:
            index = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except Exception:
            continue
        if pending_model:
            model_coords.append((index, (x, y, z)))
        else:
            last_coords.append((index, (x, y, z)))
    if not last_coords:
        last_coords = model_coords
    return last_coords


def _ca_contacts_from_coords(
    coords: list[tuple[int, tuple[float, float, float]]],
    cutoff_angstrom: float,
) -> tuple[tuple[int, int, float], ...]:
    if not coords:
        return ()
    coord_map = dict(coords)
    ordered_indexes = sorted(coord_map)
    cutoff_sq = cutoff_angstrom * cutoff_angstrom
    out = []
    for left_i, index_left in enumerate(ordered_indexes[:-1], start=0):
        left_x, left_y, left_z = coord_map[index_left]
        for index_right in ordered_indexes[left_i + 1 :]:
            if index_right - index_left <= 2:
                continue
            right_x, right_y, right_z = coord_map[index_right]
            dx = left_x - right_x
            dy = left_y - right_y
            dz = left_z - right_z
            dist_sq = dx * dx + dy * dy + dz * dz
            if dist_sq <= cutoff_sq:
                out.append((index_left, index_right, math.sqrt(dist_sq)))
    return tuple(sorted({(left, right, round(distance, 6)) for left, right, distance in out}))


def _evaluate_and_report(
    row: RealCoordinateVisualRow,
    predicted_contacts: tuple[tuple[int, int, float], ...],
    metric_out_path: Path,
) -> dict:
    predicted_pairs = {(left, right) for left, right, _ in predicted_contacts}
    metric = evaluate_contact_prediction(
        native_pairs=row.native_contact_pairs(),
        predicted_pairs=predicted_pairs,
        long_range_threshold=24,
        short_range_threshold=12,
    )
    payload = metric.to_dict()
    payload["source_accession"] = row.source_accession
    payload["row_id"] = row.row_id
    payload["sequence_length"] = row.sequence_length
    payload["predicted_contacts_from_last_frame"] = len(predicted_pairs)
    metric_out_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def _run_openmm(
    row: RealCoordinateVisualRow,
    anchors: Sequence[AnchorRecord],
    out_dir: Path,
    start_mode: str,
    steps: int,
    timestep: float,
    temperature_kelvin: float,
    seed: int,
    platform_name: str = "auto",
    cpu_threads: int = 12,
    reporter_interval_steps: int = 0,
    write_trajectory: bool = True,
) -> dict:
    try:
        import openmm as mm
        from openmm import app
        from openmm.unit import amu, kelvin, kilojoules_per_mole, nanometer, picosecond
    except ImportError:
        return {
            "success": False,
            "error": "openmm_not_available",
            "message": "OpenMM module is not installed in this environment.",
        }

    topology = app.Topology()
    chain = topology.addChain("A")
    residues = []
    for index, aa in enumerate(row.sequence, start=1):
        residue = topology.addResidue("ALA", chain)
        topology.addAtom("CA", app.Element.getBySymbol("C"), residue)
        residues.append((index, aa))
    if start_mode == "random":
        positions = _build_random_coil_positions(row.sequence_length, seed=seed)
    else:
        positions = _build_linear_positions(row.sequence_length)

    system = mm.System()
    for _ in range(row.sequence_length):
        system.addParticle(12.0 * amu)

    chain_force = mm.HarmonicBondForce()
    for left in range(row.sequence_length - 1):
        chain_force.addBond(
            left,
            left + 1,
            CA_DISTANCE_METERS * nanometer,
            800.0 * kilojoules_per_mole / nanometer ** 2,
        )
    system.addForce(chain_force)

    restraint_force = mm.CustomBondForce(
        "k*step(r-(r0+w))*(r-(r0+w))^2 + k*step((r0-w)-r)*( (r0-w)-r )^2"
    )
    restraint_force.addGlobalParameter(
        "k",
        100.0 * kilojoules_per_mole / (nanometer ** 2),
    )
    restraint_force.addGlobalParameter("w", 0.2 * nanometer)
    restraint_force.addPerBondParameter("r0")
    for anchor in anchors:
        r0 = anchor.distance_angstrom / 10.0 * nanometer
        restraint_force.addBond(
            anchor.i - 1,
            anchor.j - 1,
            [r0],
        )
    system.addForce(restraint_force)

    # masses are already set above via System.addParticle

    integrator = mm.LangevinIntegrator(
        temperature_kelvin * kelvin,
        1.0 / picosecond,
        timestep * picosecond,
    )
    integrator.setRandomNumberSeed(seed)

    def _make_platform(name: str) -> mm.Platform:
        if name == "auto":
            try:
                return _make_platform("opencl")
            except Exception:
                return _make_platform("cpu")
        if name == "opencl":
            return mm.Platform.getPlatformByName("OpenCL")
        if name == "reference":
            return mm.Platform.getPlatformByName("Reference")
        platform = mm.Platform.getPlatformByName("CPU")
        if cpu_threads > 0:
            platform.setPropertyDefaultValue("Threads", str(cpu_threads))
        return platform

    try:
        platform = _make_platform(platform_name)
        simulation = app.Simulation(topology, system, integrator, platform)
    except Exception as exc:
        if platform_name in {"auto", "opencl"}:
            simulation = app.Simulation(
                topology,
                system,
                integrator,
                _make_platform("cpu"),
            )
        else:
            raise SystemExit(f"platform init failed: {exc}")

    simulation.context.setPositions([(x, y, z) * nanometer for x, y, z in positions])
    simulation.minimizeEnergy(maxIterations=500)

    if reporter_interval_steps and reporter_interval_steps > 0:
        reporter_interval = reporter_interval_steps
    else:
        reporter_interval = max(1, min(10000, max(1, steps // 10000)))
    if write_trajectory:
        simulation.reporters.append(app.PDBReporter(str(trajectory_path), reporter_interval))
    simulation.reporters.append(app.StateDataReporter(
        str(out_dir / "openmm_dca_openmm.log"),
        reporter_interval,
        step=True,
        potentialEnergy=True,
        temperature=True,
    ))
    simulation.step(steps)

    state = simulation.context.getState(getPositions=True)
    final_positions = state.getPositions(asNumpy=False)
    final_pdb = out_dir / "openmm_dca_restrained_final.pdb"
    with final_pdb.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(topology, final_positions, handle)

    trajectory_path = out_dir / "openmm_dca_restrained_trajectory.pdb"
    return {
        "success": True,
        "trajectory_pdb": str(trajectory_path) if write_trajectory and trajectory_path.is_file() else None,
        "final_pdb": str(final_pdb),
        "anchors_applied": len(anchors),
        "log": str(out_dir / "openmm_dca_openmm.log"),
    }


def _run_fallback_dependency_free(
    row: RealCoordinateVisualRow,
    anchors: Sequence[AnchorRecord],
    out_dir: Path,
    md_steps: int,
    max_restraints: int,
    coupling_file: Path,
) -> dict:
    # Use existing dependency-free C-alpha geometry runner with the same anchor
    # source to provide deterministic fallback artifacts.
    dataset = load_coupling_dataset(coupling_file)
    row_constraints = tuple(
        constraint
        for constraint in dataset.constraints
        if constraint.source_accession == row.source_accession
        and constraint.row_id == row.row_id
    )
    packet = run_coarse_grain_md_geometry_packet(
        (row,),
        constraints=(
            row_constraints
        ),
        source_mode=EXTERNAL_DCA_MD_MODE,
        md_steps=max(8, min(72, md_steps // 1000)),
        restarts=1,
        max_restraints=max(10, min(max_restraints, len(anchors))),
    )
    row_report = packet.rows[0] if packet.rows else None
    report_path = out_dir / "dependency_free_fallback_report.json"
    report_path.write_text(json.dumps(packet.to_dict(), indent=2, sort_keys=True), encoding="utf-8")
    return {
        "success": True,
        "fallback": True,
        "report": str(report_path),
        "fallback_metric": (
            row_report.metric_after_native_audit.to_dict()
            if row_report is not None
            else None
        ),
        "fallback_row_precision": packet.mean_native_contact_precision_after_audit,
        "fallback_row_recall": packet.mean_native_contact_recall_after_audit,
        "fallback_long_range": packet.mean_long_range_contact_recall_after_audit,
        "fallback_f1": packet.mean_contact_map_f1_after_audit,
        "fallback_row_contacts": row_report.extracted_contact_count
        if row_report is not None
        else 0,
    }


def _build_metric_rows(
    row: RealCoordinateVisualRow,
    anchors: Sequence[AnchorRecord],
    metric: dict,
    openmm_results: dict,
) -> dict:
    return {
        "kind": "openmm_dca_closure_v0",
        "row_id": row.row_id,
        "source_accession": row.source_accession,
        "sequence_length": row.sequence_length,
        "anchor_count": len(anchors),
        "anchor_distance_target_angstrom": 6.0,
        "anchor_distance_window_angstrom": 2.0,
        "anchor_force_constant_kj_mol_nm2": 100.0,
        "mean_precision_after_audit": metric["native_contact_precision"],
        "mean_recall_after_audit": metric["native_contact_recall"],
        "mean_long_range_recall_after_audit": metric["long_range_contact_recall"],
        "mean_f1_after_audit": metric["contact_map_f1"],
        "openmm_success": bool(openmm_results.get("success")) and not bool(
            openmm_results.get("fallback")
        ),
        "fallback_used": bool(openmm_results.get("fallback")),
        "openmm_engine": (
            "dependency_free_fallback"
            if openmm_results.get("fallback")
            else "OpenMM"
        ),
        "trajectory_pdb": openmm_results.get("trajectory_pdb"),
        "final_pdb": openmm_results.get("final_pdb"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run C-alpha MD closure from DCA anchors.")
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--external-coupling-file",
        default=str(DEFAULT_EXTERNAL_COUPLING_FILE),
    )
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--top-anchor-count", type=int, default=50)
    parser.add_argument("--min-anchor-confidence", type=float, default=0.0)
    parser.add_argument("--anchor-min-separation", type=int, default=24)
    parser.add_argument("--start-mode", choices=("linear", "random"), default="linear")
    parser.add_argument("--steps", type=int, default=500000)
    parser.add_argument("--timestep-ps", type=float, default=0.002)
    parser.add_argument("--temperature-kelvin", type=float, default=300.0)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument(
        "--fallback",
        action="store_true",
        help="Run dependency-free fallback if OpenMM is unavailable.",
    )
    parser.add_argument(
        "--platform",
        default="auto",
        choices=("auto", "cpu", "opencl", "reference"),
        help=(
            "OpenMM platform to use. auto tries OpenCL first and falls back to CPU."
        ),
    )
    parser.add_argument(
        "--cpu-threads",
        type=int,
        default=12,
        help="CPU threads for OpenMM CPU platform (for faster CPU execution).",
    )
    parser.add_argument(
        "--reporter-interval-steps",
        type=int,
        default=0,
        help=(
            "Custom reporter interval (steps). Set to 0 for automatic interval; "
            "larger values reduce I/O and can improve throughput."
        ),
    )
    parser.add_argument(
        "--no-trajectory",
        action="store_true",
        help="Skip writing trajectory PDB file for speed.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    benchmark_file = Path(args.benchmark_file)
    coupling_file = Path(args.external_coupling_file)
    row = _load_row(benchmark_file, args.source_accession)
    anchors = _extract_top_anchors(
        row=row,
        external_coupling_path=coupling_file,
        top_anchor_count=args.top_anchor_count,
        min_anchor_confidence=args.min_anchor_confidence,
        minimum_sequence_separation=args.anchor_min_separation,
    )
    if not anchors:
        raise SystemExit("no anchors selected for this run")

    anchors_path = out_dir / "dca_anchors.csv"
    _write_csv(
        anchors_path,
        [
            {
                "i": anchor.i,
                "j": anchor.j,
                "confidence": anchor.confidence,
                "distance_angstrom": anchor.distance_angstrom,
                "force_constant_kj_mol_nm2": anchor.force_constant_kj_per_mol_nm2,
            }
            for anchor in anchors
        ],
    )

    openmm_results = _run_openmm(
        row=row,
        anchors=anchors,
        out_dir=out_dir,
        start_mode=args.start_mode,
        steps=args.steps,
        timestep=args.timestep_ps,
        temperature_kelvin=args.temperature_kelvin,
        seed=args.seed,
        platform_name=args.platform,
        cpu_threads=args.cpu_threads,
        reporter_interval_steps=args.reporter_interval_steps,
        write_trajectory=not args.no_trajectory,
    )
    if not openmm_results.get("success") and args.fallback:
        openmm_results = _run_fallback_dependency_free(
            row=row,
            anchors=anchors,
            out_dir=out_dir,
            md_steps=args.steps,
            max_restraints=args.top_anchor_count,
            coupling_file=coupling_file,
        )

    final_pdb_path = Path(openmm_results["final_pdb"]) if openmm_results.get("final_pdb") else None
    if final_pdb_path is not None and not final_pdb_path.is_file():
        final_pdb_path = None
    if final_pdb_path is not None:
        contacts = _ca_contacts_from_coords(
            _parse_last_ca_coords(final_pdb_path),
            cutoff_angstrom=CONTACT_CUTOFF_ANGSTROM,
        )
        _write_csv(
            out_dir / "openmm_final_frame_contacts.csv",
            [
                {"i": left, "j": right, "distance_angstrom": distance}
                for left, right, distance in contacts
            ],
        )
        predicted_contacts = contacts
    else:
        _write_csv(out_dir / "openmm_final_frame_contacts.csv", [])
        predicted_contacts = ()

    fallback_metric = openmm_results.get("fallback_metric")
    if openmm_results.get("fallback") and isinstance(fallback_metric, dict):
        metric = dict(fallback_metric)
        metric["source_accession"] = row.source_accession
        metric["row_id"] = row.row_id
        metric["sequence_length"] = row.sequence_length
        metric["predicted_contacts_from_last_frame"] = openmm_results.get(
            "fallback_row_contacts",
            0,
        )
        metric_out_path = out_dir / "openmm_contact_metric.json"
        metric_out_path.write_text(json.dumps(metric, indent=2, sort_keys=True), encoding="utf-8")
    else:
        metric = _evaluate_and_report(
            row=row,
            predicted_contacts=predicted_contacts,
            metric_out_path=out_dir / "openmm_contact_metric.json",
        )


    summary = _build_metric_rows(row, anchors, metric, openmm_results)
    summary["openmm_results"] = openmm_results
    certificate_path = out_dir / "openmm_dca_closure_md_certificate_v0.json"
    certificate_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
