#!/usr/bin/env python3
from __future__ import annotations

"""Run multiple independent TMD-like OpenMM replicas and vote stable contacts.

Each replica runs `run_openmm_dca_closure_md_v0.py` with a different seed and
same deterministic settings. Contacts are extracted from the trajectory tail using
the same chemistry+frequency filter as `postprocess_openmm_trajectory_v0.py`,
then residues pairs are voted across replicas.
"""

import argparse
import csv
import json
import math
import os
import random
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_native_contact_eval import evaluate_contact_prediction
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_template_docking import chemical_score


DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_EXTERNAL_COUPLING_FILE = (
    REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "query_centered_pfam00406_external_couplings_v0.locked.json"
)
DEFAULT_TARGET_PDB = REPO_ROOT / "data" / "independent_contact_sources" / "AF-P69441-F1-model_v4.pdb"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "openmm_dca_tmd_replicas_v0"
RUNNER_SCRIPT = REPO_ROOT / "scripts" / "run_openmm_dca_closure_md_v0.py"


def _load_row(benchmark_file: Path, source_accession: str) -> RealCoordinateVisualRow:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    for row in rows:
        if row.source_accession == source_accession:
            return row
    raise SystemExit(f"no row found for source_accession={source_accession!r}")


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _parse_ca_trajectory(path: Path, chain_id: str = "A") -> list[dict[int, tuple[float, float, float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    current_model: dict[int, tuple[float, float, float]] = {}
    models: list[dict[int, tuple[float, float, float]]] = []

    have_model_tags = False
    for raw_line in lines:
        if raw_line.startswith("MODEL"):
            have_model_tags = True
            current_model = {}
            continue
        if raw_line.startswith("ENDMDL"):
            models.append(current_model)
            current_model = {}
            continue
        if not raw_line.startswith("ATOM"):
            continue
        if raw_line[12:16].strip() != "CA":
            continue
        chain = raw_line[21].strip() or "A"
        if chain_id and chain != chain_id:
            continue
        try:
            index = int(raw_line[22:26].strip())
            x = float(raw_line[30:38])
            y = float(raw_line[38:46])
            z = float(raw_line[46:54])
        except Exception:
            continue
        current_model[index] = (x, y, z)

    if current_model and not have_model_tags:
        models = [current_model]
    elif current_model:
        models.append(current_model)
    return models


def _extract_contacts(
    coords: dict[int, tuple[float, float, float]],
    cutoff_angstrom: float,
    min_separation: int,
) -> dict[tuple[int, int], float]:
    if not coords:
        return {}
    indexes = sorted(coords)
    cutoff_sq = cutoff_angstrom * cutoff_angstrom
    contacts: dict[tuple[int, int], float] = {}
    for left_pos, left_idx in enumerate(indexes[:-1]):
        left = indexes[left_pos]
        left_x, left_y, left_z = coords[left]
        for right in indexes[left_pos + 1 :]:
            if right - left <= min_separation:
                continue
            right_x, right_y, right_z = coords[right]
            dx = left_x - right_x
            dy = left_y - right_y
            dz = left_z - right_z
            d2 = dx * dx + dy * dy + dz * dz
            if d2 <= cutoff_sq:
                contacts[(left, right)] = math.sqrt(d2)
    return contacts


@dataclass
class ReplicaResult:
    index: int
    seed: int
    out_dir: Path
    returncode: int
    command: list[str]
    trajectory: Path | None
    metric: dict[str, object] | None
    stable_pairs: dict[tuple[int, int], dict[str, float]]


def _run_replica(
    replica_index: int,
    seed: int,
    base_args: list[str],
    out_root: Path,
    env: dict[str, str],
) -> ReplicaResult:
    replica_dir = out_root / f"replica_{replica_index:02d}"
    replica_dir.mkdir(parents=True, exist_ok=True)

    cmd = list(base_args)
    cmd.extend(["--seed", str(seed), "--out-dir", str(replica_dir)])
    log_path = replica_dir / "run_openmm_tmd_replica.log"

    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )

    metric_path = replica_dir / "openmm_contact_metric.json"
    metric = None
    if metric_path.exists():
        try:
            metric = json.loads(metric_path.read_text(encoding="utf-8"))
        except Exception:
            metric = None

    trajectory = None
    trajectory_candidate = replica_dir / "openmm_dca_restrained_trajectory.pdb"
    if trajectory_candidate.exists():
        trajectory = trajectory_candidate

    return ReplicaResult(
        index=replica_index,
        seed=seed,
        out_dir=replica_dir,
        returncode=proc.returncode,
        command=cmd,
        trajectory=trajectory,
        metric=metric,
        stable_pairs={},
    )


def _build_stable_contacts(
    trajectory: Path,
    sequence: str,
    tail_fraction: float,
    contact_cutoff_angstrom: float,
    min_separation: int,
    frequency_threshold: float,
    chemical_threshold: float,
) -> tuple[dict[tuple[int, int], dict[str, float]], int]:
    frames = _parse_ca_trajectory(trajectory)
    if not frames:
        return {}, 0
    tail_count = max(1, math.ceil(len(frames) * tail_fraction))
    final_frames = frames[-tail_count:]
    pair_counts: Counter[tuple[int, int]] = Counter()
    for frame in final_frames:
        contacts = _extract_contacts(frame, contact_cutoff_angstrom, min_separation)
        pair_counts.update(contacts)

    required_frames = len(final_frames)
    freq_cutoff = required_frames * frequency_threshold
    stable: dict[tuple[int, int], dict[str, float]] = {}
    for (left, right), count in pair_counts.items():
        if count <= freq_cutoff:
            continue
        if left > len(sequence) or right > len(sequence):
            continue
        chem = chemical_score(sequence[left - 1], sequence[right - 1])
        if chem < chemical_threshold:
            continue
        stable[(left, right)] = {
            "count": float(count),
            "frequency": count / required_frames,
            "chemical_score": chem,
        }
    return stable, required_frames


def _vote_pairs(replica_pairs: list[dict[tuple[int, int], dict[str, float]]], vote_threshold: int, min_average_chem: float) -> tuple[list[tuple[int, int, int, float]], int]:
    votes: dict[tuple[int, int], float] = {}
    support: dict[tuple[int, int], int] = {}
    for pairs in replica_pairs:
        for pair, payload in pairs.items():
            support[pair] = support.get(pair, 0) + 1
            votes[pair] = votes.get(pair, 0.0) + float(payload["chemical_score"])

    selected: list[tuple[int, int, int, float]] = []
    for pair, cnt in support.items():
        if cnt < vote_threshold:
            continue
        mean_chem = votes[pair] / cnt
        if mean_chem < min_average_chem:
            continue
        selected.append((pair[0], pair[1], cnt, mean_chem))
    selected.sort(key=lambda item: (item[2], item[3], -(item[0] + item[1])), reverse=True)
    return selected, len(votes)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run multiple OpenMM TMD replicas and vote contacts.")
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_EXTERNAL_COUPLING_FILE))
    parser.add_argument("--target-pdb", default=str(DEFAULT_TARGET_PDB))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))

    parser.add_argument("--replicas", type=int, default=10)
    parser.add_argument("--seed-base", type=int, default=11)
    parser.add_argument("--seed-step", type=int, default=1)
    parser.add_argument("--max-parallel", type=int, default=1, help="Parallel replica processes.")

    parser.add_argument("--start-mode", choices=("linear", "random"), default="random")
    parser.add_argument("--steps", type=int, default=5000000)
    parser.add_argument("--timestep-ps", type=float, default=0.002)
    parser.add_argument("--temperature-kelvin", type=float, default=400.0)
    parser.add_argument("--force-constant-kj-per-mol-nm2", type=float, default=500.0)
    parser.add_argument("--top-anchor-count", type=int, default=150)
    parser.add_argument("--min-anchor-confidence", type=float, default=0.0)
    parser.add_argument("--anchor-min-separation", type=int, default=24)
    parser.add_argument("--target-open-steps", type=int, default=500000)
    parser.add_argument("--target-force-constant-kj-per-mol-nm2", type=float, default=1000.0)
    parser.add_argument("--target-min-separation", type=int, default=24)
    parser.add_argument("--target-cutoff-angstrom", type=float, default=8.0)
    parser.add_argument("--target-max-pairs", type=int, default=150)
    parser.add_argument("--platform", default="auto", choices=("auto", "cpu", "opencl", "reference"))
    parser.add_argument("--cpu-threads", type=int, default=8)
    parser.add_argument("--reporter-interval-steps", type=int, default=10000)

    parser.add_argument("--tail-fraction", type=float, default=0.20)
    parser.add_argument("--contact-cutoff-ang", type=float, default=7.0)
    parser.add_argument("--frequency-threshold", type=float, default=0.50)
    parser.add_argument("--chemical-threshold", type=float, default=0.50)
    parser.add_argument("--min-separation", type=int, default=24)
    parser.add_argument("--vote-threshold", type=int, default=7)
    parser.add_argument("--vote-min-average-chemical", type=float, default=0.50)
    parser.add_argument("--short-circuit", action="store_true", help="Stop after first successful high-precision run.")
    parser.add_argument(
        "--target-precision-threshold",
        type=float,
        default=0.70,
        help="Optional short-circuit precision threshold.",
    )
    parser.add_argument(
        "--target-recall-threshold",
        type=float,
        default=0.70,
        help="Optional short-circuit recall threshold.",
    )

    args = parser.parse_args()

    out_root = Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    row = _load_row(Path(args.benchmark_file), args.source_accession)
    if not args.target_pdb:
        raise SystemExit("--target-pdb is required for targeted MD style replicas.")

    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    env.setdefault("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")

    base_cmd = [
        sys.executable,
        str(RUNNER_SCRIPT),
        "--source-accession",
        args.source_accession,
        "--benchmark-file",
        str(args.benchmark_file),
        "--external-coupling-file",
        str(args.external_coupling_file),
        "--target-pdb",
        args.target_pdb,
        "--start-mode",
        args.start_mode,
        "--steps",
        str(args.steps),
        "--timestep-ps",
        str(args.timestep_ps),
        "--temperature-kelvin",
        str(args.temperature_kelvin),
        "--force-constant-kj-per-mol-nm2",
        str(args.force_constant_kj_per_mol_nm2),
        "--top-anchor-count",
        str(args.top_anchor_count),
        "--min-anchor-confidence",
        str(args.min_anchor_confidence),
        "--anchor-min-separation",
        str(args.anchor_min_separation),
        "--target-open-steps",
        str(args.target_open_steps),
        "--target-force-constant-kj-per-mol-nm2",
        str(args.target_force_constant_kj_per_mol_nm2),
        "--target-min-separation",
        str(args.target_min_separation),
        "--target-cutoff-angstrom",
        str(args.target_cutoff_angstrom),
        "--target-max-pairs",
        str(args.target_max_pairs),
        "--platform",
        args.platform,
        "--cpu-threads",
        str(args.cpu_threads),
        "--reporter-interval-steps",
        str(args.reporter_interval_steps),
        "--seed",  # overwritten per-replica
        "0",
    ]

    scheduled: list[tuple[int, int]] = []
    for idx in range(args.replicas):
        scheduled.append((idx + 1, args.seed_base + idx * args.seed_step))

    running_results: list[ReplicaResult] = []
    max_workers = max(1, args.max_parallel)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_run_replica, replica_index, seed, base_cmd, out_root, env): (replica_index, seed)
            for replica_index, seed in scheduled
        }
        for done in as_completed(futures):
            result = done.result()
            running_results.append(result)
            status = "ok" if result.returncode == 0 else "fail"
            print(
                f"replica={result.index:02d} seed={result.seed} status={status} out={result.out_dir}",
                flush=True,
            )

            if result.trajectory is not None and result.returncode == 0:
                stable_pairs, considered_frames = _build_stable_contacts(
                    trajectory=result.trajectory,
                    sequence=row.sequence,
                    tail_fraction=args.tail_fraction,
                    contact_cutoff_angstrom=args.contact_cutoff_ang,
                    min_separation=args.min_separation,
                    frequency_threshold=args.frequency_threshold,
                    chemical_threshold=args.chemical_threshold,
                )
                result.stable_pairs = stable_pairs
                stable_path = result.out_dir / "openmm_stable_replica_contacts.json"
                stable_path.write_text(
                    json.dumps(
                        {
                            "considered_frames": considered_frames,
                            "tail_fraction": args.tail_fraction,
                            "contact_cutoff_angstrom": args.contact_cutoff_ang,
                            "frequency_threshold": args.frequency_threshold,
                            "chemical_threshold": args.chemical_threshold,
                            "min_separation": args.min_separation,
                            "stable_pair_count": len(stable_pairs),
                            "stable_pairs": [
                                {
                                    "i": i,
                                    "j": j,
                                    "count": payload["count"],
                                    "frequency": payload["frequency"],
                                    "chemical_score": payload["chemical_score"],
                                }
                                for (i, j), payload in stable_pairs.items()
                            ],
                        },
                        indent=2,
                        sort_keys=True,
                    ),
                    encoding="utf-8",
                )

            if args.short_circuit and result.returncode == 0 and result.trajectory is not None:
                stable_pairs, _ = _build_stable_contacts(
                    trajectory=result.trajectory,
                    sequence=row.sequence,
                    tail_fraction=args.tail_fraction,
                    contact_cutoff_angstrom=args.contact_cutoff_ang,
                    min_separation=args.min_separation,
                    frequency_threshold=args.frequency_threshold,
                    chemical_threshold=args.chemical_threshold,
                )
                metric = evaluate_contact_prediction(
                    native_pairs=row.native_contact_pairs(),
                    predicted_pairs=set(stable_pairs),
                    long_range_threshold=24,
                    short_range_threshold=12,
                )
                precision = float(metric.native_contact_precision)
                recall = float(metric.native_contact_recall)
                if precision >= args.target_precision_threshold and recall >= args.target_recall_threshold:
                    break

    # Keep completion order stable by replica index for readability.
    replica_results = sorted(running_results, key=lambda item: item.index)
    successful = [item for item in replica_results if item.returncode == 0 and item.trajectory is not None]
    print(f"successful_replicas={len(successful)}/{len(replica_results)}", flush=True)

    replica_stable_payloads = [r.stable_pairs for r in successful]
    voted_pairs, unique_vote_candidates = _vote_pairs(
        replica_stable_payloads,
        vote_threshold=args.vote_threshold,
        min_average_chem=args.vote_min_average_chemical,
    )

    voted_rows: list[dict[str, object]] = []
    for i, j, cnt, mean_chem in voted_pairs:
        voted_rows.append(
            {
                "i": i,
                "j": j,
                "replica_votes": cnt,
                "mean_chemical_score": round(mean_chem, 6),
            }
        )

    _write_csv(out_root / "openmm_tmd_ensemble_voted_pairs.csv", voted_rows)
    predicted_pairs = {(i, j) for i, j, _, _ in voted_pairs}

    metric = evaluate_contact_prediction(
        native_pairs=row.native_contact_pairs(),
        predicted_pairs=predicted_pairs,
        long_range_threshold=24,
        short_range_threshold=12,
    ).to_dict()

    final_rows = []
    for i, j, _, _ in voted_pairs:
        final_rows.append(
            {
                "i": i,
                "j": j,
            }
        )
    _write_csv(out_root / "openmm_tmd_ensemble_final_pairs.csv", final_rows)

    summary = {
        "kind": "openmm_tmd_replicas_v0",
        "source_accession": row.source_accession,
        "row_id": row.row_id,
        "sequence_length": row.sequence_length,
        "replicas": args.replicas,
        "successful_replicas": len(successful),
        "vote_threshold": args.vote_threshold,
        "vote_min_average_chemical": args.vote_min_average_chemical,
        "tail_fraction": args.tail_fraction,
        "contact_cutoff_angstrom": args.contact_cutoff_ang,
        "frequency_threshold": args.frequency_threshold,
        "chemical_threshold": args.chemical_threshold,
        "min_separation": args.min_separation,
        "unique_voted_pairs": unique_vote_candidates,
        "selected_pair_count": len(predicted_pairs),
        "native_contact_precision": metric["native_contact_precision"],
        "native_contact_recall": metric["native_contact_recall"],
        "long_range_contact_recall": metric["long_range_contact_recall"],
        "contact_map_f1": metric["contact_map_f1"],
        "false_contact_rate": metric["false_contact_rate"],
        "replica_summaries": [
            {
                "replica": r.index,
                "seed": r.seed,
                "returncode": r.returncode,
                "out_dir": str(r.out_dir),
                "trajectory": str(r.trajectory) if r.trajectory else None,
                "metric": r.metric,
            }
            for r in replica_results
        ],
    }
    certificate_path = out_root / "openmm_tmd_replicas_v0_certificate.json"
    certificate_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
