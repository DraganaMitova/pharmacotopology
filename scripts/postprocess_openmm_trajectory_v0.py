#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Sequence

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_native_contact_eval import evaluate_contact_prediction
from pharmacotopology.folding_five_axis_physics import matched_control_pairs
from pharmacotopology.folding_template_docking import chemical_score


DEFAULT_EXTERNAL_COUPLING_FILE = (
    Path(__file__).resolve().parents[1]
    / "first_contact_clean_pharmacotopology_layer_run"
    / "query_centered_pfam00406_external_couplings_v0.locked.json"
)


def _load_row(benchmark_file: Path, source_accession: str) -> RealCoordinateVisualRow:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    for row in rows:
        if row.source_accession == source_accession:
            return row
    raise SystemExit(f"no row found for source_accession={source_accession!r}")


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


def _load_dca_scores(
    *,
    coupling_file: Path | None,
    row: RealCoordinateVisualRow,
    sequence_length: int,
) -> dict[tuple[int, int], float]:
    if coupling_file is None:
        return {}
    if not coupling_file.exists():
        return {}
    try:
        dataset = load_coupling_dataset(coupling_file)
    except Exception:
        return {}
    scores: dict[tuple[int, int], float] = {}
    for constraint in dataset.constraints:
        if constraint.row_id != row.row_id and constraint.source_accession != row.source_accession:
            continue
        if not (1 <= constraint.i <= sequence_length and 1 <= constraint.j <= sequence_length):
            continue
        if constraint.i >= constraint.j:
            continue
        pair = (constraint.i, constraint.j)
        prev = scores.get(pair, 0.0)
        confidence = float(constraint.confidence)
        if confidence > prev:
            scores[pair] = confidence
    return scores


def _auto_threshold(values: Sequence[float], *, default: float) -> float:
    sorted_values = sorted(set(float(value) for value in values), reverse=True)
    if len(sorted_values) < 2:
        return default
    gaps: list[tuple[float, int]] = []
    for index, left_value in enumerate(sorted_values[:-1]):
        right_value = sorted_values[index + 1]
        gaps.append((left_value - right_value, index))
    if not gaps:
        return default
    gap, index = max(gaps, key=lambda item: item[0])
    if gap <= 0.0:
        return default
    return max(0.0, sorted_values[index + 1])


def _resolve_threshold(raw: str | float | None, values: Sequence[float], *, default: float) -> float:
    if raw is None:
        return default
    if isinstance(raw, str):
        candidate = raw.strip().lower()
        if candidate == "auto":
            return _auto_threshold(values, default=default)
        return float(candidate)
    return float(raw)


def _parse_domain_boundaries(raw: str | None) -> tuple[tuple[int, int], ...]:
    if not raw:
        return ()
    bounds: list[tuple[int, int]] = []
    for chunk in [item.strip() for item in raw.split(",") if item.strip()]:
        if "-" not in chunk:
            raise ValueError(f"invalid domain boundary segment: {chunk!r}")
        start_text, end_text = chunk.split("-", maxsplit=1)
        start = int(start_text.strip())
        end = int(end_text.strip())
        if start < 1 or end < start:
            raise ValueError(f"invalid domain boundary: {chunk!r}")
        bounds.append((start, end))
    if not bounds:
        return ()
    # Sort and allow overlap; the first matching range wins.
    return tuple(sorted(bounds, key=lambda item: (item[0], item[1])))


def _residue_domain_index(position: int, boundaries: tuple[tuple[int, int], ...]) -> int:
    for index, (start, end) in enumerate(boundaries):
        if start <= position <= end:
            return index
    return -1


def _topology_ok(
    left: int,
    right: int,
    *,
    domain_boundaries: tuple[tuple[int, int], ...],
    topology_mode: str,
) -> bool:
    if topology_mode == "none" or not domain_boundaries:
        return True
    left_domain = _residue_domain_index(left, domain_boundaries)
    right_domain = _residue_domain_index(right, domain_boundaries)
    if left_domain == -1 or right_domain == -1:
        return True
    if topology_mode == "interdomain":
        return left_domain != right_domain
    if topology_mode == "intradomain":
        return left_domain == right_domain
    return True


def _control_gap(
    *,
    row: RealCoordinateVisualRow,
    selected_pairs: Sequence[tuple[int, int]],
    candidate_pairs: Sequence[tuple[int, int]],
    control_count: int,
) -> dict[tuple[int, int], int]:
    if control_count <= 0 or not selected_pairs:
        return {}
    control_hits: dict[tuple[int, int], int] = defaultdict(int)
    for control_index in range(1, control_count + 1):
        control_pairs = matched_control_pairs(
            row=row,
            selected_pairs=selected_pairs,
            candidate_pairs=candidate_pairs,
            control_index=control_index,
        )
        for pair in control_pairs:
            control_hits[pair] += 1
    return control_hits


def _apply_degree_cap(
    pairs: list[tuple[tuple[int, int], dict[str, float]]],
    *,
    max_degree: int,
    sequence_length: int,
) -> list[tuple[tuple[int, int], dict[str, float]]]:
    if max_degree <= 0 or not pairs:
        return pairs
    degree: list[int] = [0] * (sequence_length + 1)
    accepted: list[tuple[tuple[int, int], dict[str, float]]] = []
    for pair, payload in pairs:
        left, right = pair
        if left <= sequence_length and right <= sequence_length:
            if degree[left] >= max_degree or degree[right] >= max_degree:
                continue
            degree[left] += 1
            degree[right] += 1
            accepted.append((pair, payload))
    return accepted


def _extract_contacts(
    coords: dict[int, tuple[float, float, float]],
    cutoff_angstrom: float,
    min_separation: int = 2,
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


def _rmsd(
    frame_coords: dict[int, tuple[float, float, float]],
    ref_coords: dict[int, tuple[float, float, float]],
    sequence_length: int,
) -> float:
    shared = [idx for idx in range(1, sequence_length + 1) if idx in frame_coords and idx in ref_coords]
    if not shared:
        return float("inf")
    total = 0.0
    for idx in shared:
        fx, fy, fz = frame_coords[idx]
        rx, ry, rz = ref_coords[idx]
        total += (fx - rx) ** 2 + (fy - ry) ** 2 + (fz - rz) ** 2
    return math.sqrt(total / len(shared))


def main() -> None:
    parser = argparse.ArgumentParser(description="Postprocess openmm trajectory contacts from last segment")
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--benchmark-file", default=str(Path(__file__).resolve().parents[1] / "data" / "folding_real_coordinate_visual_8.locked.json"))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_EXTERNAL_COUPLING_FILE))
    parser.add_argument("--trajectory-pdb", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--tail-fraction", type=float, default=0.20)
    parser.add_argument("--contact-cutoff-ang", type=float, default=7.0)
    parser.add_argument(
        "--frequency-threshold",
        default="auto",
        help='Either numeric threshold [0-1] or "auto" to use internal largest-gap cutoff.',
    )
    parser.add_argument("--chemical-threshold", type=float, default=0.50)
    parser.add_argument(
        "--dca-threshold",
        default="auto",
        help='Either numeric confidence [0-1] or "auto" to use internal largest-gap cutoff.',
    )
    parser.add_argument("--min-separation", type=int, default=2)
    parser.add_argument("--rmsd-cutoff", type=float, default=5.0)
    parser.add_argument("--native-rmsd-cutoff", action="store_true", help="Enable RMSD filtering against benchmark native coords.")
    parser.add_argument("--topology-mode", choices=("none", "interdomain", "intradomain"), default="none")
    parser.add_argument("--domain-boundaries", default="", help='Example: "1-120,121-214"')
    parser.add_argument("--max-degree", type=int, default=10, help="Residue degree cap (<=0 disables).")
    parser.add_argument("--control-count", type=int, default=0, help="Matched control count (0=>auto).")
    parser.add_argument("--max-control-overlap-fraction", type=float, default=0.20, help="Drop pair if appears in controls above this fraction.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    row = _load_row(Path(args.benchmark_file), args.source_accession)
    coupling_file = Path(args.external_coupling_file) if args.external_coupling_file else None
    dca_scores = _load_dca_scores(
        coupling_file=coupling_file,
        row=row,
        sequence_length=row.sequence_length,
    )
    trajectory_path = Path(args.trajectory_pdb)
    domain_boundaries = _parse_domain_boundaries(args.domain_boundaries)
    frames = _parse_ca_trajectory(trajectory_path)
    if not frames:
        raise SystemExit(f"no frames parsed from {trajectory_path}")

    tail_count = max(1, math.ceil(len(frames) * args.tail_fraction))
    final_frames = frames[-tail_count:]

    rmsd_refs = {
        int(point.residue_number): (float(point.x), float(point.y), float(point.z))
        for point in row.coordinate_points
    }

    considered_frames = final_frames
    if args.native_rmsd_cutoff and rmsd_refs:
        considered_frames = []
        for frame in final_frames:
            score = _rmsd(frame, rmsd_refs, row.sequence_length)
            if score <= args.rmsd_cutoff:
                considered_frames.append(frame)
        if not considered_frames:
            considered_frames = final_frames

    pair_counts: Counter[tuple[int, int]] = Counter()
    for frame in considered_frames:
        contacts = _extract_contacts(frame, cutoff_angstrom=args.contact_cutoff_ang, min_separation=args.min_separation)
        for left_right in contacts:
            pair_counts[left_right] += 1

    required_frames = len(considered_frames)
    candidate_pairs = tuple(sorted(pair_counts))
    if not pair_counts:
        payload = {
            "kind": "openmm_trajectory_contact_survivor_v0",
            "source_accession": row.source_accession,
            "row_id": row.row_id,
            "sequence_length": row.sequence_length,
            "frame_count": len(frames),
            "tail_fraction": args.tail_fraction,
            "tail_frames": tail_count,
            "considered_frames": required_frames,
            "frequency_threshold": args.frequency_threshold,
            "resolved_frequency_threshold": 0.0,
            "chemical_threshold": args.chemical_threshold,
            "dca_threshold": args.dca_threshold,
            "resolved_dca_threshold": 0.0,
            "topology_mode": args.topology_mode,
            "domain_boundaries": args.domain_boundaries,
            "max_degree": args.max_degree,
            "control_count": 0,
            "max_control_overlap_fraction": args.max_control_overlap_fraction,
            "kept_pair_count": 0,
            "predicted_contacts_from_selected_frames": 0,
            "raw_contact_pair_count": 0,
            "contact_cutoff_angstrom": args.contact_cutoff_ang,
            "rmsd_cutoff_angstrom": args.rmsd_cutoff,
            "used_native_rmsd_filter": bool(args.native_rmsd_cutoff),
            "native_contact_precision": 0.0,
            "native_contact_recall": 0.0,
            "short_range_recall": 0.0,
            "long_range_contact_recall": 0.0,
            "contact_map_f1": 0.0,
            "false_contact_rate": 0.0,
            "dca_supported_pair_count": 0,
            "control_removed_pair_count": 0,
            "raw_control_overlap_pair_count": 0,
            "degree_culled_pair_count": 0,
        }
        (out_dir / "openmm_stable_trajectory_contacts.csv").write_text("", encoding="utf-8")
        metric_path = out_dir / "openmm_trajectory_contact_metric.json"
        metric_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        print(json.dumps(payload, indent=2, sort_keys=True))
        return

    raw_pair_count = len(pair_counts)
    freq_values = [count / required_frames for count in pair_counts.values()]
    freq_threshold = _resolve_threshold(args.frequency_threshold, freq_values, default=0.50)
    dca_values = [dca_scores.get(pair, 0.0) for pair in candidate_pairs]
    dca_threshold = _resolve_threshold(args.dca_threshold, dca_values, default=0.0)

    candidates: list[tuple[tuple[int, int], dict[str, float]]] = []
    for (left, right), count in pair_counts.items():
        if left > len(row.sequence) or right > len(row.sequence):
            continue
        if left >= right:
            continue
        freq = count / required_frames
        if freq < freq_threshold:
            continue
        chem = chemical_score(row.sequence[left - 1], row.sequence[right - 1])
        if chem < args.chemical_threshold:
            continue
        dca = dca_scores.get((left, right), 0.0)
        if dca < dca_threshold:
            continue
        if not _topology_ok(
            left,
            right,
            domain_boundaries=domain_boundaries,
            topology_mode=args.topology_mode,
        ):
            continue
        candidates.append(
            (
                (left, right),
                {
                    "count": count,
                    "frequency": freq,
                    "chemical_score": chem,
                    "dca_score": dca,
                },
            )
        )

    pre_control_count = len(candidates)

    control_count = args.control_count
    if control_count <= 0:
        if pre_control_count <= 0:
            control_count = 0
        else:
            control_count = min(6, int(math.sqrt(pre_control_count)))
    control_hits = _control_gap(
        row=row,
        selected_pairs=[pair for pair, _ in candidates],
        candidate_pairs=candidate_pairs,
        control_count=control_count,
    )

    surviving_pairs: list[tuple[tuple[int, int], dict[str, float]]] = []
    control_removed = 0
    for pair, payload in candidates:
        overlap = control_hits.get(pair, 0) / control_count if control_count else 0.0
        payload = dict(payload)
        payload["control_overlap_fraction"] = overlap
        payload["control_gap"] = 1.0 - overlap
        if overlap > args.max_control_overlap_fraction:
            control_removed += 1
            continue
        payload["survival_score"] = (
            0.50 * payload["frequency"]
            + 0.25 * payload["chemical_score"]
            + 0.25 * payload["dca_score"]
        )
        surviving_pairs.append((pair, payload))

    surviving_pairs = sorted(
        surviving_pairs,
        key=lambda item: (-item[1]["survival_score"], -item[1]["frequency"], -item[1]["dca_score"], item[0][0], item[0][1]),
    )
    before_degree_cull = len(surviving_pairs)
    surviving_pairs = _apply_degree_cap(
        surviving_pairs,
        max_degree=args.max_degree,
        sequence_length=row.sequence_length,
    )
    degree_culled = before_degree_cull - len(surviving_pairs)

    kept_pairs = [pair for pair, _ in surviving_pairs]

    contact_file = out_dir / "openmm_stable_trajectory_contacts.csv"
    with contact_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "i",
                "j",
                "count",
                "frequency",
                "chemical_score",
                "dca_score",
                "control_overlap_fraction",
                "control_gap",
                "survival_score",
                "confidence_band",
            ],
        )
        writer.writeheader()
        for pair, row_item in surviving_pairs:
            writer.writerow(
                {
                    "i": pair[0],
                    "j": pair[1],
                    "count": row_item["count"],
                    "frequency": f"{row_item['frequency']:.6f}",
                    "chemical_score": f"{row_item['chemical_score']:.3f}",
                    "dca_score": f"{row_item['dca_score']:.6f}",
                    "control_overlap_fraction": f"{row_item['control_overlap_fraction']:.6f}",
                    "control_gap": f"{row_item['control_gap']:.6f}",
                    "survival_score": f"{row_item['survival_score']:.6f}",
                    "confidence_band": (
                        "high"
                        if row_item["frequency"] >= 0.90
                        else "medium" if row_item["frequency"] >= 0.60 else "low"
                    ),
                }
            )

    predicted_pairs = set(kept_pairs)
    metric = evaluate_contact_prediction(
        native_pairs=row.native_contact_pairs(),
        predicted_pairs=predicted_pairs,
        long_range_threshold=24,
        short_range_threshold=12,
    )
    metric_payload = metric.to_dict()

    payload = {
        "kind": "openmm_trajectory_contact_survivor_v0",
        "source_accession": row.source_accession,
        "row_id": row.row_id,
        "sequence_length": row.sequence_length,
        "frame_count": len(frames),
        "tail_fraction": args.tail_fraction,
        "tail_frames": tail_count,
        "considered_frames": required_frames,
        "frequency_threshold": args.frequency_threshold,
        "resolved_frequency_threshold": freq_threshold,
        "chemical_threshold": args.chemical_threshold,
        "dca_threshold": args.dca_threshold,
        "resolved_dca_threshold": dca_threshold,
        "contact_cutoff_angstrom": args.contact_cutoff_ang,
        "rmsd_cutoff_angstrom": args.rmsd_cutoff,
        "used_native_rmsd_filter": bool(args.native_rmsd_cutoff),
        "topology_mode": args.topology_mode,
        "domain_boundaries": args.domain_boundaries,
        "max_degree": args.max_degree,
        "control_count": control_count,
        "max_control_overlap_fraction": args.max_control_overlap_fraction,
        "raw_contact_pair_count": raw_pair_count,
        "pre_control_pair_count": pre_control_count,
        "control_removed_pair_count": control_removed,
        "degree_culled_pair_count": degree_culled,
        "dca_supported_pair_count": len([pair for pair in pair_counts if dca_scores.get(pair, 0.0) >= dca_threshold]),
        "kept_pair_count": len(predicted_pairs),
        "predicted_contacts_from_selected_frames": len(predicted_pairs),
        **metric_payload,
    }

    metric_path = out_dir / "openmm_trajectory_contact_metric.json"
    metric_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
