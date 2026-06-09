#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_native_contact_eval import evaluate_contact_prediction
from pharmacotopology.folding_template_docking import chemical_score


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
    parser.add_argument("--trajectory-pdb", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--tail-fraction", type=float, default=0.20)
    parser.add_argument("--contact-cutoff-ang", type=float, default=7.0)
    parser.add_argument("--frequency-threshold", type=float, default=0.50)
    parser.add_argument("--chemical-threshold", type=float, default=0.50)
    parser.add_argument("--min-separation", type=int, default=2)
    parser.add_argument("--rmsd-cutoff", type=float, default=5.0)
    parser.add_argument("--native-rmsd-cutoff", action="store_true", help="Enable RMSD filtering against benchmark native coords.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    row = _load_row(Path(args.benchmark_file), args.source_accession)
    trajectory_path = Path(args.trajectory_pdb)
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
    contact_rows = []
    for frame in considered_frames:
        contacts = _extract_contacts(frame, cutoff_angstrom=args.contact_cutoff_ang, min_separation=args.min_separation)
        for left_right in contacts:
            pair_counts[left_right] += 1

    required_frames = len(considered_frames)
    freq_cutoff = required_frames * args.frequency_threshold

    kept_pairs: list[tuple[int, int]] = []
    confidence_rows = []
    for (left, right), count in pair_counts.items():
        if count <= freq_cutoff:
            continue
        chem = chemical_score(row.sequence[left - 1], row.sequence[right - 1]) if left <= len(row.sequence) and right <= len(row.sequence) else 0.1
        if chem < args.chemical_threshold:
            continue
        freq = count / required_frames
        confidence_rows.append(
            {
                "i": left,
                "j": right,
                "count": count,
                "frequency": freq,
                "chemical_score": chem,
            }
        )
        kept_pairs.append((left, right))

    contact_file = out_dir / "openmm_stable_trajectory_contacts.csv"
    with contact_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["i", "j", "count", "frequency", "chemical_score", "confidence_band"],
        )
        writer.writeheader()
        for row_item in sorted(
            confidence_rows,
            key=lambda item: (-item["frequency"], -item["chemical_score"], item["i"], item["j"]),
        ):
            writer.writerow(
                {
                    "i": row_item["i"],
                    "j": row_item["j"],
                    "count": row_item["count"],
                    "frequency": f"{row_item['frequency']:.6f}",
                    "chemical_score": f"{row_item['chemical_score']:.3f}",
                    "confidence_band": (
                        "high" if row_item["frequency"] >= 0.90 else "medium"
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
        "source_accession": row.source_accession,
        "row_id": row.row_id,
        "sequence_length": row.sequence_length,
        "frame_count": len(frames),
        "tail_fraction": args.tail_fraction,
        "tail_frames": tail_count,
        "considered_frames": required_frames,
        "frequency_threshold": args.frequency_threshold,
        "chemical_threshold": args.chemical_threshold,
        "contact_cutoff_angstrom": args.contact_cutoff_ang,
        "rmsd_cutoff_angstrom": args.rmsd_cutoff,
        "used_native_rmsd_filter": bool(args.native_rmsd_cutoff),
        "kept_pair_count": len(predicted_pairs),
        "predicted_contacts_from_selected_frames": len(predicted_pairs),
        "raw_contact_pair_count": len(pair_counts),
        **metric_payload,
    }

    metric_path = out_dir / "openmm_trajectory_contact_metric.json"
    metric_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
