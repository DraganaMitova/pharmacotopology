#!/usr/bin/env python3
"""Render a leakage-labeled visual proof GIF for the 4AKE AlphaFold ensemble run.

The GIF is not a scoring shortcut.  Native contacts are read only for the final
visual/evaluation frame after the ensemble selection has already been written.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Iterable, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from PIL import Image, ImageDraw, ImageFont  # noqa: E402

from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


DEFAULT_REPORT = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "independent_contact_ensemble_4ake_alphafold_v0.json"
)
DEFAULT_EVIDENCE = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "independent_contact_ensemble_4ake_alphafold_evidence_v0.json"
)
DEFAULT_SELECTED = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "independent_contact_ensemble_4ake_alphafold_selected_contacts_v0.csv"
)
DEFAULT_BENCHMARK = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "visual_proofs"
DEFAULT_GIF = DEFAULT_OUTPUT_DIR / "4ake_alphafold_ensemble_visual_proof_v0.gif"
DEFAULT_MANIFEST = DEFAULT_OUTPUT_DIR / "4ake_alphafold_ensemble_visual_proof_manifest_v0.json"

RGB_WHITE = (255, 255, 255)
RGB_TEXT = (22, 22, 22)
RGB_SUBTLE = (230, 230, 230)
RGB_GRID = (244, 244, 244)
RGB_CANDIDATE = (172, 172, 172)
RGB_DCA = (30, 105, 210)
RGB_INDEPENDENT = (235, 135, 35)
RGB_FINAL = (118, 60, 180)
RGB_TP = (20, 155, 75)
RGB_FP = (205, 45, 45)
RGB_FN = (160, 160, 160)
RGB_BLACK = (0, 0, 0)


def _repo_path(path: Path) -> Path:
    return path if path.is_absolute() else (REPO_ROOT / path)


def _safe_rel(path: Path) -> str:
    resolved = _repo_path(path).resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return str(resolved)


def _read_json(path: Path) -> Mapping[str, object]:
    return json.loads(_repo_path(path).read_text(encoding="utf-8"))


def _pair_from_mapping(row: Mapping[str, object]) -> tuple[int, int]:
    i = int(row["i"])
    j = int(row["j"])
    if j < i:
        i, j = j, i
    return i, j


def _load_evidence_pairs(path: Path) -> dict[str, set[tuple[int, int]]]:
    payload = _read_json(path)
    rows = payload.get("contacts", []) if isinstance(payload, Mapping) else []
    output: dict[str, set[tuple[int, int]]] = {
        "candidate_region": set(),
        "external_coupling": set(),
        "independent_structure": set(),
    }
    if not isinstance(rows, list):
        raise ValueError(f"evidence JSON has no contacts list: {path}")
    for item in rows:
        if not isinstance(item, Mapping):
            continue
        family = str(item.get("source_family", "unknown"))
        output.setdefault(family, set()).add(_pair_from_mapping(item))
    return output


def _load_selected_pairs(path: Path) -> set[tuple[int, int]]:
    path = _repo_path(path)
    if not path.exists() or path.stat().st_size == 0:
        return set()
    output: set[tuple[int, int]] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            output.add(_pair_from_mapping(row))
    return output


def _row_by_accession(benchmark_file: Path, source_accession: str):
    rows = load_real_coordinate_visual_rows(_repo_path(benchmark_file))
    for row in rows:
        if row.source_accession == source_accession:
            return row
    raise ValueError(f"source accession not found: {source_accession}")


def _project(pair: tuple[int, int], sequence_length: int, plot_left: int, plot_top: int, plot_size: int) -> tuple[int, int]:
    i, j = pair
    x = plot_left + round((j - 1) * (plot_size - 1) / max(1, sequence_length - 1))
    y = plot_top + round((i - 1) * (plot_size - 1) / max(1, sequence_length - 1))
    return x, y


def _draw_pairs(
    draw: ImageDraw.ImageDraw,
    *,
    pairs: Iterable[tuple[int, int]],
    sequence_length: int,
    plot_left: int,
    plot_top: int,
    plot_size: int,
    color: tuple[int, int, int],
    radius: int = 1,
) -> None:
    for pair in pairs:
        x, y = _project(pair, sequence_length, plot_left, plot_top, plot_size)
        if radius <= 1:
            draw.point((x, y), fill=color)
            draw.point((y - plot_top + plot_left, x - plot_left + plot_top), fill=color)
        else:
            draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=color)
            sx = y - plot_top + plot_left
            sy = x - plot_left + plot_top
            draw.ellipse((sx - radius, sy - radius, sx + radius, sy + radius), fill=color)


def _font(size: int = 14):
    try:
        return ImageFont.truetype("DejaVuSans.ttf", size)
    except OSError:
        return ImageFont.load_default()


def _draw_legend(draw: ImageDraw.ImageDraw, items: Sequence[tuple[str, tuple[int, int, int]]], x: int, y: int) -> None:
    font = _font(14)
    for label, color in items:
        draw.rectangle((x, y + 4, x + 16, y + 16), fill=color)
        draw.text((x + 24, y), label, fill=RGB_TEXT, font=font)
        y += 24


def _base_frame(
    *,
    title: str,
    subtitle: str,
    sequence_length: int,
    report: Mapping[str, object],
    safety: Mapping[str, object],
    size: int,
) -> tuple[Image.Image, ImageDraw.ImageDraw, int, int, int]:
    width = size
    height = size + 180
    img = Image.new("RGB", (width, height), RGB_WHITE)
    draw = ImageDraw.Draw(img)
    title_font = _font(22)
    subtitle_font = _font(14)
    small_font = _font(12)
    draw.text((24, 16), title, fill=RGB_TEXT, font=title_font)
    draw.text((24, 48), subtitle, fill=RGB_TEXT, font=subtitle_font)

    plot_left = 58
    plot_top = 92
    plot_size = size - 116
    draw.rectangle(
        (plot_left, plot_top, plot_left + plot_size, plot_top + plot_size),
        outline=RGB_BLACK,
        fill=RGB_WHITE,
    )
    for fraction in (0.25, 0.50, 0.75):
        gx = plot_left + int(plot_size * fraction)
        gy = plot_top + int(plot_size * fraction)
        draw.line((gx, plot_top, gx, plot_top + plot_size), fill=RGB_GRID)
        draw.line((plot_left, gy, plot_left + plot_size, gy), fill=RGB_GRID)
    draw.line((plot_left, plot_top, plot_left + plot_size, plot_top + plot_size), fill=RGB_SUBTLE)
    draw.text((plot_left, plot_top + plot_size + 8), f"residue index 1..{sequence_length}", fill=RGB_TEXT, font=small_font)

    metric_text = (
        f"LR precision={report.get('long_range_precision')}  "
        f"LR recall={report.get('long_range_recall')}  "
        f"contact precision={report.get('contact_precision')}  "
        f"contact recall={report.get('contact_recall')}"
    )
    safety_text = (
        f"claim={safety.get('benchmark_claim_allowed')}  "
        f"coordinate_truth_before_selection={safety.get('coordinate_truth_used_before_selection')}  "
        f"native_truth_before_selection={safety.get('native_truth_used_before_selection')}  "
        f"raw_sequence_exposed={safety.get('raw_sequence_exposed')}"
    )
    draw.text((24, height - 58), metric_text, fill=RGB_TEXT, font=small_font)
    draw.text((24, height - 34), safety_text, fill=RGB_TEXT, font=small_font)
    return img, draw, plot_left, plot_top, plot_size


def _make_frame(
    *,
    title: str,
    subtitle: str,
    sequence_length: int,
    report: Mapping[str, object],
    safety: Mapping[str, object],
    size: int,
    pair_layers: Sequence[tuple[str, Iterable[tuple[int, int]], tuple[int, int, int], int]],
    legend: Sequence[tuple[str, tuple[int, int, int]]],
) -> Image.Image:
    img, draw, plot_left, plot_top, plot_size = _base_frame(
        title=title,
        subtitle=subtitle,
        sequence_length=sequence_length,
        report=report,
        safety=safety,
        size=size,
    )
    for _label, pairs, color, radius in pair_layers:
        _draw_pairs(
            draw,
            pairs=pairs,
            sequence_length=sequence_length,
            plot_left=plot_left,
            plot_top=plot_top,
            plot_size=plot_size,
            color=color,
            radius=radius,
        )
    _draw_legend(draw, legend, plot_left, plot_top + plot_size + 34)
    return img


def render_visual_proof(
    *,
    report_path: Path,
    evidence_path: Path,
    selected_path: Path,
    benchmark_file: Path,
    source_accession: str,
    output_gif: Path,
    output_manifest: Path,
    output_frames_dir: Path,
    size: int = 720,
    duration_ms: int = 1150,
) -> Mapping[str, object]:
    report_payload = _read_json(report_path)
    report = report_payload.get("report", {})
    safety = report_payload.get("safety", {})
    if not isinstance(report, Mapping) or not isinstance(safety, Mapping):
        raise ValueError(f"probe report has invalid report/safety payload: {report_path}")

    row = _row_by_accession(benchmark_file, source_accession)
    native_pairs = set(row.native_contact_pairs())
    native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
    evidence_pairs = _load_evidence_pairs(evidence_path)
    selected_pairs = _load_selected_pairs(selected_path)
    selected_long = {pair for pair in selected_pairs if pair[1] - pair[0] >= 24}
    tp_long = selected_long & native_long
    fp_long = selected_long - native_long
    fn_long = native_long - selected_long

    frames: list[Image.Image] = []
    frames.append(
        _make_frame(
            title="4AKE visual proof: candidate region pool",
            subtitle="Broad sequence/candidate frontier. High ceiling, too noisy alone.",
            sequence_length=row.sequence_length,
            report=report,
            safety=safety,
            size=size,
            pair_layers=(("candidate", evidence_pairs.get("candidate_region", set()), RGB_CANDIDATE, 1),),
            legend=(("candidate region", RGB_CANDIDATE),),
        )
    )
    frames.append(
        _make_frame(
            title="4AKE visual proof: external DCA/coupling signal",
            subtitle="Sparse evolutionary coupling evidence. Precise but insufficient alone for 4AKE.",
            sequence_length=row.sequence_length,
            report=report,
            safety=safety,
            size=size,
            pair_layers=(("coupling", evidence_pairs.get("external_coupling", set()), RGB_DCA, 2),),
            legend=(("external coupling", RGB_DCA),),
        )
    )
    frames.append(
        _make_frame(
            title="4AKE visual proof: AlphaFold independent source",
            subtitle="Independent predicted-structure contacts. This is not native-coordinate leakage.",
            sequence_length=row.sequence_length,
            report=report,
            safety=safety,
            size=size,
            pair_layers=(("independent", evidence_pairs.get("independent_structure", set()), RGB_INDEPENDENT, 2),),
            legend=(("independent structure", RGB_INDEPENDENT),),
        )
    )
    frames.append(
        _make_frame(
            title="4AKE visual proof: ensemble/collapse selection",
            subtitle="Only contacts supported by candidate frontier plus independent structure survive.",
            sequence_length=row.sequence_length,
            report=report,
            safety=safety,
            size=size,
            pair_layers=(
                ("candidate", evidence_pairs.get("candidate_region", set()), (218, 218, 218), 1),
                ("selected", selected_pairs, RGB_FINAL, 2),
            ),
            legend=(("candidate background", (218, 218, 218)), ("final ensemble", RGB_FINAL)),
        )
    )
    frames.append(
        _make_frame(
            title="4AKE visual proof: post-selection evaluation",
            subtitle="Native truth is attached here only after selection: green=TP, red=FP, gray=missed native LR.",
            sequence_length=row.sequence_length,
            report=report,
            safety=safety,
            size=size,
            pair_layers=(
                ("missed native LR", fn_long, RGB_FN, 1),
                ("FP LR", fp_long, RGB_FP, 2),
                ("TP LR", tp_long, RGB_TP, 2),
            ),
            legend=(("TP long-range", RGB_TP), ("FP long-range", RGB_FP), ("missed native LR", RGB_FN)),
        )
    )

    output_gif = _repo_path(output_gif)
    output_manifest = _repo_path(output_manifest)
    output_frames_dir = _repo_path(output_frames_dir)
    output_frames_dir.mkdir(parents=True, exist_ok=True)
    frame_paths: list[str] = []
    for index, frame in enumerate(frames, start=1):
        path = output_frames_dir / f"4ake_visual_proof_frame_{index:02d}.png"
        frame.save(path)
        frame_paths.append(_safe_rel(path))

    output_gif.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        output_gif,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
        optimize=False,
    )

    manifest = {
        "kind": "4ake_alphafold_visual_proof_gif_v0",
        "source_accession": source_accession,
        "row_id": row.row_id,
        "sequence_length": row.sequence_length,
        "report": {
            "benchmark_claim_allowed": report.get("benchmark_claim_allowed"),
            "claim_rejection_reason": report.get("claim_rejection_reason"),
            "long_range_precision": report.get("long_range_precision"),
            "long_range_recall": report.get("long_range_recall"),
            "contact_precision": report.get("contact_precision"),
            "contact_recall": report.get("contact_recall"),
            "final_pair_count": report.get("final_pair_count"),
            "final_long_range_pair_count": report.get("final_long_range_pair_count"),
            "true_positive_long_range_contacts": report.get("true_positive_long_range_contacts"),
        },
        "safety": safety,
        "visual_boundary": {
            "native_truth_used_for_rendering_only_after_selection": True,
            "native_truth_used_before_selection": safety.get("native_truth_used_before_selection"),
            "coordinate_truth_used_before_selection": safety.get("coordinate_truth_used_before_selection"),
        },
        "pair_counts": {
            "candidate_region": len(evidence_pairs.get("candidate_region", set())),
            "external_coupling": len(evidence_pairs.get("external_coupling", set())),
            "independent_structure": len(evidence_pairs.get("independent_structure", set())),
            "selected_pairs": len(selected_pairs),
            "selected_long_range_pairs": len(selected_long),
            "true_positive_long_range_pairs": len(tp_long),
            "false_positive_long_range_pairs": len(fp_long),
            "missed_native_long_range_pairs": len(fn_long),
        },
        "outputs": {
            "gif": _safe_rel(output_gif),
            "frames": frame_paths,
        },
    }
    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    output_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Render the 4AKE AlphaFold ensemble visual proof GIF.")
    parser.add_argument("--report", default=str(DEFAULT_REPORT))
    parser.add_argument("--evidence", default=str(DEFAULT_EVIDENCE))
    parser.add_argument("--selected", default=str(DEFAULT_SELECTED))
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK))
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--output-gif", default=str(DEFAULT_GIF))
    parser.add_argument("--output-manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--output-frames-dir", default=str(DEFAULT_OUTPUT_DIR / "frames"))
    parser.add_argument("--size", type=int, default=720)
    parser.add_argument("--duration-ms", type=int, default=1150)
    args = parser.parse_args()

    manifest = render_visual_proof(
        report_path=Path(args.report),
        evidence_path=Path(args.evidence),
        selected_path=Path(args.selected),
        benchmark_file=Path(args.benchmark_file),
        source_accession=args.source_accession,
        output_gif=Path(args.output_gif),
        output_manifest=Path(args.output_manifest),
        output_frames_dir=Path(args.output_frames_dir),
        size=args.size,
        duration_ms=args.duration_ms,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
