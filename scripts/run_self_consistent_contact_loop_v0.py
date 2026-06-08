#!/usr/bin/env python3
"""Run native-free self-consistent contact folding probe.

This is not an AlphaFold/PDB replacement claim.  It asks a narrower question:
can the internal coupling + closure field generate a stable self-source that
beats matched negative controls without reading native/coordinate truth?
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_coupling_nucleus_selector import (  # noqa: E402
    build_coupling_nucleus_context,
    select_coupling_trace_loop_aggressive_relative_frontier_events,
    select_coupling_trace_loop_self_deciding_frontier_generated_events,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset  # noqa: E402
from pharmacotopology.folding_independent_contact_evidence import write_contact_evidence_json  # noqa: E402
from pharmacotopology.folding_native_contact_eval import ContactPair  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_self_consistent_contact_loop import (  # noqa: E402
    SELF_CONSISTENT_CONTACT_LOOP_KIND,
    run_self_consistent_contact_loop,
)


DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_COUPLING_FILE = (
    REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
)
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "self_consistent_contact_loop_v0"


def _row_by_accession(rows: Sequence[RealCoordinateVisualRow], accession: str) -> RealCoordinateVisualRow:
    matches = [row for row in rows if row.source_accession == accession]
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one row for {accession}, found {len(matches)}")
    return matches[0]


def _csv_write(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _selected_rows(pairs: Sequence[ContactPair]) -> list[dict[str, object]]:
    return [
        {
            "i": left,
            "j": right,
            "sequence_separation": right - left,
        }
        for left, right in pairs
    ]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run a native-free self-consistent folding/contact loop and audit it "
            "against matched negative controls."
        )
    )
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument(
        "--event-source",
        default="aggressive_relative_frontier",
        choices=("aggressive_relative_frontier", "self_deciding_generated_frontier"),
    )
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = parser.parse_args()

    rows = tuple(load_real_coordinate_visual_rows(Path(args.benchmark_file)))
    row = _row_by_accession(rows, args.source_accession)
    coupling_dataset = load_coupling_dataset(Path(args.external_coupling_file))
    constraints_by_row = coupling_dataset.constraints_by_row_id()
    row_constraints = tuple(constraints_by_row.get(row.row_id, ()))
    context = build_coupling_nucleus_context(rows=rows, coupling_dataset=coupling_dataset)
    if args.event_source == "self_deciding_generated_frontier":
        all_events = select_coupling_trace_loop_self_deciding_frontier_generated_events(context)
    else:
        all_events = select_coupling_trace_loop_aggressive_relative_frontier_events(context)
    events = tuple(event for event in all_events if event.row_id == row.row_id)

    packet = run_self_consistent_contact_loop(
        row=row,
        events=events,
        constraints=row_constraints,
        event_source=args.event_source,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "self_consistent_contact_loop_report.json"
    iterations_path = out_dir / "self_consistent_contact_loop_iterations.csv"
    controls_path = out_dir / "self_consistent_contact_loop_controls.csv"
    selected_path = out_dir / "self_consistent_contact_loop_selected_contacts.csv"
    evidence_path = out_dir / "self_consistent_contact_loop_evidence.json"

    payload = {
        "kind": SELF_CONSISTENT_CONTACT_LOOP_KIND,
        "report": packet.report.to_dict(),
        "safety": {
            "coordinate_truth_used_before_selection": packet.report.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": packet.report.native_truth_used_before_selection,
            "native_truth_attached_after_selection_for_evaluation": True,
            "raw_sequence_exposed": packet.report.raw_sequence_exposed,
            "external_independent_claim_allowed": packet.report.external_independent_claim_allowed,
            "global_folding_claim_allowed": packet.report.global_folding_claim_allowed,
            "folding_problem_solved": packet.report.folding_problem_solved,
        },
        "selected_contact_map_hash": packet.report.final_contact_map_hash,
        "outputs": {
            "iterations": str(iterations_path),
            "controls": str(controls_path),
            "selected_contacts": str(selected_path),
            "evidence": str(evidence_path),
        },
    }
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _csv_write(iterations_path, [item.to_dict() for item in packet.iterations])
    _csv_write(controls_path, [item.to_dict() for item in packet.controls])
    _csv_write(selected_path, _selected_rows(packet.selected_pairs))
    write_contact_evidence_json(
        evidence_path,
        evidence=packet.self_evidence,
        source_kind="self_consistent_internal_evidence_bundle_v0",
    )

    print(json.dumps(packet.report.to_dict(), indent=2, sort_keys=True))
    print(report_path)


if __name__ == "__main__":
    main()
