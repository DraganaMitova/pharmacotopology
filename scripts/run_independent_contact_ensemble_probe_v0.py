from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Mapping, Sequence

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
from pharmacotopology.folding_independent_contact_evidence import (  # noqa: E402
    INDEPENDENT_CONTACT_ENSEMBLE_KIND,
    candidate_region_evidence_from_events,
    contact_evidence_from_predicted_pdb,
    coupling_contact_evidence,
    evaluate_ensemble_contacts,
    load_contact_evidence_json,
    native_coordinate_positive_control_evidence,
    write_contact_evidence_json,
)
from pharmacotopology.folding_native_contact_eval import contact_map_hash  # noqa: E402
from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_COUPLING_FILE = (
    REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
)
OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_REPORT_OUTPUT = OUTPUT_DIR / "independent_contact_ensemble_4ake_probe_v0.json"
DEFAULT_DECISIONS_OUTPUT = OUTPUT_DIR / "independent_contact_ensemble_4ake_decisions_v0.csv"
DEFAULT_SELECTED_OUTPUT = OUTPUT_DIR / "independent_contact_ensemble_4ake_selected_contacts_v0.csv"
DEFAULT_EVIDENCE_OUTPUT = OUTPUT_DIR / "independent_contact_ensemble_4ake_evidence_v0.json"


def _csv_write(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _row_by_accession(
    rows: Sequence[RealCoordinateVisualRow],
    source_accession: str,
) -> RealCoordinateVisualRow:
    for row in rows:
        if row.source_accession == source_accession:
            return row
    raise ValueError(f"source accession not found in benchmark: {source_accession}")


def _events_for_source(
    *,
    row: RealCoordinateVisualRow,
    context,
    event_source: str,
) -> tuple[NucleusClosureEvent, ...]:
    if event_source == "candidate_pool":
        return tuple(
            event for event in context.physical_context.candidate_events if event.row_id == row.row_id
        )
    if event_source == "competitive_pool":
        return tuple(event for event in context.competitive_events if event.row_id == row.row_id)
    if event_source == "aggressive_relative_frontier":
        return tuple(
            event
            for event in select_coupling_trace_loop_aggressive_relative_frontier_events(context)
            if event.row_id == row.row_id
        )
    if event_source == "self_deciding_generated_frontier":
        return tuple(
            event
            for event in select_coupling_trace_loop_self_deciding_frontier_generated_events(context)
            if event.row_id == row.row_id
        )
    raise ValueError(f"unknown event source: {event_source}")


def _family_counts(evidence) -> dict[str, int]:
    output: dict[str, set[tuple[int, int]]] = {}
    for item in evidence:
        output.setdefault(item.source_family, set()).add(item.pair())
    return {family: len(pairs) for family, pairs in sorted(output.items())}


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Probe whether an independent contact source can rescue a broad 4AKE "
            "candidate frontier without pretending native coordinates are a valid source."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument(
        "--event-source",
        default="aggressive_relative_frontier",
        choices=(
            "candidate_pool",
            "competitive_pool",
            "aggressive_relative_frontier",
            "self_deciding_generated_frontier",
        ),
    )
    parser.add_argument(
        "--independent-contact-json",
        action="append",
        default=[],
        help="JSON file containing independent contact evidence rows.",
    )
    parser.add_argument(
        "--predicted-pdb",
        default=None,
        help="External predicted PDB/AlphaFold/RoseTTAFold-style structure to convert to contacts.",
    )
    parser.add_argument("--predicted-pdb-chain", default=None)
    parser.add_argument("--predicted-source-id", default="external_predicted_structure")
    parser.add_argument(
        "--native-coordinate-positive-control",
        action="store_true",
        help=(
            "Build an explicit leakage positive-control source from the benchmark "
            "coordinates. This is not allowed as a benchmark claim."
        ),
    )
    parser.add_argument("--min-votes", type=int, default=2)
    parser.add_argument(
        "--allow-without-candidate-region",
        action="store_true",
        help="Allow independent + DCA contacts even if outside the current candidate region pool.",
    )
    parser.add_argument(
        "--allow-without-independent-structure",
        action="store_true",
        help="Allow candidate + DCA contacts without an independent structure source. Default is safer: require independent structure support.",
    )
    parser.add_argument("--report-output", default=str(DEFAULT_REPORT_OUTPUT))
    parser.add_argument("--decisions-output", default=str(DEFAULT_DECISIONS_OUTPUT))
    parser.add_argument("--selected-output", default=str(DEFAULT_SELECTED_OUTPUT))
    parser.add_argument("--evidence-output", default=str(DEFAULT_EVIDENCE_OUTPUT))
    args = parser.parse_args()

    rows = tuple(load_real_coordinate_visual_rows(Path(args.benchmark_file)))
    row = _row_by_accession(rows, args.source_accession)
    coupling_dataset = load_coupling_dataset(Path(args.external_coupling_file))
    constraints_by_row = coupling_dataset.constraints_by_row_id()
    row_constraints = tuple(constraints_by_row.get(row.row_id, ()))
    context = build_coupling_nucleus_context(rows=rows, coupling_dataset=coupling_dataset)
    events = _events_for_source(row=row, context=context, event_source=args.event_source)

    evidence = []
    evidence.extend(candidate_region_evidence_from_events(row=row, events=events))
    evidence.extend(coupling_contact_evidence(row=row, constraints=row_constraints))

    for evidence_path in args.independent_contact_json:
        loaded = load_contact_evidence_json(Path(evidence_path))
        evidence.extend(item for item in loaded if item.row_id == row.row_id)

    if args.predicted_pdb is not None:
        evidence.extend(
            contact_evidence_from_predicted_pdb(
                row=row,
                pdb_path=Path(args.predicted_pdb),
                chain_id=args.predicted_pdb_chain,
                source_id=args.predicted_source_id,
            )
        )

    if args.native_coordinate_positive_control:
        evidence.extend(native_coordinate_positive_control_evidence(row))

    report, decisions, _metric = evaluate_ensemble_contacts(
        row=row,
        evidence_sources=tuple(evidence),
        min_votes_required=args.min_votes,
        require_candidate_region_support=not args.allow_without_candidate_region,
        require_independent_structure_support=not args.allow_without_independent_structure,
    )
    selected_decisions = tuple(decision for decision in decisions if decision.selected)
    selected_pairs = tuple(decision.pair() for decision in selected_decisions)

    payload = {
        "kind": INDEPENDENT_CONTACT_ENSEMBLE_KIND,
        "row_id": row.row_id,
        "source_accession": row.source_accession,
        "event_source": args.event_source,
        "selected_event_count": len(events),
        "source_family_pair_counts": _family_counts(evidence),
        "selected_contact_map_hash": contact_map_hash(selected_pairs),
        "report": report.to_dict(),
        "selected_contacts": [decision.to_dict() for decision in selected_decisions],
        "safety": {
            "benchmark_claim_allowed": report.benchmark_claim_allowed,
            "claim_rejection_reason": report.claim_rejection_reason,
            "native_truth_attached_after_selection_for_evaluation": True,
            "coordinate_truth_used_before_selection": report.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": report.native_truth_used_before_selection,
            "raw_sequence_exposed": report.raw_sequence_exposed,
        },
    }

    report_path = Path(args.report_output)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _csv_write(Path(args.decisions_output), [decision.to_dict() for decision in decisions])
    _csv_write(Path(args.selected_output), [decision.to_dict() for decision in selected_decisions])
    write_contact_evidence_json(
        Path(args.evidence_output),
        evidence=tuple(evidence),
        source_kind="ensemble_probe_evidence_bundle",
    )

    print(json.dumps(payload["report"], indent=2, sort_keys=True))
    print(report_path)


if __name__ == "__main__":
    main()
