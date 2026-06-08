#!/usr/bin/env python3
from __future__ import annotations

"""Run MSA-free learned-model ensemble + physics/context refinement.

Inputs are already generated single-sequence structure predictions (PDB files).
This script does not call AlphaFold, MSAs, templates, or native coordinates
before contact selection.  Native contacts are attached only after selection for
benchmark audit.
"""

import argparse
import csv
import json
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_contact_law_features import contact_law_feature_rows_for_row
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_msa_free_model_ensemble import (
    MSAFreeStructureModelSpec,
    run_msa_free_model_ensemble,
    scan_model_pdb_specs,
)
from pharmacotopology.folding_nucleus_closure_search import nucleus_closure_events_for_row
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
DEFAULT_OUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "msa_free_model_ensemble_consensus_v0"
DEFAULT_SCAN_DIR = REPO_ROOT / "external_msa_free_predictors" / "tryhard_runs_live"


def _csv_write(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
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


def _row_by_accession(rows: Sequence[RealCoordinateVisualRow], accession: str) -> RealCoordinateVisualRow:
    matches = [row for row in rows if row.source_accession == accession]
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one row for {accession}, found {len(matches)}")
    return matches[0]


def _parse_model_pdb_arg(value: str, *, default_chain_id: str | None) -> MSAFreeStructureModelSpec:
    """Parse SOURCE_ID=PATH[#CHAIN] or PATH[#CHAIN]."""

    if "=" in value:
        source_id, rest = value.split("=", 1)
    else:
        rest = value
        source_id = Path(rest.split("#", 1)[0]).stem.replace(" ", "_")
    if "#" in rest:
        path_text, chain = rest.rsplit("#", 1)
        chain_id = chain or default_chain_id
    else:
        path_text = rest
        chain_id = default_chain_id
    return MSAFreeStructureModelSpec(
        source_id=source_id.strip() or Path(path_text).stem,
        pdb_path=str(Path(path_text).expanduser()),
        chain_id=chain_id,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="MSA-free learned-model ensemble consensus + physics refinement")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument(
        "--model-pdb",
        action="append",
        default=[],
        help="Add predicted PDB as SOURCE_ID=/path/model.pdb#CHAIN or /path/model.pdb#CHAIN. Repeatable.",
    )
    parser.add_argument(
        "--predicted-structure-dir",
        action="append",
        default=[],
        help="Directory to scan recursively for .pdb predictions.",
    )
    parser.add_argument("--default-chain", default=None)
    parser.add_argument("--allow-alphafold-source", action="store_true")
    parser.add_argument("--allow-msa-or-template-source", action="store_true")
    parser.add_argument("--iterations", type=int, default=4)
    parser.add_argument("--selection-threshold", type=float, default=0.62)
    parser.add_argument("--max-selected-contacts", type=int, default=None)
    parser.add_argument("--solved-precision-threshold", type=float, default=0.70)
    parser.add_argument("--solved-recall-threshold", type=float, default=0.70)
    parser.add_argument("--no-candidate-region-context", action="store_true")
    parser.add_argument("--no-coupling-context", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = tuple(load_real_coordinate_visual_rows(Path(args.benchmark_file)))
    row = _row_by_accession(rows, args.source_accession)

    specs: list[MSAFreeStructureModelSpec] = []
    for item in args.model_pdb:
        spec = _parse_model_pdb_arg(item, default_chain_id=args.default_chain)
        specs.append(
            MSAFreeStructureModelSpec(
                source_id=spec.source_id,
                pdb_path=spec.pdb_path,
                chain_id=spec.chain_id,
                allow_alphafold_source=args.allow_alphafold_source,
                allow_msa_or_template_source=args.allow_msa_or_template_source,
            )
        )
    scan_dirs = [Path(item).expanduser() for item in args.predicted_structure_dir]
    if not specs and not scan_dirs and DEFAULT_SCAN_DIR.exists():
        scan_dirs = [DEFAULT_SCAN_DIR]
    for directory in scan_dirs:
        specs.extend(
            scan_model_pdb_specs(
                root=directory,
                default_chain_id=args.default_chain,
                allow_alphafold_source=args.allow_alphafold_source,
                allow_msa_or_template_source=args.allow_msa_or_template_source,
            )
        )

    coupling_pairs = []
    if not args.no_coupling_context:
        coupling_dataset = load_coupling_dataset(Path(args.external_coupling_file))
        constraints = tuple(
            item for item in coupling_dataset.constraints if item.source_accession == row.source_accession
        )
        coupling_pairs = [constraint.pair() for constraint in constraints]

    candidate_region_pairs = []
    if not args.no_candidate_region_context:
        features = contact_law_feature_rows_for_row(row)
        events = nucleus_closure_events_for_row(row, features)
        for event in events:
            candidate_region_pairs.extend(event.candidate_region_pairs())

    packet = run_msa_free_model_ensemble(
        row=row,
        model_specs=tuple(specs),
        coupling_pairs=coupling_pairs,
        candidate_region_pairs=candidate_region_pairs,
        iteration_count=args.iterations,
        selection_threshold=args.selection_threshold,
        solved_precision_threshold=args.solved_precision_threshold,
        solved_recall_threshold=args.solved_recall_threshold,
        max_selected_contacts=args.max_selected_contacts,
    )

    payload = packet.to_dict()
    report_path = out_dir / "msa_free_model_ensemble_consensus_report.json"
    decisions_path = out_dir / "msa_free_model_ensemble_consensus_decisions.csv"
    selected_path = out_dir / "msa_free_model_ensemble_consensus_selected_contacts.csv"
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    _csv_write(decisions_path, [decision.to_dict() for decision in packet.decisions])
    _csv_write(selected_path, [decision.to_dict() for decision in packet.decisions if decision.selected])

    summary = {
        "kind": "msa_free_model_ensemble_consensus_summary_v0",
        "source_accession": packet.source_accession,
        "usable_model_count": packet.usable_model_count,
        "folding_problem_solved": packet.folding_problem_solved,
        "folding_solution_mode": packet.folding_solution_mode,
        "direct_best_model_source_id": packet.direct_best_model_source_id,
        "direct_precision": packet.direct_best_model_metric.native_contact_precision,
        "direct_recall": packet.direct_best_model_metric.native_contact_recall,
        "consensus_precision": packet.consensus_metric.native_contact_precision,
        "consensus_recall": packet.consensus_metric.native_contact_recall,
        "claim_rejection_reason": packet.claim_rejection_reason,
        "report_path": str(report_path),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
