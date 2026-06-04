from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.artifact_io import write_csv_rows  # noqa: E402
from pharmacotopology.folding_evolutionary_constraints import (  # noqa: E402
    COUPLING_CONSTRAINT_KIND,
    EVOLUTIONARY_COUPLING_LAYER_KIND,
)
from pharmacotopology.folding_external_coupling_importer import (  # noqa: E402
    ExternalCouplingImportResult,
    import_external_coupling_dataset,
)
from pharmacotopology.folding_external_coupling_sources import (  # noqa: E402
    EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
    SERIOUS_EXTERNAL_COUPLING_POLICY,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_TARGET_MANIFEST = Path("data/external_coupling_target_manifest_v0.locked.json")
DEFAULT_OUTPUT = Path(
    "data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
)
DEFAULT_RUN_MANIFEST_OUTPUT = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_target_manifest_v0.json"
)
DEFAULT_BUILD_LOG_OUTPUT = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_build_log_v0.csv"
)


def _pdb_chain(source_accession: str) -> tuple[str, str]:
    pdb_id, chain_id = source_accession.split(":", 1)
    return pdb_id, chain_id


def target_manifest_rows(
    rows: tuple[RealCoordinateVisualRow, ...],
    *,
    status_by_row: Optional[dict[str, str]] = None,
) -> list[dict[str, object]]:
    status_by_row = status_by_row or {}
    output: list[dict[str, object]] = []
    for row in rows:
        pdb_id, chain_id = _pdb_chain(row.source_accession)
        output.append(
            {
                "row_id": row.row_id,
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "source_accession": row.source_accession,
                "sequence_length": row.sequence_length,
                "uniprot_accession": None,
                "sifts_mapping_available": False,
                "mapping_source_kind": "not_resolved_in_offline_v0",
                "mapping_sha256": None,
                "sequence_source_kind": "locked_rcsb_pdb_fasta_from_benchmark",
                "sequence_sha256": f"sha256:{row.sequence_sha256}",
                "msa_attempted": False,
                "msa_source_kind": "pfam_or_uniref_or_evcouplings",
                "external_coupling_status": status_by_row.get(row.row_id, "pending"),
            }
        )
    return output


def _build_status_from_import_status(status: object) -> str:
    row_status = getattr(status, "row_external_status")
    reason = getattr(status, "rejection_reason")
    raw_constraint_count = int(getattr(status, "raw_constraint_count"))
    if row_status == "external_couplings_available":
        return "external_couplings_available"
    if row_status == "external_couplings_rejected_coordinate_taint":
        return "external_couplings_rejected_coordinate_taint"
    if row_status == "external_couplings_rejected_low_coverage":
        return "external_couplings_rejected_low_coverage"
    if row_status == "external_couplings_rejected_mapping_ambiguous":
        return "external_couplings_rejected_position_mapping_ambiguous"
    if row_status == "external_couplings_rejected_low_depth":
        if raw_constraint_count == 0 or "no_external_couplings" in reason:
            return "external_couplings_rejected_no_sequence_mapping"
        return "external_couplings_rejected_low_msa_depth"
    return "external_couplings_rejected_tool_failed"


def _empty_build_log(rows: tuple[RealCoordinateVisualRow, ...]) -> list[dict[str, object]]:
    return [
        {
            "row_id": row.row_id,
            "source_accession": row.source_accession,
            "external_coupling_status": "external_couplings_rejected_no_sequence_mapping",
            "rejection_reason": "no_raw_external_coupling_file_supplied",
            "raw_constraint_count": 0,
            "accepted_constraint_count": 0,
            "duplicate_count_dropped": 0,
            "target_coverage": 0.0,
            "focus_sequence_mapping_confidence": 0.0,
            "effective_sequence_count_over_length": 0.0,
            "top_l_couplings_available": False,
            "raw_sequence_exposed": False,
        }
        for row in rows
    ]


def _build_log_from_import(
    result: ExternalCouplingImportResult,
) -> list[dict[str, object]]:
    return [
        {
            "row_id": status.row_id,
            "source_accession": status.source_accession,
            "external_coupling_status": _build_status_from_import_status(status),
            "rejection_reason": status.rejection_reason,
            "raw_constraint_count": status.raw_constraint_count,
            "accepted_constraint_count": status.accepted_constraint_count,
            "duplicate_count_dropped": 0,
            "target_coverage": status.target_coverage,
            "focus_sequence_mapping_confidence": (
                status.focus_sequence_mapping_confidence
            ),
            "effective_sequence_count_over_length": (
                status.effective_sequence_count_over_length
            ),
            "top_l_couplings_available": status.top_l_couplings_available,
            "raw_sequence_exposed": status.raw_sequence_exposed,
        }
        for status in result.row_statuses
    ]


def _constraints_from_import(
    result: ExternalCouplingImportResult,
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    available_rows = {
        status.row_id
        for status in result.row_statuses
        if status.row_external_status == "external_couplings_available"
    }
    for audit in result.constraint_audits:
        if audit.row_id not in available_rows:
            continue
        output.append(
            {
                "row_id": audit.row_id,
                "source_accession": audit.source_accession,
                "constraint_id": audit.constraint_id,
                "i": audit.i,
                "j": audit.j,
                "sequence_separation": audit.sequence_separation,
                "normalized_separation": audit.normalized_separation,
                "confidence": audit.confidence,
                "raw_score": audit.raw_score,
                "apc_corrected_score": audit.apc_corrected_score,
                "rank": audit.rank,
                "rank_fraction": audit.rank_fraction,
                "constraint_class": "external_dca_coupling",
                "source_kind": audit.source_kind,
                "msa_source_kind": audit.msa_source_kind,
                "msa_sha256": audit.msa_sha256,
                "msa_depth": audit.msa_depth,
                "effective_sequence_count": audit.effective_sequence_count,
                "effective_sequence_count_over_length": (
                    audit.effective_sequence_count_over_length
                ),
                "target_coverage": audit.target_coverage,
                "focus_sequence_mapping_confidence": (
                    audit.focus_sequence_mapping_confidence
                ),
                "coordinate_truth_used_to_build_constraint": (
                    audit.coordinate_truth_used_to_build_constraint
                ),
                "native_truth_used_before_coupling_selection": (
                    audit.native_truth_used_before_coupling_selection
                ),
                "structure_model_used": audit.structure_model_used,
                "raw_sequence_exposed": audit.raw_sequence_exposed,
            }
        )
    return output


def _write_json(path: Path, payload: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def build_real_external_coupling_file_v0(
    *,
    benchmark_file: Path,
    raw_external_coupling_file: Optional[Path],
    output: Path,
    run_manifest_output: Path,
    build_log_output: Path,
    target_manifest_path: Path,
) -> tuple[Path, Path, Path]:
    rows = tuple(load_real_coordinate_visual_rows(benchmark_file))
    result: Optional[ExternalCouplingImportResult] = None
    if raw_external_coupling_file is not None:
        result = import_external_coupling_dataset(
            rows=rows,
            external_coupling_file=raw_external_coupling_file,
            policy=SERIOUS_EXTERNAL_COUPLING_POLICY,
        )

    build_log = _empty_build_log(rows) if result is None else _build_log_from_import(result)
    status_by_row = {
        str(row["row_id"]): str(row["external_coupling_status"])
        for row in build_log
    }
    manifest = target_manifest_rows(rows, status_by_row=status_by_row)
    constraints = [] if result is None else _constraints_from_import(result)
    coupling_source_kind = (
        "external_msa_dca_plmc_v1"
        if result is None
        else result.dataset.coupling_source_kind
    )
    payload = {
        "layer_kind": EVOLUTIONARY_COUPLING_LAYER_KIND,
        "batch_id": EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
        "constraint_kind": COUPLING_CONSTRAINT_KIND,
        "coupling_source_kind": coupling_source_kind,
        "external_evolutionary_couplings_used": True,
        "coordinate_truth_used_to_build_constraints": (
            False if result is None else result.dataset.coordinate_truth_tainted
        ),
        "native_truth_used_before_coupling_selection": (
            False if result is None else result.dataset.native_truth_tainted
        ),
        "oracle_constraint_control": (
            False if result is None else result.dataset.oracle_constraint_control
        ),
        "raw_sequence_exposed": False,
        "source_benchmark_file": str(benchmark_file),
        "target_manifest_file": str(target_manifest_path),
        "benchmark_row_ids_preregistered": [row.row_id for row in rows],
        "reject_duplicate_coupling_pairs": True,
        "duplicate_count_dropped": 0,
        "constraints": constraints,
    }
    return (
        _write_json(output, payload),
        _write_json(run_manifest_output, manifest),
        write_csv_rows(build_log, build_log_output),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build REAL_EXTERNAL_COUPLING_FILE_BUILD_V0 artifacts. The script "
            "does not infer or fabricate couplings; without a raw external file "
            "it emits a zero-usable-row external coupling file and honest build log."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--raw-external-coupling-file",
        default=None,
        help="Optional provenance-locked raw external MSA/DCA coupling JSON.",
    )
    parser.add_argument("--target-manifest", default=str(DEFAULT_TARGET_MANIFEST))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument(
        "--run-manifest-output",
        default=str(DEFAULT_RUN_MANIFEST_OUTPUT),
    )
    parser.add_argument("--build-log-output", default=str(DEFAULT_BUILD_LOG_OUTPUT))
    args = parser.parse_args()

    raw_external = (
        None
        if args.raw_external_coupling_file is None
        else Path(args.raw_external_coupling_file)
    )
    outputs = build_real_external_coupling_file_v0(
        benchmark_file=Path(args.benchmark_file),
        raw_external_coupling_file=raw_external,
        output=Path(args.output),
        run_manifest_output=Path(args.run_manifest_output),
        build_log_output=Path(args.build_log_output),
        target_manifest_path=Path(args.target_manifest),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
