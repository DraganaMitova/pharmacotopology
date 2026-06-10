#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from pharmacotopology.external_coupling_json_writer import (
    write_external_coupling_json,
)
from pharmacotopology.external_dca_runner import (
    PfamDomainMapping,
    run_pfam_apc_covariation_for_row,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_ALIGNMENT_FILE = (
    REPO_ROOT / "external_msa" / "4ake_pfam00406" / "PF00406_full.sto"
)
DEFAULT_OUTPUT = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "4ake_query_centered_pfam00406_external_couplings.v0.locked.json"
)
DEFAULT_SOURCE_ACCESSION = "4AKE:A"
DEFAULT_PFAM_ID = "PF00406"
DEFAULT_DOMAIN_START = 5
DEFAULT_DOMAIN_END = 187
DEFAULT_MINIMUM_SEQUENCE_SEPARATION = 24
DEFAULT_MAX_RECORDS = 4000

COUPLING_SOURCE_KIND = "external_query_centered_pfam00406_mi_apc_v1"


def _constraint_rows(
    *,
    row: RealCoordinateVisualRow,
    pairs: Sequence,
    msa_sha256: str,
    msa_depth: int,
    effective_sequence_count: float,
) -> list[dict[str, object]]:
    selected = tuple(pairs)
    max_positive_score = max((max(0.0, pair.apc_corrected_score) for pair in selected), default=0.0)
    constraints: list[dict[str, object]] = []
    for rank, pair in enumerate(selected, start=1):
        confidence = (
            max(0.01, max(0.0, pair.apc_corrected_score) / max_positive_score)
            if max_positive_score
            else 0.01
        )
        constraints.append(
            {
                "row_id": row.row_id,
                "source_accession": row.source_accession,
                "constraint_id": f"pfam_query_centered_{row.row_id}_{pair.i}_{pair.j}_{rank}",
                "i": pair.i,
                "j": pair.j,
                "sequence_separation": pair.j - pair.i,
                "normalized_separation": round((pair.j - pair.i) / row.sequence_length, 6),
                "confidence": round(confidence, 6),
                "raw_score": pair.raw_score,
                "apc_corrected_score": pair.apc_corrected_score,
                "rank": rank,
                "rank_fraction": round(rank / row.sequence_length, 6),
                "constraint_class": "external_query_centered_pfam00406_apc_coupling",
                "source_kind": COUPLING_SOURCE_KIND,
                "msa_source_kind": "interpro_pfam_full_stockholm_query_stratified_mi_apc",
                "coordinate_truth_used_to_build_constraint": False,
                "native_truth_used_before_coupling_selection": False,
                "structure_model_used": False,
                "raw_sequence_exposed": False,
                "msa_sha256": msa_sha256,
                "msa_depth": msa_depth,
                "effective_sequence_count": float(effective_sequence_count),
                "effective_sequence_count_over_length": round(
                    effective_sequence_count / row.sequence_length,
                    6,
                ),
                "target_coverage": 1.0,
                "focus_sequence_mapping_confidence": 1.0,
                "top_L_couplings_available": len(selected) >= row.sequence_length,
                "pfam_id": pair.pfam_id,
            }
        )
    return constraints


def _row_by_source_accession(
    rows: Sequence[RealCoordinateVisualRow],
    source_accession: str,
) -> RealCoordinateVisualRow:
    for row in rows:
        if row.source_accession == source_accession:
            return row
    raise SystemExit(f"no row found for source_accession={source_accession!r}")


def build_query_centered_pfam00406_apc_couplings_v0(
    *,
    benchmark_file: Path,
    source_accession: str,
    alignment_file: Path,
    output: Path,
    pfam_id: str,
    domain_start: int,
    domain_end: int,
    max_records: int,
    minimum_sequence_separation: int,
) -> Path:
    rows = tuple(load_real_coordinate_visual_rows(benchmark_file))
    row = _row_by_source_accession(rows=rows, source_accession=source_accession)

    alignment_dir = Path("/tmp/pharmacotopology_query_centered_pfam00406")
    alignment_dir.mkdir(parents=True, exist_ok=True)
    alignment_cache_path = alignment_dir / f"{pfam_id}.full.sto.gz"
    alignment_cache_path.write_bytes(alignment_file.read_bytes())

    result = run_pfam_apc_covariation_for_row(
        row=row,
        mappings=(
            PfamDomainMapping(
                pfam_id=pfam_id,
                name=pfam_id,
                description=f"query-centered synthetic alignment domain from {pfam_id}",
                start=domain_start,
                end=domain_end,
                coverage=1.0,
            ),
        ),
        alignment_dir=alignment_dir,
        max_records=max_records,
        minimum_sequence_separation=minimum_sequence_separation,
        sample_strategy="query_stratified",
    )
    constraints = _constraint_rows(
        row=row,
        pairs=result.pairs,
        msa_sha256=result.msa_sha256,
        msa_depth=result.sample_depth,
        effective_sequence_count=float(result.sample_depth),
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    write_external_coupling_json(
        rows=(row,),
        constraints=constraints,
        output_path=output,
        coupling_source_kind=COUPLING_SOURCE_KIND,
        source_benchmark_file=benchmark_file,
        build_metadata={
            "build_batch_id": "query_centered_pfam00406_apc_couplings_v0",
            "dca_tool": "internal_query_stratified_pfam_apc",
            "covariation_method": "interpro_pfam_full_alignment_query_stratified_mi_apc",
            "sample_strategy": "query_stratified",
            "external_data_sources": (
                "Local alignment: InterPro/PDBe Pfam full Stockholm data (query-stratified sample)",
            ),
            "sample_depth": result.sample_depth,
            "total_depth_seen": result.total_depth_seen,
            "pfam_id": pfam_id,
            "domain_start": domain_start,
            "domain_end": domain_end,
        },
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build query-centered PF00406 MI/APC couplings from a local Pfam Stockholm alignment."
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--source-accession",
        default=DEFAULT_SOURCE_ACCESSION,
        help="target benchmark row source accession (default 4AKE:A)",
    )
    parser.add_argument(
        "--alignment-file",
        default=str(DEFAULT_ALIGNMENT_FILE),
        help="local Pfam full Stockholm for PF00406",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT),
        help="output coupling JSON path",
    )
    parser.add_argument("--pfam-id", default=DEFAULT_PFAM_ID)
    parser.add_argument("--domain-start", type=int, default=DEFAULT_DOMAIN_START)
    parser.add_argument("--domain-end", type=int, default=DEFAULT_DOMAIN_END)
    parser.add_argument("--max-records", type=int, default=DEFAULT_MAX_RECORDS)
    parser.add_argument(
        "--minimum-sequence-separation",
        type=int,
        default=DEFAULT_MINIMUM_SEQUENCE_SEPARATION,
    )
    args = parser.parse_args()

    benchmark_file = Path(args.benchmark_file)
    alignment_file = Path(args.alignment_file)
    output = Path(args.output)

    if not alignment_file.exists():
        raise SystemExit(f"alignment file does not exist: {alignment_file}")

    output = build_query_centered_pfam00406_apc_couplings_v0(
        benchmark_file=benchmark_file,
        source_accession=args.source_accession,
        alignment_file=alignment_file,
        output=output,
        pfam_id=args.pfam_id,
        domain_start=args.domain_start,
        domain_end=args.domain_end,
        max_records=args.max_records,
        minimum_sequence_separation=args.minimum_sequence_separation,
    )
    print(output)


if __name__ == "__main__":
    main()
