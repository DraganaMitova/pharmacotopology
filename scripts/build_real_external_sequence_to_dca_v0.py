from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.artifact_io import write_csv_rows  # noqa: E402
from pharmacotopology.external_coupling_json_writer import (  # noqa: E402
    write_external_coupling_json,
)
from pharmacotopology.external_dca_runner import (  # noqa: E402
    PfamDomainMapping,
    run_pfam_apc_covariation_for_row,
)
from pharmacotopology.folding_external_coupling_sources import (  # noqa: E402
    REAL_EXTERNAL_SEQUENCE_TO_DCA_BUILD_BATCH_ID,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_OUTPUT = Path(
    "data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
)
DEFAULT_MAPPING_LOG = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_sequence_mapping_v0.csv"
)
DEFAULT_MSA_LOG = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_msa_build_log_v0.csv"
)
DEFAULT_DCA_LOG = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_dca_build_log_v0.csv"
)
DEFAULT_WORK_DIR = Path("/private/tmp/pharmacotopology_external_sequence_to_dca_v0")
PDBE_PFAM_MAPPING_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb_id}"
INTERPRO_PFAM_FULL_ALIGNMENT_URL = (
    "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam_id}/"
    "?annotation=alignment:full"
)
COUPLING_SOURCE_KIND = "external_evcouplings_sequence_covariation_v1"


def _download(url: str, output: Path) -> Path:
    if output.exists() and output.stat().st_size > 0:
        return output
    output.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=180) as response:
        output.write_bytes(response.read())
    return output


def _mapping_path(work_dir: Path, pdb_id: str) -> Path:
    return work_dir / f"pfam_{pdb_id.lower()}.json"


def _alignment_path(work_dir: Path, pfam_id: str) -> Path:
    return work_dir / f"{pfam_id}.full.sto.gz"


def _load_pdbe_pfam_mapping(
    *,
    row: RealCoordinateVisualRow,
    work_dir: Path,
) -> tuple[PfamDomainMapping, ...]:
    pdb_id, chain_id = row.source_accession.split(":", 1)
    path = _download(
        PDBE_PFAM_MAPPING_URL.format(pdb_id=pdb_id.lower()),
        _mapping_path(work_dir, pdb_id),
    )
    parsed = json.loads(path.read_text(encoding="utf-8"))
    entry = parsed.get(pdb_id.lower(), {})
    pfams = entry.get("Pfam", {}) if isinstance(entry, Mapping) else {}
    mappings: list[PfamDomainMapping] = []
    for pfam_id, info in pfams.items():
        if not isinstance(info, Mapping):
            continue
        for mapping in info.get("mappings", []):
            if not isinstance(mapping, Mapping):
                continue
            if mapping.get("chain_id") != chain_id and mapping.get("struct_asym_id") != chain_id:
                continue
            start = mapping.get("start", {})
            end = mapping.get("end", {})
            if not isinstance(start, Mapping) or not isinstance(end, Mapping):
                continue
            mappings.append(
                PfamDomainMapping(
                    pfam_id=str(pfam_id),
                    name=str(info.get("name", "")),
                    description=str(info.get("description", "")),
                    start=int(start["residue_number"]),
                    end=int(end["residue_number"]),
                    coverage=float(mapping.get("coverage", 0.0)),
                )
            )
    return tuple(mappings)


def _download_alignments(
    *,
    mappings_by_row: Mapping[str, Sequence[PfamDomainMapping]],
    work_dir: Path,
) -> None:
    pfam_ids = sorted(
        {
            mapping.pfam_id
            for mappings in mappings_by_row.values()
            for mapping in mappings
        }
    )
    for pfam_id in pfam_ids:
        _download(
            INTERPRO_PFAM_FULL_ALIGNMENT_URL.format(pfam_id=pfam_id),
            _alignment_path(work_dir, pfam_id),
        )


def _constraint_rows(
    *,
    row: RealCoordinateVisualRow,
    result: object,
) -> list[dict[str, object]]:
    pairs = tuple(getattr(result, "pairs"))
    constraints: list[dict[str, object]] = []
    top_l_available = len(pairs) >= row.sequence_length
    positive_scores = [max(0.0, float(pair.apc_corrected_score)) for pair in pairs]
    max_positive_score = max(positive_scores) if positive_scores else 0.0
    for rank, pair in enumerate(pairs, start=1):
        confidence = (
            max(0.01, max(0.0, float(pair.apc_corrected_score)) / max_positive_score)
            if max_positive_score
            else 0.01
        )
        constraints.append(
            {
                "row_id": row.row_id,
                "source_accession": row.source_accession,
                "constraint_id": f"pfam_apc_{row.row_id}_{pair.i}_{pair.j}_{rank}",
                "i": pair.i,
                "j": pair.j,
                "sequence_separation": pair.j - pair.i,
                "normalized_separation": round((pair.j - pair.i) / row.sequence_length, 6),
                "confidence": round(confidence, 6),
                "raw_score": pair.raw_score,
                "apc_corrected_score": pair.apc_corrected_score,
                "rank": rank,
                "rank_fraction": round(rank / row.sequence_length, 6),
                "constraint_class": "external_pfam_full_alignment_apc_covariation",
                "source_kind": COUPLING_SOURCE_KIND,
                "msa_source_kind": "interpro_pfam_full_stockholm_alignment_apc_covariation",
                "msa_sha256": getattr(result, "msa_sha256"),
                "msa_depth": getattr(result, "sample_depth"),
                "effective_sequence_count": float(getattr(result, "sample_depth")),
                "effective_sequence_count_over_length": round(
                    float(getattr(result, "sample_depth")) / row.sequence_length,
                    6,
                ),
                "target_coverage": 1.0,
                "focus_sequence_mapping_confidence": 1.0,
                "top_L_couplings_available": top_l_available,
                "coordinate_truth_used_to_build_constraint": False,
                "native_truth_used_before_coupling_selection": False,
                "structure_model_used": False,
                "raw_sequence_exposed": False,
                "pfam_id": pair.pfam_id,
            }
        )
    return constraints


def build_real_external_sequence_to_dca_v0(
    *,
    benchmark_file: Path,
    output: Path,
    mapping_log: Path,
    msa_log: Path,
    dca_log: Path,
    work_dir: Path,
    max_records: int,
) -> tuple[Path, Path, Path, Path]:
    rows = tuple(load_real_coordinate_visual_rows(benchmark_file))
    work_dir.mkdir(parents=True, exist_ok=True)
    mappings_by_row = {
        row.row_id: _load_pdbe_pfam_mapping(row=row, work_dir=work_dir)
        for row in rows
    }
    _download_alignments(mappings_by_row=mappings_by_row, work_dir=work_dir)

    mapping_rows: list[dict[str, object]] = []
    msa_rows: list[dict[str, object]] = []
    dca_rows: list[dict[str, object]] = []
    constraints: list[dict[str, object]] = []
    for row in rows:
        mappings = mappings_by_row[row.row_id]
        mapping_rows.append(
            {
                "row_id": row.row_id,
                "source_accession": row.source_accession,
                "mapping_status": (
                    "mapping_resolved"
                    if mappings
                    else "mapping_failed_no_sifts_uniprot_match"
                ),
                "failure_kind": "" if mappings else "mapping_failed",
                "rejection_reason": "" if mappings else "no_pdbe_pfam_mapping",
                "pfam_ids": ",".join(mapping.pfam_id for mapping in mappings),
                "mapping_source_kind": "pdbe_sifts_pfam_mapping_api",
                "coordinate_truth_used_to_build_constraints": False,
                "native_truth_used_before_coupling_selection": False,
                "structure_model_used": False,
                "raw_sequence_exposed": False,
            }
        )
        result = run_pfam_apc_covariation_for_row(
            row=row,
            mappings=mappings,
            alignment_dir=work_dir,
            max_records=max_records,
        )
        msa_rows.append(
            {
                "row_id": row.row_id,
                "source_accession": row.source_accession,
                "msa_status": "msa_attempted" if result.sample_depth else "msa_failed",
                "failure_kind": "" if result.sample_depth else result.failure_kind,
                "rejection_reason": "" if result.sample_depth else result.rejection_reason,
                "pfam_ids": result.pfam_ids,
                "msa_source_kind": "interpro_pfam_full_stockholm_alignment",
                "sample_depth": result.sample_depth,
                "total_depth_seen": result.total_depth_seen,
                "msa_sha256": result.msa_sha256,
                "raw_sequence_exposed": False,
            }
        )
        row_constraints = _constraint_rows(row=row, result=result)
        constraints.extend(row_constraints)
        dca_row = result.to_dict()
        dca_row["top_l_couplings_available"] = len(row_constraints) >= row.sequence_length
        dca_rows.append(dca_row)

    output_path = write_external_coupling_json(
        rows=rows,
        constraints=constraints,
        output_path=output,
        coupling_source_kind=COUPLING_SOURCE_KIND,
        source_benchmark_file=benchmark_file,
        build_metadata={
            "build_batch_id": REAL_EXTERNAL_SEQUENCE_TO_DCA_BUILD_BATCH_ID,
            "dca_tool_used": False,
            "dca_tool_missing": True,
            "covariation_fallback_used": True,
            "covariation_method": "interpro_pfam_full_alignment_mi_apc",
            "external_data_sources": (
                "PDBe SIFTS-derived Pfam mapping API",
                "InterPro Pfam full Stockholm alignments",
            ),
        },
    )
    return (
        output_path,
        write_csv_rows(mapping_rows, mapping_log),
        write_csv_rows(msa_rows, msa_log),
        write_csv_rows(dca_rows, dca_log),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build REAL_EXTERNAL_SEQUENCE_TO_DCA_BUILD_V0 artifacts from real "
            "external PDBe/Pfam/InterPro sequence-family data. If local plmc/"
            "ccmpred are unavailable, V0 emits an honest Pfam full-alignment "
            "MI/APC covariation channel rather than fabricating DCA output."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--mapping-log-output", default=str(DEFAULT_MAPPING_LOG))
    parser.add_argument("--msa-log-output", default=str(DEFAULT_MSA_LOG))
    parser.add_argument("--dca-log-output", default=str(DEFAULT_DCA_LOG))
    parser.add_argument("--work-dir", default=str(DEFAULT_WORK_DIR))
    parser.add_argument("--max-records", type=int, default=4000)
    args = parser.parse_args()

    outputs = build_real_external_sequence_to_dca_v0(
        benchmark_file=Path(args.benchmark_file),
        output=Path(args.output),
        mapping_log=Path(args.mapping_log_output),
        msa_log=Path(args.msa_log_output),
        dca_log=Path(args.dca_log_output),
        work_dir=Path(args.work_dir),
        max_records=args.max_records,
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
