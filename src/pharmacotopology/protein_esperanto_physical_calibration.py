from __future__ import annotations

"""Real physical calibration inputs for Protein Esperanto.

These inputs do not make the coarse operator engine a physics simulator.  They
give the engine real coordinate-derived calibration context while preserving the
truth boundary: target-native contacts stay out of prediction.
"""

import json
from hashlib import sha256
from pathlib import Path
from statistics import mean
from typing import Any

from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows


CALIBRATION_KIND = "PROTEIN_ESPERANTO_REAL_PHYSICAL_CALIBRATION_INPUTS_v0"


def _stable_hash(data: Any) -> str:
    payload = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return sha256(payload).hexdigest()


def _observed_mean(values: list[float]) -> float:
    return round(mean(values), 6) if values else 0.0


def _contact_order(native_pairs: list[tuple[int, int]], sequence_length: int) -> float:
    if not native_pairs or sequence_length <= 0:
        return 0.0
    return round(sum(right - left for left, right in native_pairs) / (len(native_pairs) * sequence_length), 6)


def build_real_physical_calibration_inputs(
    benchmark_file: Path,
    *,
    excluded_target_accession: str | None = None,
) -> dict[str, Any]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    excluded = excluded_target_accession or None
    calibration_rows: list[dict[str, Any]] = []
    for row in rows:
        if excluded and row.source_accession == excluded:
            continue
        native_pairs = [(left, right) for left, right in row.native_contact_pairs()]
        calibration_rows.append({
            "row_id": row.row_id,
            "source_accession": row.source_accession,
            "source_database": row.source_database,
            "reference_structure_source": row.reference_structure_source,
            "reference_fold_class": row.reference_fold_class,
            "sequence_sha256": row.sequence_sha256,
            "sequence_length": row.sequence_length,
            "coordinate_trace_hash": row.coordinate_trace_hash,
            "native_contact_map_hash": row.native_contact_map_hash,
            "coordinate_coverage": row.coordinate_coverage,
            "native_contact_count": len(native_pairs),
            "native_contact_order": _contact_order(native_pairs, row.sequence_length),
        })
    fold_classes = sorted({row["reference_fold_class"] for row in calibration_rows})
    contact_counts = [float(row["native_contact_count"]) for row in calibration_rows]
    contact_orders = [float(row["native_contact_order"]) for row in calibration_rows]
    coverage = [float(row["coordinate_coverage"]) for row in calibration_rows]
    manifest = {
        "kind": CALIBRATION_KIND,
        "source_dataset": str(benchmark_file),
        "source_dataset_sha256": _stable_hash(json.loads(benchmark_file.read_text(encoding="utf-8"))),
        "source_coordinate_database": "RCSB_PDB",
        "calibration_input_type": "coordinate_derived_native_contact_summary",
        "row_count": len(calibration_rows),
        "excluded_target_accession": excluded,
        "target_native_excluded_from_calibration": True,
        "target_native_contacts_used_before_prediction": False,
        "coordinate_truth_used_as_prediction_input": False,
        "leave_one_target_out_calibration": True,
        "universal_physical_law_claim_allowed": False,
        "folding_problem_solved": False,
        "observable_families": [
            "coordinate_coverage",
            "native_contact_count",
            "native_contact_order",
            "native_contact_map_hash",
            "coordinate_trace_hash",
        ],
        "fold_class_coverage": fold_classes,
        "observed_contact_count_mean": _observed_mean(contact_counts),
        "observed_native_contact_order_mean": _observed_mean(contact_orders),
        "observed_coordinate_coverage_mean": _observed_mean(coverage),
        "calibration_rows": calibration_rows,
    }
    manifest["calibration_hash"] = _stable_hash({key: value for key, value in manifest.items() if key != "calibration_hash"})
    return manifest


def write_real_physical_calibration_inputs(
    benchmark_file: Path,
    output_file: Path,
    *,
    excluded_target_accession: str | None = None,
) -> dict[str, Any]:
    manifest = build_real_physical_calibration_inputs(
        benchmark_file,
        excluded_target_accession=excluded_target_accession,
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
