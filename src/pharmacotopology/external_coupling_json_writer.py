from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_evolutionary_constraints import (
    COUPLING_CONSTRAINT_KIND,
    EVOLUTIONARY_COUPLING_LAYER_KIND,
)
from pharmacotopology.folding_external_coupling_sources import (
    EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
)


def write_external_coupling_json(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    constraints: Sequence[Mapping[str, object]],
    output_path: Path,
    coupling_source_kind: str,
    source_benchmark_file: Path,
    build_metadata: Mapping[str, object],
) -> Path:
    payload = {
        "layer_kind": EVOLUTIONARY_COUPLING_LAYER_KIND,
        "batch_id": EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
        "constraint_kind": COUPLING_CONSTRAINT_KIND,
        "coupling_source_kind": coupling_source_kind,
        "external_evolutionary_couplings_used": bool(constraints),
        "external_coupling_build_attempted": True,
        "external_constraint_count": len(constraints),
        "coordinate_truth_used_to_build_constraints": False,
        "native_truth_used_before_coupling_selection": False,
        "oracle_constraint_control": False,
        "raw_sequence_exposed": False,
        "source_benchmark_file": str(source_benchmark_file),
        "benchmark_row_ids_preregistered": [row.row_id for row in rows],
        "constraints": list(constraints),
        **dict(build_metadata),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path
