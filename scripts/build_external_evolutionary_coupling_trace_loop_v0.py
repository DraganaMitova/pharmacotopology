from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.artifact_io import write_csv_rows  # noqa: E402
from pharmacotopology.folding_external_coupling_importer import (  # noqa: E402
    import_external_coupling_dataset,
    write_imported_external_coupling_dataset,
)
from pharmacotopology.folding_external_coupling_sources import (  # noqa: E402
    SERIOUS_EXTERNAL_COUPLING_POLICY,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_OUTPUT = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_imported.json"
)
DEFAULT_ROW_STATUS_OUTPUT = Path(
    "first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_row_status.csv"
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Validate and normalize a provenance-locked external MSA/DCA "
            "coupling file for EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_V0."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", required=True)
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument(
        "--row-status-output",
        default=str(DEFAULT_ROW_STATUS_OUTPUT),
    )
    args = parser.parse_args()

    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=Path(args.external_coupling_file),
        policy=SERIOUS_EXTERNAL_COUPLING_POLICY,
    )
    outputs = (
        write_imported_external_coupling_dataset(result, Path(args.output)),
        write_csv_rows(
            [status.to_dict() for status in result.row_statuses],
            Path(args.row_status_output),
        ),
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
