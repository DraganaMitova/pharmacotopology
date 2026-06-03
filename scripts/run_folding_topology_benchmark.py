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

from pharmacotopology.folding_metrics import summarize_benchmark  # noqa: E402
from pharmacotopology.folding_topology import (  # noqa: E402
    FoldingTopologyComparison,
    comparison_to_dict,
    run_folding_topology_benchmark,
)


DEFAULT_REPORT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/folding_topology_benchmark_report.json"
)
DEFAULT_CSV_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/folding_topology_benchmark.csv"
)


def _signature_json(value: object) -> str:
    if isinstance(value, Mapping):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return str(value)


def _csv_rows(
    comparisons: Sequence[FoldingTopologyComparison],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for comparison in comparisons:
        row = comparison_to_dict(comparison)
        rows.append(
            {
                "protein_id": row["protein_id"],
                "sequence_length": row["sequence_length"],
                "reference_structure_source": row["reference_structure_source"],
                "predicted_topology_signature": _signature_json(
                    row["predicted_topology_signature"]
                ),
                "reference_topology_signature": _signature_json(
                    row["reference_topology_signature"]
                ),
                "predicted_fold_class": row["predicted_fold_class"],
                "reference_fold_class": row["reference_fold_class"],
                "contact_map_similarity": row["contact_map_similarity"],
                "fold_class_match": row["fold_class_match"],
                "uncertainty_radius": row["uncertainty_radius"],
                "evidence_readiness": row["evidence_readiness"],
                "failure_reason": row["failure_reason"],
            }
        )
    return rows


def write_folding_benchmark_outputs(
    comparisons: Sequence[FoldingTopologyComparison],
    report_path: Path = DEFAULT_REPORT_PATH,
    csv_path: Path = DEFAULT_CSV_PATH,
) -> tuple[Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(
            summarize_benchmark(comparisons),
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    rows = _csv_rows(comparisons)
    with csv_path.open("w", encoding="utf-8", newline="") as file:
        if rows:
            writer = csv.DictWriter(file, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        else:
            file.write("")

    return report_path, csv_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the protein folding topology benchmark shell. The default "
            "references are placeholders and do not validate folding claims."
        )
    )
    parser.add_argument(
        "--report-output",
        default=str(DEFAULT_REPORT_PATH),
        help="Path for the benchmark JSON report.",
    )
    parser.add_argument(
        "--csv-output",
        default=str(DEFAULT_CSV_PATH),
        help="Path for the benchmark CSV rows.",
    )
    args = parser.parse_args()

    comparisons = run_folding_topology_benchmark()
    report_path, csv_path = write_folding_benchmark_outputs(
        comparisons,
        Path(args.report_output),
        Path(args.csv_output),
    )
    print(report_path)
    print(csv_path)


if __name__ == "__main__":
    main()
