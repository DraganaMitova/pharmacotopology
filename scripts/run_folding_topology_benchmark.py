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

from pharmacotopology.folding_metrics import (  # noqa: E402
    confusion_matrix,
    summarize_benchmark,
)
from pharmacotopology.folding_reference_loader import (  # noqa: E402
    FoldingReferenceDatasetValidation,
    load_folding_reference_dataset,
)
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
                "source_database": row["source_database"],
                "source_accession": row["source_accession"],
                "reference_structure_source": row["reference_structure_source"],
                "reference_label_source": row["reference_label_source"],
                "predicted_topology_signature": _signature_json(
                    row["predicted_topology_signature"]
                ),
                "reference_topology_signature": _signature_json(
                    row["reference_topology_signature"]
                ),
                "predicted_internal_fold_class": row[
                    "predicted_internal_fold_class"
                ],
                "reference_source_fold_class": row["reference_source_fold_class"],
                "predicted_fold_class": row["predicted_fold_class"],
                "reference_fold_class": row["reference_fold_class"],
                "contact_map_similarity": row["contact_map_similarity"],
                "fold_class_match": row["fold_class_match"],
                "uncertainty_radius": row["uncertainty_radius"],
                "evidence_readiness": row["evidence_readiness"],
                "failure_reason": row["failure_reason"],
                "is_external_reference": row["is_external_reference"],
                "reference_topology_signature_kind": row[
                    "reference_topology_signature_kind"
                ],
                "curation_notes": ";".join(
                    str(note) for note in row.get("curation_notes", ())
                ),
            }
        )
    return rows


def write_folding_benchmark_outputs(
    comparisons: Sequence[FoldingTopologyComparison],
    report_path: Path = DEFAULT_REPORT_PATH,
    csv_path: Path = DEFAULT_CSV_PATH,
    reference_dataset_validation: FoldingReferenceDatasetValidation | None = None,
    reference_dataset_metadata: Mapping[str, object] | None = None,
) -> tuple[Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    report = summarize_benchmark(comparisons)
    if reference_dataset_validation is not None:
        report["reference_dataset_validation"] = (
            reference_dataset_validation.to_dict()
        )
        if (
            reference_dataset_validation.valid
            and reference_dataset_validation.references_loaded > 0
            and reference_dataset_validation.placeholder_reference_count == 0
        ):
            report["benchmark_kind"] = "real_external_folding_topology_benchmark"
            report["practical_use"] = "external_folding_topology_alignment_review"
            report["external_validation_required"] = False
    if reference_dataset_metadata:
        for key in (
            "benchmark_sources",
            "target_benchmark_size",
            "locked_after_generation",
            "no_retuning_flag",
            "lock_certificate",
            "stratification_targets",
            "stratification_counts",
        ):
            if key in reference_dataset_metadata:
                report[key] = reference_dataset_metadata[key]
    report_path.write_text(
        json.dumps(
            report,
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


def _derived_output_path(report_path: Path, suffix: str) -> Path:
    name = report_path.name
    if name.endswith("_report.json"):
        return report_path.with_name(name.replace("_report.json", suffix))
    return report_path.with_name(report_path.stem + suffix)


def write_certificate_output(
    report: Mapping[str, object],
    output_path: Path,
) -> Path:
    payload = {
        "benchmark_kind": report.get("benchmark_kind"),
        "benchmark_size": report.get("benchmark_size"),
        "external_rows": report.get("external_rows"),
        "locked_after_generation": report.get("locked_after_generation"),
        "no_retuning_flag": report.get("no_retuning_flag"),
        "folding_solution_claim_created": report.get(
            "folding_solution_claim_created",
            False,
        ),
        "folding_problem_solved": report.get("folding_problem_solved", False),
        "lock_certificate": report.get("lock_certificate", {}),
        "reference_dataset_validation": report.get(
            "reference_dataset_validation",
            {},
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path


def write_failures_output(
    comparisons: Sequence[FoldingTopologyComparison],
    output_path: Path,
) -> Path:
    rows = [
        {
            "protein_id": comparison.protein_id,
            "source_database": comparison.source_database,
            "source_accession": comparison.source_accession,
            "reference_fold_class": comparison.reference_fold_class,
            "predicted_fold_class": comparison.predicted_fold_class,
            "predicted_internal_fold_class": comparison.predicted_internal_fold_class,
            "contact_map_similarity": comparison.contact_map_similarity,
            "uncertainty_radius": comparison.uncertainty_radius,
            "failure_reason": comparison.failure_reason,
        }
        for comparison in comparisons
        if comparison.failure_reason or not comparison.fold_class_match
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as file:
        fieldnames = [
            "protein_id",
            "source_database",
            "source_accession",
            "reference_fold_class",
            "predicted_fold_class",
            "predicted_internal_fold_class",
            "contact_map_similarity",
            "uncertainty_radius",
            "failure_reason",
        ]
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def write_confusion_matrix_output(
    comparisons: Sequence[FoldingTopologyComparison],
    output_path: Path,
) -> Path:
    matrix = confusion_matrix(comparisons)
    predicted_classes = sorted(
        {
            predicted
            for predictions in matrix.values()
            for predicted in predictions
        }
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=["reference_fold_class"] + predicted_classes,
        )
        writer.writeheader()
        for actual, predictions in matrix.items():
            row = {"reference_fold_class": actual}
            row.update(
                {
                    predicted: predictions.get(predicted, 0)
                    for predicted in predicted_classes
                }
            )
            writer.writerow(row)
    return output_path


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
    parser.add_argument(
        "--certificate-output",
        default="",
        help="Optional certificate JSON path. Defaults beside the report.",
    )
    parser.add_argument(
        "--failures-output",
        default="",
        help="Optional failures CSV path. Defaults beside the report.",
    )
    parser.add_argument(
        "--confusion-output",
        default="",
        help="Optional confusion-matrix CSV path. Defaults beside the report.",
    )
    parser.add_argument(
        "--benchmark-file",
        default="",
        help=(
            "Optional JSON file with externally derived folding reference rows. "
            "When omitted, bundled placeholder references are used."
        ),
    )
    parser.add_argument(
        "--require-external",
        action="store_true",
        help="Reject benchmark files with placeholder/example reference sources.",
    )
    args = parser.parse_args()

    validation = None
    references = None
    metadata = None
    if args.benchmark_file:
        try:
            dataset = load_folding_reference_dataset(
                Path(args.benchmark_file),
                require_external=args.require_external,
            )
        except ValueError as exc:
            parser.error(str(exc))
        references = dataset.references
        validation = dataset.validation
        metadata = dataset.metadata
    elif args.require_external:
        parser.error("--require-external requires --benchmark-file")

    comparisons = (
        run_folding_topology_benchmark(references)
        if references is not None
        else run_folding_topology_benchmark()
    )
    report_path, csv_path = write_folding_benchmark_outputs(
        comparisons,
        Path(args.report_output),
        Path(args.csv_output),
        reference_dataset_validation=validation,
        reference_dataset_metadata=metadata,
    )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    certificate_path = (
        Path(args.certificate_output)
        if args.certificate_output
        else _derived_output_path(report_path, "_certificate.json")
    )
    failures_path = (
        Path(args.failures_output)
        if args.failures_output
        else _derived_output_path(report_path, "_failures.csv")
    )
    confusion_path = (
        Path(args.confusion_output)
        if args.confusion_output
        else _derived_output_path(report_path, "_confusion_matrix.csv")
    )
    write_certificate_output(report, certificate_path)
    write_failures_output(comparisons, failures_path)
    write_confusion_matrix_output(comparisons, confusion_path)
    print(report_path)
    print(csv_path)
    print(certificate_path)
    print(failures_path)
    print(confusion_path)


if __name__ == "__main__":
    main()
