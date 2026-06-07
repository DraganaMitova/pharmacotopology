from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from pharmacotopology.folding_topology import (
    FOLDING_TOPOLOGY_DIMENSIONS,
    FoldingReferenceExample,
    KNOWN_FOLD_CLASSES,
    classify_broad_fold_topology,
    normalize_fold_class,
    normalize_sequence,
    reference_source_is_placeholder,
    reference_from_mapping,
    signature_to_dict,
)


EXTERNAL_REFERENCE_SOURCE_PREFIXES: tuple[str, ...] = (
    "pdb:",
    "afdb:",
    "alphafold:",
    "casp:",
    "cath:",
    "scop:",
    "disprot:",
    "external:",
)

@dataclass(frozen=True)
class FoldingReferenceDatasetValidation:
    dataset_path: str
    references_loaded: int
    external_reference_count: int
    placeholder_reference_count: int
    require_external: bool
    valid: bool
    violations: tuple[str, ...]
    warnings: tuple[str, ...]
    clinical_use_allowed: bool = False
    drug_design_created: bool = False
    molecule_generated: bool = False
    protein_sequence_design_created: bool = False
    folding_solution_claim_created: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class FoldingReferenceDataset:
    references: tuple[FoldingReferenceExample, ...]
    validation: FoldingReferenceDatasetValidation
    metadata: dict[str, Any]


def reference_source_is_external(source: str) -> bool:
    normalized = source.strip().lower()
    if not normalized:
        return False
    if reference_source_is_placeholder(source):
        return False
    return normalized.startswith(EXTERNAL_REFERENCE_SOURCE_PREFIXES)


def _json_rows(data: object) -> list[Mapping[str, object]]:
    raw_rows: object
    if isinstance(data, list):
        raw_rows = data
    elif isinstance(data, Mapping):
        raw_rows = data.get("references")
    else:
        raise ValueError("Benchmark file must contain a list or an object.")

    if not isinstance(raw_rows, list):
        raise ValueError("Benchmark file must include a references list.")
    rows: list[Mapping[str, object]] = []
    for index, row in enumerate(raw_rows, start=1):
        if not isinstance(row, Mapping):
            raise ValueError(f"Reference row {index} must be an object.")
        rows.append(row)
    return rows


def _load_json(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {path}: {exc}") from exc


def _reference_from_row(row: Mapping[str, object], index: int) -> FoldingReferenceExample:
    try:
        reference = reference_from_mapping(row)
        normalize_sequence(reference.sequence)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid folding reference row {index}: {exc}") from exc
    return reference


def load_folding_references(path: Path) -> tuple[FoldingReferenceExample, ...]:
    rows = _json_rows(_load_json(path))
    return tuple(
        _reference_from_row(row, index)
        for index, row in enumerate(rows, start=1)
    )


def validate_folding_references(
    references: Sequence[FoldingReferenceExample],
    *,
    dataset_path: Path,
    require_external: bool = False,
) -> FoldingReferenceDatasetValidation:
    violations: list[str] = []
    warnings: list[str] = []
    external_count = 0
    placeholder_count = 0

    if not references:
        violations.append("no_reference_rows_loaded")

    for index, reference in enumerate(references, start=1):
        source = reference.reference_structure_source
        if reference_source_is_external(source):
            external_count += 1
        else:
            placeholder_count += 1
            if require_external:
                violations.append(
                    f"row[{index}].reference_structure_source_not_external"
                )
        if require_external and not reference.is_external_reference:
            violations.append(f"row[{index}].is_external_reference_not_true")
        if require_external and not reference.source_database:
            violations.append(f"row[{index}].source_database_empty")
        if require_external and not reference.source_accession:
            violations.append(f"row[{index}].source_accession_empty")
        if require_external and not reference.reference_label_source:
            violations.append(f"row[{index}].reference_label_source_empty")

        normalized_sequence = normalize_sequence(reference.sequence)
        if (
            reference.sequence_length
            and reference.sequence_length != len(normalized_sequence)
        ):
            violations.append(f"row[{index}].sequence_length_mismatch")

        if reference.reference_fold_class not in KNOWN_FOLD_CLASSES:
            available = ",".join(KNOWN_FOLD_CLASSES)
            violations.append(
                f"row[{index}].unknown_reference_fold_class:{available}"
            )

        signature = signature_to_dict(reference.reference_topology_signature)
        if tuple(signature) != FOLDING_TOPOLOGY_DIMENSIONS:
            violations.append(f"row[{index}].topology_signature_schema_mismatch")
        for dimension, value in signature.items():
            if not 0.0 <= value <= 1.0:
                violations.append(
                    f"row[{index}].{dimension}_outside_unit_interval"
                )

        signature_kind = reference.reference_topology_signature_kind.lower()
        reference_class_from_signature = normalize_fold_class(
            classify_broad_fold_topology(reference.reference_topology_signature)
        )
        if (
            "prototype" not in signature_kind
            and reference_class_from_signature
            != normalize_fold_class(reference.reference_fold_class)
        ):
            warnings.append(
                f"row[{index}].reference_fold_class_signature_mismatch"
            )

    if require_external and external_count == 0:
        violations.append("no_external_reference_rows_loaded")

    return FoldingReferenceDatasetValidation(
        dataset_path=str(dataset_path),
        references_loaded=len(references),
        external_reference_count=external_count,
        placeholder_reference_count=placeholder_count,
        require_external=require_external,
        valid=not violations,
        violations=tuple(violations),
        warnings=tuple(warnings),
    )


def load_folding_reference_dataset(
    path: Path,
    *,
    require_external: bool = False,
) -> FoldingReferenceDataset:
    data = _load_json(path)
    references = tuple(
        _reference_from_row(row, index)
        for index, row in enumerate(_json_rows(data), start=1)
    )
    validation = validate_folding_references(
        references,
        dataset_path=path,
        require_external=require_external,
    )
    if not validation.valid:
        joined = "; ".join(validation.violations)
        raise ValueError(f"Invalid folding benchmark dataset: {joined}")
    metadata = dict(data) if isinstance(data, Mapping) else {}
    metadata.pop("references", None)
    return FoldingReferenceDataset(
        references=references,
        validation=validation,
        metadata=metadata,
    )
