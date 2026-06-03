from __future__ import annotations

import json
from collections import Counter
from dataclasses import asdict, dataclass
from hashlib import sha256
from typing import Any, Mapping, Sequence

from pharmacotopology.folding_real_sources import (
    REAL_FOLDING_BENCHMARK_SOURCES,
    STRATIFIED_500_TARGETS,
    FoldingBenchmarkSource,
    source_catalog_to_dict,
)
from pharmacotopology.folding_topology import (
    FOLDING_TOPOLOGY_DIMENSIONS,
    FoldingReferenceExample,
    FoldingTopologyComparison,
    comparison_to_dict,
    signature_to_dict,
)


FOLD_CLASS_TO_STRATUM: dict[str, str] = {
    "alpha_rich": "mostly_alpha_or_compact",
    "beta_rich": "mostly_beta_or_long_range",
    "alpha_beta_mixed": "alpha_beta_mixed",
    "multidomain_boundary": "multidomain_or_boundary_sensitive",
    "disordered_flexible": "disordered_or_flexible",
    "compact_folded": "mostly_alpha_or_compact",
    "long_range_contact_rich": "mostly_beta_or_long_range",
    "mixed_uncertain": "alpha_beta_mixed",
    "boundary_sensitive": "multidomain_or_boundary_sensitive",
    "entanglement_rich": "multidomain_or_boundary_sensitive",
    "flexible_or_disordered": "disordered_or_flexible",
}


@dataclass(frozen=True)
class FoldingBenchmarkLockCertificate:
    benchmark_kind: str
    target_benchmark_size: int
    benchmark_size: int
    locked_after_generation: bool
    no_retuning_flag: bool
    dataset_hash: str
    recipe_commit_hash: str
    source_manifest_hash: str
    benchmark_sources: tuple[str, ...]
    stratification_targets: dict[str, int]
    stratification_counts: dict[str, int]
    lock_blockers: tuple[str, ...]
    clinical_use_allowed: bool = False
    drug_design_created: bool = False
    molecule_generated: bool = False
    protein_sequence_design_created: bool = False
    folding_solution_claim_created: bool = False
    folding_problem_solved: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def canonical_json(data: object) -> str:
    return json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_json(data: object) -> str:
    return "sha256:" + sha256(canonical_json(data).encode("utf-8")).hexdigest()


def reference_to_dict(reference: FoldingReferenceExample) -> dict[str, object]:
    return {
        "protein_id": reference.protein_id,
        "source_database": reference.source_database,
        "source_accession": reference.source_accession,
        "sequence": reference.sequence,
        "sequence_length": reference.sequence_length or len(reference.sequence),
        "reference_structure_source": reference.reference_structure_source,
        "reference_label_source": reference.reference_label_source,
        "reference_fold_class": reference.reference_fold_class,
        "reference_topology_signature_kind": (
            reference.reference_topology_signature_kind
        ),
        "reference_topology_signature": signature_to_dict(
            reference.reference_topology_signature
        ),
        "is_external_reference": reference.is_external_reference,
        "curation_notes": list(reference.curation_notes),
    }


def references_to_dicts(
    references: Sequence[FoldingReferenceExample],
) -> list[dict[str, object]]:
    return [reference_to_dict(reference) for reference in references]


def benchmark_source_ids_from_references(
    references: Sequence[FoldingReferenceExample],
    sources: Sequence[FoldingBenchmarkSource] = REAL_FOLDING_BENCHMARK_SOURCES,
) -> tuple[str, ...]:
    used: list[str] = []
    for source in sources:
        source_id = source.source_id.lower()
        for reference in references:
            source_database = reference.source_database.lower()
            source_fields = (
                reference.reference_structure_source.lower(),
                reference.reference_label_source.lower(),
                reference.source_accession.lower(),
            )
            if source_id in source_database:
                used.append(source.source_id)
                break
            if any(
                field.startswith(prefix)
                for field in source_fields
                for prefix in source.accepted_prefixes
            ):
                used.append(source.source_id)
                break
    return tuple(dict.fromkeys(used))


def fold_class_counts_from_references(
    references: Sequence[FoldingReferenceExample],
) -> dict[str, int]:
    return dict(
        sorted(Counter(reference.reference_fold_class for reference in references).items())
    )


def fold_class_counts_from_comparisons(
    comparisons: Sequence[FoldingTopologyComparison],
) -> dict[str, int]:
    return dict(
        sorted(Counter(comparison.reference_fold_class for comparison in comparisons).items())
    )


def stratum_for_fold_class(fold_class: str) -> str:
    return FOLD_CLASS_TO_STRATUM.get(fold_class, "alpha_beta_mixed")


def stratification_counts(
    references: Sequence[FoldingReferenceExample],
) -> dict[str, int]:
    counts = {key: 0 for key in STRATIFIED_500_TARGETS}
    for reference in references:
        counts[stratum_for_fold_class(reference.reference_fold_class)] += 1
    return counts


def stratification_deficits(
    counts: Mapping[str, int],
    targets: Mapping[str, int] = STRATIFIED_500_TARGETS,
) -> dict[str, int]:
    return {
        key: max(int(target) - int(counts.get(key, 0)), 0)
        for key, target in targets.items()
    }


def stratification_targets_for_size(target_size: int) -> dict[str, int]:
    if target_size <= 0:
        raise ValueError("target_size must be positive")
    strata = tuple(STRATIFIED_500_TARGETS)
    base = target_size // len(strata)
    remainder = target_size % len(strata)
    return {
        stratum: base + (1 if index < remainder else 0)
        for index, stratum in enumerate(strata)
    }


def confusion_matrix(
    comparisons: Sequence[FoldingTopologyComparison],
) -> dict[str, dict[str, int]]:
    matrix: dict[str, dict[str, int]] = {}
    for comparison in comparisons:
        actual = comparison.reference_fold_class
        predicted = comparison.predicted_fold_class
        matrix.setdefault(actual, {})
        matrix[actual][predicted] = matrix[actual].get(predicted, 0) + 1
    return {actual: dict(sorted(predicted.items())) for actual, predicted in sorted(matrix.items())}


def per_class_accuracy(
    comparisons: Sequence[FoldingTopologyComparison],
) -> dict[str, float]:
    totals: Counter[str] = Counter()
    matches: Counter[str] = Counter()
    for comparison in comparisons:
        totals[comparison.reference_fold_class] += 1
        if comparison.fold_class_match:
            matches[comparison.reference_fold_class] += 1
    return {
        fold_class: round(matches[fold_class] / total, 6)
        for fold_class, total in sorted(totals.items())
    }


def similarity_distribution(
    comparisons: Sequence[FoldingTopologyComparison],
) -> dict[str, int]:
    buckets = {
        "0.00-0.49": 0,
        "0.50-0.69": 0,
        "0.70-0.84": 0,
        "0.85-0.94": 0,
        "0.95-1.00": 0,
    }
    for comparison in comparisons:
        value = comparison.contact_map_similarity
        if value < 0.50:
            buckets["0.00-0.49"] += 1
        elif value < 0.70:
            buckets["0.50-0.69"] += 1
        elif value < 0.85:
            buckets["0.70-0.84"] += 1
        elif value < 0.95:
            buckets["0.85-0.94"] += 1
        else:
            buckets["0.95-1.00"] += 1
    return buckets


def worst_mismatches(
    comparisons: Sequence[FoldingTopologyComparison],
    *,
    limit: int = 25,
) -> list[dict[str, object]]:
    rows = [
        comparison_to_dict(comparison)
        for comparison in sorted(
            comparisons,
            key=lambda item: (item.fold_class_match, item.contact_map_similarity),
        )
        if not comparison.fold_class_match or comparison.contact_map_similarity < 0.70
    ]
    return rows[:limit]


def topology_signature_delta(row: Mapping[str, object]) -> dict[str, float]:
    predicted = row.get("predicted_topology_signature", {})
    reference = row.get("reference_topology_signature", {})
    if not isinstance(predicted, Mapping) or not isinstance(reference, Mapping):
        return {dimension: 0.0 for dimension in FOLDING_TOPOLOGY_DIMENSIONS}
    return {
        dimension: round(
            float(predicted.get(dimension, 0.0))
            - float(reference.get(dimension, 0.0)),
            6,
        )
        for dimension in FOLDING_TOPOLOGY_DIMENSIONS
    }


def lock_blockers(
    references: Sequence[FoldingReferenceExample],
    *,
    target_size: int,
    lock_requested: bool,
    stratification_targets: Mapping[str, int] | None = None,
) -> tuple[str, ...]:
    blockers: list[str] = []
    if not lock_requested:
        blockers.append("lock_not_requested")
    if len(references) == 0:
        blockers.append("no_external_reference_rows_loaded")
    if len(references) != target_size:
        blockers.append("benchmark_size_mismatch")
    targets = stratification_targets or stratification_targets_for_size(target_size)
    counts = stratification_counts(references)
    deficits = stratification_deficits(counts, targets)
    for stratum, deficit in deficits.items():
        if deficit:
            blockers.append(f"{stratum}_deficit:{deficit}")
    return tuple(blockers)


def build_lock_certificate(
    references: Sequence[FoldingReferenceExample],
    *,
    target_size: int,
    lock_requested: bool,
    recipe_commit_hash: str = "",
    sources: Sequence[FoldingBenchmarkSource] = REAL_FOLDING_BENCHMARK_SOURCES,
    stratification_targets: Mapping[str, int] | None = None,
) -> FoldingBenchmarkLockCertificate:
    reference_rows = references_to_dicts(references)
    source_manifest = source_catalog_to_dict(sources)
    targets = dict(stratification_targets or stratification_targets_for_size(target_size))
    blockers = lock_blockers(
        references,
        target_size=target_size,
        lock_requested=lock_requested,
        stratification_targets=targets,
    )
    locked = lock_requested and not blockers
    return FoldingBenchmarkLockCertificate(
        benchmark_kind="real_external_folding_topology_benchmark",
        target_benchmark_size=target_size,
        benchmark_size=len(references),
        locked_after_generation=locked,
        no_retuning_flag=locked,
        dataset_hash=sha256_json(reference_rows),
        recipe_commit_hash=recipe_commit_hash,
        source_manifest_hash=sha256_json(source_manifest),
        benchmark_sources=benchmark_source_ids_from_references(references, sources),
        stratification_targets=targets,
        stratification_counts=stratification_counts(references),
        lock_blockers=blockers,
    )


def build_locked_benchmark_payload(
    references: Sequence[FoldingReferenceExample],
    *,
    target_size: int,
    lock_requested: bool,
    recipe_commit_hash: str = "",
    sources: Sequence[FoldingBenchmarkSource] = REAL_FOLDING_BENCHMARK_SOURCES,
    stratification_targets: Mapping[str, int] | None = None,
) -> dict[str, Any]:
    certificate = build_lock_certificate(
        references,
        target_size=target_size,
        lock_requested=lock_requested,
        recipe_commit_hash=recipe_commit_hash,
        sources=sources,
        stratification_targets=stratification_targets,
    )
    return {
        "benchmark_kind": "real_external_folding_topology_benchmark",
        "benchmark_size": len(references),
        "target_benchmark_size": target_size,
        "external_validation_required": not certificate.locked_after_generation,
        "folding_solution_claim_created": False,
        "folding_problem_solved": False,
        "locked_after_generation": certificate.locked_after_generation,
        "no_retuning_flag": certificate.no_retuning_flag,
        "benchmark_sources": list(certificate.benchmark_sources),
        "stratification_targets": certificate.stratification_targets,
        "stratification_counts": certificate.stratification_counts,
        "lock_certificate": certificate.to_dict(),
        "boundary": {
            "clinical_use_allowed": False,
            "drug_design_created": False,
            "molecule_generated": False,
            "protein_sequence_design_created": False,
            "folding_solution_claim_created": False,
            "folding_problem_solved": False,
        },
        "references": references_to_dicts(references),
    }
