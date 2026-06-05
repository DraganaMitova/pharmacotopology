from __future__ import annotations

from collections import Counter
from typing import Any, Sequence

from pharmacotopology.folding_topology import (
    FoldingTopologyComparison,
    comparison_to_dict,
)


def mean_contact_map_similarity(
    comparisons: Sequence[FoldingTopologyComparison],
) -> float:
    if not comparisons:
        return 0.0
    return round(
        sum(comparison.contact_map_similarity for comparison in comparisons)
        / len(comparisons),
        6,
    )


def fold_class_match_rate(comparisons: Sequence[FoldingTopologyComparison]) -> float:
    if not comparisons:
        return 0.0
    return round(
        sum(1 for comparison in comparisons if comparison.fold_class_match)
        / len(comparisons),
        6,
    )


def match_count(comparisons: Sequence[FoldingTopologyComparison]) -> int:
    return sum(1 for comparison in comparisons if comparison.fold_class_match)


def mismatch_count(comparisons: Sequence[FoldingTopologyComparison]) -> int:
    return len(comparisons) - match_count(comparisons)


def external_reference_count(
    comparisons: Sequence[FoldingTopologyComparison],
) -> int:
    return sum(1 for comparison in comparisons if comparison.is_external_reference)


def mean_uncertainty_radius(
    comparisons: Sequence[FoldingTopologyComparison],
) -> float:
    if not comparisons:
        return 0.0
    return round(
        sum(comparison.uncertainty_radius for comparison in comparisons)
        / len(comparisons),
        6,
    )


def evidence_readiness_summary(
    comparisons: Sequence[FoldingTopologyComparison],
) -> dict[str, int]:
    return dict(
        sorted(
            Counter(comparison.evidence_readiness for comparison in comparisons).items()
        )
    )


def confusion_matrix(
    comparisons: Sequence[FoldingTopologyComparison],
) -> dict[str, dict[str, int]]:
    matrix: dict[str, dict[str, int]] = {}
    for comparison in comparisons:
        matrix.setdefault(comparison.reference_fold_class, {})
        predicted = comparison.predicted_fold_class
        matrix[comparison.reference_fold_class][predicted] = (
            matrix[comparison.reference_fold_class].get(predicted, 0) + 1
        )
    return {actual: dict(sorted(predicted.items())) for actual, predicted in sorted(matrix.items())}


def per_class_accuracy(
    comparisons: Sequence[FoldingTopologyComparison],
) -> dict[str, float]:
    totals: dict[str, int] = {}
    matches: dict[str, int] = {}
    for comparison in comparisons:
        actual = comparison.reference_fold_class
        totals[actual] = totals.get(actual, 0) + 1
        if comparison.fold_class_match:
            matches[actual] = matches.get(actual, 0) + 1
    return {
        actual: round(matches.get(actual, 0) / total, 6)
        for actual, total in sorted(totals.items())
    }


def benchmark_rows(
    comparisons: Sequence[FoldingTopologyComparison],
) -> list[dict[str, object]]:
    return [comparison_to_dict(comparison) for comparison in comparisons]


def summarize_benchmark(
    comparisons: Sequence[FoldingTopologyComparison],
) -> dict[str, Any]:
    failures = [
        comparison
        for comparison in comparisons
        if comparison.failure_reason
    ]
    perfect_matches = [
        comparison
        for comparison in comparisons
        if comparison.fold_class_match and comparison.contact_map_similarity >= 0.95
    ]
    return {
        "benchmark_kind": "protein_folding_topology_hypothesis_benchmark",
        "benchmark_size": len(comparisons),
        "practical_use": "folding_topology_benchmark_shell",
        "simulation_only": True,
        "hypothesis_numbers_only": True,
        "external_validation_required": True,
        "locked_after_generation": False,
        "no_retuning_flag": False,
        "benchmark_sources": [],
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "folding_solution_claim_created": False,
        "folding_problem_solved": False,
        "comparisons_reviewed": len(comparisons),
        "external_rows": external_reference_count(comparisons),
        "failure_count": len(failures),
        "perfect_matches": len(perfect_matches),
        "match_count": match_count(comparisons),
        "mismatch_count": mismatch_count(comparisons),
        "mismatches": mismatch_count(comparisons),
        "accuracy": fold_class_match_rate(comparisons),
        "mean_contact_map_similarity": mean_contact_map_similarity(comparisons),
        "fold_class_match_rate": fold_class_match_rate(comparisons),
        "per_class_accuracy": per_class_accuracy(comparisons),
        "confusion_matrix": confusion_matrix(comparisons),
        "mean_uncertainty_radius": mean_uncertainty_radius(comparisons),
        "evidence_readiness_summary": evidence_readiness_summary(comparisons),
        "boundary_statement": (
            "This benchmark shell tests whether abstract topology signatures "
            "from protein sequences can be compared with folding-relevant "
            "reference summaries. The default references are placeholders, "
            "not external validation."
        ),
        "next_calibration_steps": [
            "replace placeholder references with external structure benchmark rows",
            "derive reference topology signatures from contact maps and fold-class labels",
            "report confidence intervals across benchmark families",
            "keep molecule generation and clinical interpretation outside the boundary",
        ],
        "comparisons": benchmark_rows(comparisons),
    }
