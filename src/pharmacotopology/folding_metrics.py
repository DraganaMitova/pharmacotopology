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
    return {
        "benchmark_kind": "protein_folding_topology_hypothesis_benchmark",
        "practical_use": "folding_topology_benchmark_shell",
        "simulation_only": True,
        "hypothesis_numbers_only": True,
        "external_validation_required": True,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "folding_solution_claim_created": False,
        "folding_problem_solved": False,
        "comparisons_reviewed": len(comparisons),
        "failure_count": len(failures),
        "mean_contact_map_similarity": mean_contact_map_similarity(comparisons),
        "fold_class_match_rate": fold_class_match_rate(comparisons),
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
