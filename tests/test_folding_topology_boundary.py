from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_metrics import summarize_benchmark
from pharmacotopology.folding_topology import (
    DEFAULT_FOLDING_BENCHMARKS,
    FOLDING_TOPOLOGY_DIMENSIONS,
    run_folding_topology_benchmark,
    signature_to_dict,
)


def test_folding_benchmark_report_keeps_safety_boundary() -> None:
    comparisons = run_folding_topology_benchmark()
    report = summarize_benchmark(comparisons)

    assert report["benchmark_kind"] == "protein_folding_topology_hypothesis_benchmark"
    assert report["simulation_only"] is True
    assert report["hypothesis_numbers_only"] is True
    assert report["external_validation_required"] is True
    assert report["clinical_use_allowed"] is False
    assert report["drug_design_created"] is False
    assert report["molecule_generated"] is False
    assert report["protein_sequence_design_created"] is False
    assert report["folding_solution_claim_created"] is False
    assert report["folding_problem_solved"] is False
    assert report["comparisons_reviewed"] == len(DEFAULT_FOLDING_BENCHMARKS)


def test_folding_benchmark_rows_include_required_fields() -> None:
    report = summarize_benchmark(run_folding_topology_benchmark())
    required = {
        "protein_id",
        "sequence_length",
        "reference_structure_source",
        "predicted_topology_signature",
        "reference_topology_signature",
        "contact_map_similarity",
        "fold_class_match",
        "uncertainty_radius",
        "evidence_readiness",
        "failure_reason",
    }

    for row in report["comparisons"]:
        assert required.issubset(row)
        assert row["evidence_readiness"] == "benchmark_shell_only"
        assert row["failure_reason"] == "external_structure_benchmark_not_attached"
        assert row["clinical_use_allowed"] is False
        assert row["drug_design_created"] is False
        assert row["molecule_generated"] is False
        assert row["protein_sequence_design_created"] is False


def test_folding_signatures_use_bounded_dimension_schema() -> None:
    comparison = run_folding_topology_benchmark()[0]
    predicted = signature_to_dict(comparison.predicted_topology_signature)
    reference = signature_to_dict(comparison.reference_topology_signature)

    assert tuple(predicted) == FOLDING_TOPOLOGY_DIMENSIONS
    assert tuple(reference) == FOLDING_TOPOLOGY_DIMENSIONS
    assert all(0.0 <= value <= 1.0 for value in predicted.values())
    assert all(0.0 <= value <= 1.0 for value in reference.values())
