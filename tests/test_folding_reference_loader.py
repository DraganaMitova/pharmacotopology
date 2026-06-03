import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from pharmacotopology.folding_reference_loader import (  # noqa: E402
    load_folding_reference_dataset,
    reference_source_is_external,
    validate_folding_references,
)
from pharmacotopology.folding_topology import (  # noqa: E402
    DEFAULT_FOLDING_BENCHMARKS,
    compare_to_reference,
    reference_from_mapping,
    run_folding_topology_benchmark,
)
from run_folding_topology_benchmark import (  # noqa: E402
    write_folding_benchmark_outputs,
)


def _external_reference_row() -> dict[str, object]:
    reference = DEFAULT_FOLDING_BENCHMARKS[0]
    return {
        "protein_id": "unit_test_external_alpha",
        "sequence": reference.sequence,
        "reference_structure_source": "pdb:UNITTEST_ALPHA",
        "reference_fold_class": reference.reference_fold_class,
        "reference_topology_signature": {
            "sequence_complexity": 0.61,
            "secondary_structure_balance": 0.74,
            "contact_map_closure": 0.72,
            "hydrophobic_core_closure": 0.65,
            "loop_disorder_pressure": 0.18,
            "domain_boundary_stability": 0.76,
            "long_range_contact_order": 0.48,
            "conformational_flexibility": 0.28,
            "knot_or_entanglement_signature": 0.18,
            "uncertainty_radius": 0.26,
        },
    }


def test_reference_source_external_detection_rejects_placeholders() -> None:
    assert reference_source_is_external("pdb:1ABC")
    assert reference_source_is_external("external:curated_contact_map")
    assert not reference_source_is_external("benchmark_placeholder:alpha")
    assert not reference_source_is_external("pdb:REPLACE_ME")


def test_template_source_stays_benchmark_shell_only() -> None:
    row = _external_reference_row()
    row["reference_structure_source"] = "pdb:REPLACE_ME"

    comparison = compare_to_reference(reference_from_mapping(row))

    assert comparison.evidence_readiness == "benchmark_shell_only"
    assert comparison.failure_reason == "external_structure_benchmark_not_attached"


def test_loader_accepts_external_benchmark_rows(tmp_path: Path) -> None:
    path = tmp_path / "real_folding_rows.json"
    path.write_text(
        json.dumps({"references": [_external_reference_row()]}),
        encoding="utf-8",
    )

    dataset = load_folding_reference_dataset(path, require_external=True)

    assert len(dataset.references) == 1
    assert dataset.validation.valid is True
    assert dataset.validation.external_reference_count == 1
    assert dataset.validation.placeholder_reference_count == 0
    assert dataset.validation.clinical_use_allowed is False
    assert dataset.validation.drug_design_created is False
    assert dataset.validation.protein_sequence_design_created is False


def test_loader_rejects_placeholder_when_external_required() -> None:
    validation = validate_folding_references(
        DEFAULT_FOLDING_BENCHMARKS[:1],
        dataset_path=Path("placeholder"),
        require_external=True,
    )

    assert validation.valid is False
    assert validation.external_reference_count == 0
    assert "row[1].reference_structure_source_not_external" in validation.violations
    assert "no_external_reference_rows_loaded" in validation.violations


def test_benchmark_report_includes_reference_dataset_validation(
    tmp_path: Path,
) -> None:
    path = tmp_path / "real_folding_rows.json"
    path.write_text(
        json.dumps({"references": [_external_reference_row()]}),
        encoding="utf-8",
    )
    dataset = load_folding_reference_dataset(path, require_external=True)
    comparisons = run_folding_topology_benchmark(dataset.references)
    report_path, _ = write_folding_benchmark_outputs(
        comparisons,
        tmp_path / "report.json",
        tmp_path / "rows.csv",
        reference_dataset_validation=dataset.validation,
    )
    report = json.loads(report_path.read_text(encoding="utf-8"))

    assert report["reference_dataset_validation"]["valid"] is True
    assert report["reference_dataset_validation"]["external_reference_count"] == 1
    assert report["clinical_use_allowed"] is False
    assert report["folding_problem_solved"] is False
