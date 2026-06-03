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

from build_real_folding_benchmark_500 import (  # noqa: E402
    build_real_folding_benchmark_payload,
    write_real_folding_benchmark_payload,
)
from pharmacotopology.folding_real_sources import (  # noqa: E402
    STRATIFIED_500_TARGETS,
    accepted_reference_prefixes,
    validate_source_catalog,
)
from pharmacotopology.folding_structure_features import (  # noqa: E402
    build_locked_benchmark_payload,
)


def test_real_folding_source_catalog_keeps_boundary_closed() -> None:
    review = validate_source_catalog()

    assert review.valid is True
    assert review.sources_reviewed >= 5
    assert "pdb:" in accepted_reference_prefixes()
    assert "disprot:" in accepted_reference_prefixes()
    assert review.clinical_use_allowed is False
    assert review.drug_design_created is False
    assert review.molecule_generated is False
    assert review.protein_sequence_design_created is False
    assert review.folding_solution_claim_created is False


def test_stratified_500_targets_are_balanced() -> None:
    assert sum(STRATIFIED_500_TARGETS.values()) == 500
    assert set(STRATIFIED_500_TARGETS.values()) == {100}


def test_empty_real_500_shell_is_not_locked_or_evidence() -> None:
    payload = build_locked_benchmark_payload(
        (),
        target_size=500,
        lock_requested=True,
        recipe_commit_hash="unit-test",
    )
    certificate = payload["lock_certificate"]

    assert payload["benchmark_size"] == 0
    assert payload["external_validation_required"] is True
    assert payload["locked_after_generation"] is False
    assert payload["folding_problem_solved"] is False
    assert certificate["lock_blockers"]
    assert "no_external_reference_rows_loaded" in certificate["lock_blockers"]
    assert payload["references"] == []


def test_builder_writes_empty_target_shell(tmp_path: Path) -> None:
    payload = build_real_folding_benchmark_payload(
        input_path=None,
        target_size=500,
        lock_requested=True,
        recipe_commit_hash="unit-test",
    )
    output = write_real_folding_benchmark_payload(
        payload,
        tmp_path / "folding_benchmarks_real_500.locked.json",
    )
    parsed = json.loads(output.read_text(encoding="utf-8"))

    assert parsed["status"] == "empty_locked_dataset_shell_no_external_rows_attached"
    assert parsed["locked_after_generation"] is False
    assert parsed["target_benchmark_size"] == 500
    assert parsed["boundary"]["folding_solution_claim_created"] is False


def test_tracked_real_500_shell_is_empty_and_blocked() -> None:
    path = ROOT / "data" / "folding_benchmarks_real_500.locked.json"
    parsed = json.loads(path.read_text(encoding="utf-8"))

    assert parsed["benchmark_size"] == 0
    assert parsed["references"] == []
    assert parsed["locked_after_generation"] is False
    assert parsed["external_validation_required"] is True
    assert parsed["folding_problem_solved"] is False
    assert "no_external_reference_rows_loaded" in (
        parsed["lock_certificate"]["lock_blockers"]
    )
