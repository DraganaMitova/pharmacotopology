from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "run_v41_protein_folding_claim_gate_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_v41_baseline_allows_c3_and_blocks_c5() -> None:
    runner = _load(RUNNER, "v41_runner_pass")
    cert = runner.build_v41()
    assert cert["control_status"] == "V41_MECHANISM_CLAIM_ALLOWED_C5_BLOCKED"
    assert cert["max_allowed_claim_level"] == "C3_FALSIFIABLE_OPERATOR_MECHANISM"
    assert cert["mechanism_claim_allowed"] is True
    assert cert["c5_claim_allowed"] is False
    assert cert["protein_folding_solved_claim_allowed"] is False
    assert cert["folding_problem_solved"] is False
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False


def test_v41_controls_all_pass_and_attack_surface_is_nonempty() -> None:
    runner = _load(RUNNER, "v41_runner_controls")
    cert = runner.build_v41()
    assert cert["control_count"] >= 12
    assert cert["passed_control_count"] == cert["control_count"]
    assert cert["scientist_attack_surface_count"] > 0
    assert cert["coordinate_derived_source_count_for_claim"] == 0
    assert cert["internal_runtime_source_count_for_claim"] == 0
    assert cert["native_metrics_used_for_claim"] is False
    assert cert["md_used_for_claim"] is False


def test_claim_text_is_bounded_and_forbidden_claims_are_explicit() -> None:
    runner = _load(RUNNER, "v41_runner_text")
    cert = runner.build_v41()
    allowed = cert["allowed_public_claim"].lower()
    assert "not solved protein folding" in allowed
    assert "not a de novo protein-structure predictor" in allowed
    assert "we solved protein folding" in cert["forbidden_public_claims"]
    assert "we outperform AlphaFold" in cert["forbidden_public_claims"]
    assert "this is a validated drug-discovery engine" in cert["forbidden_public_claims"]


def test_v41_blocks_prior_claim_allowed_true() -> None:
    runner = _load(RUNNER, "v41_runner_block_claim")
    certs = runner.load_certificates()
    mutated = copy.deepcopy(certs)
    mutated["V40"]["claim_allowed"] = True
    cert = runner.build_v41(mutated)
    assert cert["control_status"] == "V41_BLOCKED_CLAIM_BOUNDARY_VIOLATION"
    assert cert["mechanism_claim_allowed"] is False
    assert "all_prior_claim_allowed_flags_false" in cert["failed_checks"]


def test_v41_blocks_folding_problem_solved_true() -> None:
    runner = _load(RUNNER, "v41_runner_block_solved")
    certs = runner.load_certificates()
    mutated = copy.deepcopy(certs)
    mutated["V39"]["folding_problem_solved"] = True
    cert = runner.build_v41(mutated)
    assert cert["control_status"] == "V41_BLOCKED_CLAIM_BOUNDARY_VIOLATION"
    assert "all_prior_folding_problem_solved_flags_false" in cert["failed_checks"]


def test_v41_blocks_coordinate_counts_for_claim() -> None:
    runner = _load(RUNNER, "v41_runner_block_coord")
    certs = runner.load_certificates()
    mutated = copy.deepcopy(certs)
    mutated["V36"]["coordinate_derived_source_count"] = 1
    cert = runner.build_v41(mutated)
    assert cert["control_status"] == "V41_BLOCKED_CLAIM_BOUNDARY_VIOLATION"
    assert cert["coordinate_derived_source_count_for_claim"] > 0
    assert "no_coordinate_derived_sources_for_noncoordinate_claims" in cert["failed_checks"]


def test_writer_outputs_claim_gate_artifacts(tmp_path: Path) -> None:
    runner = _load(RUNNER, "v41_runner_writer")
    paths = runner.write_outputs(tmp_path / "out")
    for path in paths.values():
        assert path.exists()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["artifacts"]["claim_ladder"].endswith("claim_ladder.json")
    assert cert["artifacts"]["evidence_to_claim_matrix"].endswith("evidence_to_claim_matrix.csv")
    assert cert["max_allowed_claim_level"] == "C3_FALSIFIABLE_OPERATOR_MECHANISM"
    for artifact_key in [
        "claim_ladder",
        "evidence_to_claim_matrix",
        "allowed_claim_text",
        "forbidden_claim_text",
        "scientist_attack_surface",
    ]:
        assert Path(cert["artifacts"][artifact_key]).exists()
