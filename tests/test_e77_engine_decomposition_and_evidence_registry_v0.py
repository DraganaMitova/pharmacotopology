from __future__ import annotations

import json
from pathlib import Path

from pharmacotopology.protein_esperanto_engine import (
    E77_ABSTENTION_TAXONOMY,
    E77_REGISTRY_NAMES,
    evidence_registry_bundle,
)


ROOT = Path(__file__).resolve().parents[1]
E77_ROOT = ROOT / "data" / "protein_esperanto_engine" / "E77"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_e77_registry_bundle_externalizes_engine_ontology() -> None:
    bundle = evidence_registry_bundle()
    registries = bundle["registries"]
    token_registry = registries["token_evidence_registry"]

    assert bundle["engine_revision"] == "E77"
    assert bundle["baseline_engine_revision"] == "E76"
    assert bundle["registry_names"] == E77_REGISTRY_NAMES
    assert set(registries) == set(E77_REGISTRY_NAMES)
    assert bundle["abstention_taxonomy"] == E77_ABSTENTION_TAXONOMY
    assert bundle["complexity_budget"]["grammar_count"] == 19
    assert bundle["complexity_budget"]["operator_count"] == 30
    assert bundle["complexity_budget"]["sentence_clause_count"] == 7
    assert bundle["complexity_budget"]["token_count"] == token_registry["token_count"]
    assert token_registry["token_only_acceptance_forbidden"] is True
    assert all(row["allowed_use"] == "evidence_proposal_only" for row in token_registry["entries"])
    assert all(row["cannot_directly_accept"] is True for row in token_registry["entries"])
    assert {row["source"] for row in token_registry["entries"]}.issubset(
        {"definition_by_known_words", "usage_by_context", "manually_seeded_bootstrap"}
    )
    assert bundle["token_soup_firewall"]["token_hit_role"] == "evidence_proposal_only"
    assert bundle["token_soup_firewall"]["token_only_acceptance_count"] == 0
    assert bundle["physical_basis_claim_allowed"] is False
    assert bundle["protein_folding_solved"] is False


def test_e77_exported_registries_are_deterministic_and_non_accepting() -> None:
    cert = _read(E77_ROOT / "e77_engine_decomposition_and_evidence_registry_certificate.json")
    grammar = _read(E77_ROOT / "e77_grammar_registry.json")
    token = _read(E77_ROOT / "e77_token_evidence_registry.json")
    operator = _read(E77_ROOT / "e77_operator_registry.json")
    clauses = _read(E77_ROOT / "e77_sentence_clause_registry.json")

    assert cert["engine_revision"] == "E77"
    assert cert["registry_names"] == E77_REGISTRY_NAMES
    assert grammar["registry_hash"] == cert["registry_hashes"]["grammar_registry"]
    assert token["registry_hash"] == cert["registry_hashes"]["token_evidence_registry"]
    assert operator["registry_hash"] == cert["registry_hashes"]["operator_registry"]
    assert clauses["registry_hash"] == cert["registry_hashes"]["sentence_clause_registry"]
    assert token["token_only_acceptance_forbidden"] is True
    assert all(row["cannot_directly_accept"] is True for row in token["entries"])
    assert cert["token_soup_firewall"]["cannot_directly_accept"] is True
    assert cert["no_static_thresholds_used"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
