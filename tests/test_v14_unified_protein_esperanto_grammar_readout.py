from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
for _path in (REPO_ROOT / "scripts", REPO_ROOT / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from run_v14_unified_protein_esperanto_grammar_readout_v0 import (  # noqa: E402
    _coherence_checks,
    _global_status,
    normalize_1cll,
    normalize_1ubq,
    normalize_4ake,
)


def test_normalize_1ubq_positive_adaptive_chemical_core() -> None:
    payload = {
        "target_role": "single_domain_compact",
        "chemical_policy": "adaptive_soft_guard",
        "claim_lock_status": "claim_locked_pending_cross_target_validation",
        "claim_lock_failed_checks": ["dca_background_enrichment_pass"],
        "selected_frequency_band": {
            "selected_pair_count": 1,
            "selected_balanced_core": [[23, 48]],
            "support_by_selected_pair": {"23-48": 9},
            "mean_frequency_by_selected_pair": {"23-48": 0.784444},
            "chemical_score_by_selected_pair": {"23-48": 0.1},
            "dca_score_by_selected_pair": {"23-48": 0.947233},
            "dca_absolute_support_pass": True,
            "dca_background_enrichment_ratio": 0.969389,
            "dca_background_enrichment_pass": False,
            "noise_added": 0,
            "long_range_evidence_polluted": False,
            "classification_coverage_ratio": 1.0,
            "legacy_chemical_hard_gate_would_block_selected_count": 1,
        },
    }
    row = normalize_1ubq(payload)
    assert row["protein"] == "1UBQ"
    assert row["positive_evidence_found"] is True
    assert row["selected_balanced_core"] == ["23-48"]
    assert row["claim_allowed"] is False
    assert "global_chemical_gate_rejected_as_too_strict" in row["final_status"]


def test_normalize_1cll_positive_c_domain_core_not_full_hinge() -> None:
    payload = {
        "target_role": "multi_domain_composite",
        "topology_policy": "hierarchical_domain_core_plus_interdomain",
        "claim_lock_status": "claim_lock_passed_but_claim_still_disabled",
        "claim_lock_failed_checks": [],
        "selected_C_domain_core": [[97, 133]],
        "selected_pair_count": 1,
        "selected_frequency_band": {
            "selected_C_domain_core": [[97, 133]],
            "selected_interdomain_hinge": [],
            "selected_pair_count": 1,
            "support_by_selected_pair": {"97-133": 7},
            "mean_frequency_by_selected_pair": {"97-133": 0.742},
            "chemical_score_by_selected_pair": {"97-133": 0.1},
            "dca_score_by_selected_pair": {"97-133": 1.0},
            "dca_absolute_support_pass": True,
            "dca_background_enrichment_ratio": 1.162021,
            "dca_background_enrichment_pass": True,
            "noise_added": 0,
            "long_range_evidence_polluted": False,
            "classification_coverage_ratio": 1.0,
        },
    }
    row = normalize_1cll(payload)
    assert row["protein"] == "1CLL"
    assert row["positive_evidence_found"] is True
    assert row["selected_C_domain_core"] == ["97-133"]
    assert row["selected_hinge_or_interdomain"] == []
    assert "interdomain_hinge_not_yet_proven" in row["final_status"]
    assert row["claim_allowed"] is False


def test_normalize_4ake_missing_is_honest_not_claimed() -> None:
    row = normalize_4ake(None)
    assert row["protein"] == "4AKE"
    assert row["artifact_status"] == "missing"
    assert row["positive_evidence_found"] is False
    assert row["claim_allowed"] is False


def test_normalize_4ake_positive_role_aware_artifact() -> None:
    payload = {
        "run_type": "GARAGE_ROLE_AWARE_RESCUE_SELECTOR_V9",
        "selected_strict_scaffold": [[1, 30]],
        "selected_balanced_core": [[20, 60]],
        "selected_border_rescue": [[100, 150]],
        "selected_pair_count": 3,
        "noise_added": 0,
        "long_range_evidence_polluted": False,
    }
    row = normalize_4ake(payload)
    assert row["positive_evidence_found"] is True
    assert row["selected_hinge_or_interdomain"] == ["100-150"]
    assert row["claim_allowed"] is False


def test_global_status_three_positive_claim_disabled() -> None:
    rows = [
        normalize_4ake({"selected_pair_count": 1, "selected_balanced_core": [[1, 30]]}),
        normalize_1ubq({"selected_frequency_band": {"selected_pair_count": 1, "selected_balanced_core": [[23, 48]]}}),
        normalize_1cll({"selected_frequency_band": {"selected_pair_count": 1, "selected_C_domain_core": [[97, 133]]}}),
    ]
    checks = _coherence_checks(rows)
    assert checks["no_claim_allowed_anywhere"] is True
    status = _global_status(rows, checks)
    assert status == "unified_role_aware_evidence_grammar_coherent_across_three_object_types_claim_disabled"


def test_global_status_partial_when_4ake_missing_but_1ubq_1cll_positive() -> None:
    rows = [
        normalize_4ake(None),
        normalize_1ubq({"selected_frequency_band": {"selected_pair_count": 1, "selected_balanced_core": [[23, 48]]}}),
        normalize_1cll({"selected_frequency_band": {"selected_pair_count": 1, "selected_C_domain_core": [[97, 133]]}}),
    ]
    checks = _coherence_checks(rows)
    status = _global_status(rows, checks)
    assert status == "partial_unified_grammar_panel_positive_on_1UBQ_1CLL_4AKE_artifact_missing_claim_disabled"
