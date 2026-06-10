from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
for _path in (REPO_ROOT / "scripts", REPO_ROOT / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from run_v15_dynamic_separation_grammar_readout_v0 import (  # noqa: E402
    _coherence_checks,
    _global_status,
    normalize_1cll_dynamic,
    normalize_1ubq_dynamic,
    normalize_4ake_dynamic,
)


def test_1ubq_dynamic_separation_keeps_selected_core_without_fixed_cutoff() -> None:
    payload = {
        "target_role": "single_domain_compact",
        "chemical_policy": "adaptive_soft_guard",
        "sequence_length": 76,
        "selected_frequency_band": {
            "selected_pair_count": 1,
            "selected_balanced_core": [[23, 48]],
            "support_by_selected_pair": {"23-48": 9},
            "mean_frequency_by_selected_pair": {"23-48": 0.784444},
            "chemical_score_by_selected_pair": {"23-48": 0.1},
            "dca_score_by_selected_pair": {"23-48": 0.947233},
            "noise_added": 0,
            "long_range_evidence_polluted": False,
            "classification_coverage_ratio": 1.0,
        },
    }
    row = normalize_1ubq_dynamic(payload)
    assert row["positive_evidence_found"] is True
    assert row["selected_pairs"] == ["23-48"]
    assert row["fixed_residue_cutoff_used"] is False
    role = row["dynamic_pair_roles"]["23-48"]
    assert role["domain_relation"] == "single_domain"
    assert role["separation_filter_applied"] is False
    assert role["fixed_residue_cutoff_used"] is False
    assert role["evidence_class"] == "single_domain_compact_core_evidence"


def test_1cll_dynamic_separation_preserves_c_domain_core_role() -> None:
    payload = {
        "target_role": "multi_domain_composite",
        "topology_policy": "hierarchical_domain_core_plus_interdomain",
        "sequence_length": 148,
        "selected_C_domain_core": [[97, 133]],
        "selected_pair_count": 1,
        "role_by_candidate_pair": {
            "97-133": {
                "domain_relation": "intradomain_C",
                "topology_role": "C_domain_compact_core_candidate",
                "evidence_class": "C_domain_core_evidence",
                "tail_frequency_mean": 0.742,
                "chemical_score": 0.1,
                "dca_score": 1.0,
            }
        },
        "selected_frequency_band": {
            "selected_C_domain_core": [[97, 133]],
            "selected_pair_count": 1,
            "support_by_selected_pair": {"97-133": 7},
            "mean_frequency_by_selected_pair": {"97-133": 0.742},
            "chemical_score_by_selected_pair": {"97-133": 0.1},
            "dca_score_by_selected_pair": {"97-133": 1.0},
            "noise_added": 0,
            "long_range_evidence_polluted": False,
            "classification_coverage_ratio": 1.0,
        },
    }
    row = normalize_1cll_dynamic(payload)
    assert row["positive_evidence_found"] is True
    assert row["selected_C_domain_core"] == ["97-133"]
    role = row["dynamic_pair_roles"]["97-133"]
    assert role["domain_relation"] == "intradomain_C"
    assert role["evidence_class"] == "C_domain_core_evidence"
    assert role["fixed_residue_cutoff_used"] is False
    assert "interdomain_hinge_not_yet_proven" in row["final_status"]


def test_4ake_missing_remains_claim_disabled() -> None:
    row = normalize_4ake_dynamic()
    assert row["artifact_status"] == "missing_machine_readable_grammar_artifact"
    assert row["fixed_residue_cutoff_used"] is False
    assert row["claim_allowed"] is False


def test_global_status_partial_with_1ubq_1cll_positive_and_4ake_missing() -> None:
    rows = [
        normalize_4ake_dynamic(),
        normalize_1ubq_dynamic({"selected_frequency_band": {"selected_pair_count": 1, "selected_balanced_core": [[23, 48]]}}),
        normalize_1cll_dynamic({"selected_frequency_band": {"selected_pair_count": 1, "selected_C_domain_core": [[97, 133]]}}),
    ]
    checks = _coherence_checks(rows)
    assert checks["no_fixed_residue_cutoff_used_anywhere"] is True
    assert checks["no_claim_allowed_anywhere"] is True
    assert _global_status(rows) == "dynamic_separation_grammar_positive_on_1UBQ_1CLL_4AKE_missing_claim_disabled"
