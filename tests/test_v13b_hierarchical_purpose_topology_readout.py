from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
for _path in (REPO_ROOT / "scripts", REPO_ROOT / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from run_v13b_hierarchical_purpose_topology_readout_v0 import (  # noqa: E402
    _classify_pair_topology,
    _evaluate_hierarchical_threshold,
    _select_highest_passing_threshold,
)


class DummyRow:
    sequence = "A" * 160


def test_classifies_intradomain_c_core_without_interdomain_failure() -> None:
    result = _classify_pair_topology(
        (97, 133),
        domain_boundaries_raw="1-75,80-148",
        local_support_max_separation=8,
        long_range_min_separation=24,
    )
    assert result["domain_relation"] == "intradomain_C"
    assert result["topology_role"] == "C_domain_compact_core_candidate"
    assert result["evidence_class"] == "C_domain_core_evidence"


def test_classifies_local_support_separately_from_core() -> None:
    result = _classify_pair_topology(
        (22, 29),
        domain_boundaries_raw="1-75,80-148",
        local_support_max_separation=8,
        long_range_min_separation=24,
    )
    assert result["domain_relation"] == "intradomain_N"
    assert result["topology_role"] == "local_support"
    assert result["evidence_class"] == "local_support"


def test_classifies_interdomain_hinge_candidate() -> None:
    result = _classify_pair_topology(
        (60, 100),
        domain_boundaries_raw="1-75,80-148",
        local_support_max_separation=8,
        long_range_min_separation=24,
    )
    assert result["domain_relation"] == "interdomain"
    assert result["topology_role"] == "interdomain_hinge_candidate"
    assert result["evidence_class"] == "interdomain_hinge_evidence"


def test_hierarchical_threshold_can_select_c_domain_core() -> None:
    pair = (97, 133)
    effective = {"strict": set(), "balanced": {pair}, "balanced_rescue": set(), "monitor": set(), "unknown": set()}
    freqs = {pair: [0.75, 0.72, 0.75, 0.74, 0.66, 0.73, 0.75, 0.72, 0.79, 0.74]}
    row = _evaluate_hierarchical_threshold(
        0.73,
        row=DummyRow(),
        candidate_pairs={pair},
        pair_frequencies=freqs,
        dca_scores={pair: 1.0},
        effective_lane_pairs=effective,
        domain_boundaries_raw="1-75,80-148",
        vote_threshold=7,
        balanced_dca_threshold=0.80,
        local_support_max_separation=8,
        long_range_min_separation=24,
        adaptive_chemical_soft_guard=True,
        legacy_chemical_reference_threshold=0.50,
    )
    assert row["passes_hierarchical_fit"] is True
    assert row["selected_C_domain_core"] == [[97, 133]]
    assert row["selected_interdomain_hinge"] == []
    assert row["support_by_selected_pair"]["97-133"] == 7
    assert 0.0 <= row["mean_frequency_by_selected_pair"]["97-133"] <= 1.0


def test_select_highest_passing_hierarchical_threshold() -> None:
    rows = [
        {"threshold": 0.72, "passes_hierarchical_fit": True, "selected_pair_count": 1},
        {"threshold": 0.73, "passes_hierarchical_fit": True, "selected_pair_count": 1},
        {"threshold": 0.74, "passes_hierarchical_fit": False, "selected_pair_count": 0},
    ]
    selected = _select_highest_passing_threshold(rows)
    assert selected is not None
    assert selected["threshold"] == 0.73
