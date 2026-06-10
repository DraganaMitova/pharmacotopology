from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = REPO_ROOT / "scripts"
SRC_ROOT = REPO_ROOT / "src"
for _path in (SCRIPTS_ROOT, SRC_ROOT):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from run_v13a_purpose_gate_readout_v0 import (
    _build_v13_purpose_preflight,
    _classify_target_purpose,
    _select_highest_passing_threshold,
)


def test_single_domain_compact_does_not_require_strict_or_rescue() -> None:
    purpose = _classify_target_purpose(
        sequence_length=76,
        topology_mode="none",
        domain_boundaries_raw="",
        effective_lane_counts={"strict": 0, "balanced": 7, "balanced_rescue": 0, "monitor": 70},
    )
    assert purpose["target_role"] == "single_domain_compact"
    assert purpose["strict_required"] is False
    assert purpose["rescue_required"] is False

    preflight = _build_v13_purpose_preflight(
        input_preflight={
            "status": "ready",
            "checks": {
                "external_coupling_file_exists": True,
                "anchor_profile_file_exists": True,
            },
            "target_pdb_provenance": {
                "sequence_exact_match": True,
                "prepared_target_coverage_ratio": 1.0,
            },
        },
        trajectory_count=10,
        audit_pair_count=1,
        effective_lane_counts={"strict": 0, "balanced": 7, "balanced_rescue": 0, "monitor": 70},
        target_purpose=purpose,
    )
    assert preflight["status"] == "ready"
    assert preflight["strict_absence_is_fatal"] is False
    assert "effective_strict_count_positive" not in preflight["failed_checks"]


def test_low_signal_abstains_without_balanced_lane() -> None:
    purpose = _classify_target_purpose(
        sequence_length=180,
        topology_mode="none",
        domain_boundaries_raw="",
        effective_lane_counts={"strict": 0, "balanced": 0, "balanced_rescue": 0, "monitor": 0},
    )
    assert purpose["target_role"] == "low_signal_abstain"


def test_selects_highest_passing_threshold_not_hardcoded() -> None:
    rows = [
        {"threshold": 0.74, "passes_purpose_fit": True, "selected_pair_count": 2},
        {"threshold": 0.76, "passes_purpose_fit": True, "selected_pair_count": 1},
        {"threshold": 0.77, "passes_purpose_fit": False, "selected_pair_count": 0},
    ]
    selected = _select_highest_passing_threshold(rows)
    assert selected is not None
    assert selected["threshold"] == 0.76
