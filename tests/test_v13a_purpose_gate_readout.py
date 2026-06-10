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

from run_v13a_purpose_gate_readout_v0 import _evaluate_frequency_threshold, _tail_pair_frequencies


def test_tail_pair_frequencies_count_contact_presence_not_distance(tmp_path) -> None:
    def atom(serial: int, resid: int, x: float, y: float, z: float) -> str:
        return (
            f"ATOM  {serial:5d}  CA  ALA A{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        )

    lines = []
    serial = 1
    # Four models, pair 1-30 contact in exactly two tail frames.
    for model, dist in enumerate([20.0, 20.0, 5.0, 5.0], start=1):
        lines.append(f"MODEL     {model}\n")
        lines.append(atom(serial, 1, 0.0, 0.0, 0.0)); serial += 1
        lines.append(atom(serial, 30, dist, 0.0, 0.0)); serial += 1
        lines.append("ENDMDL\n")
    traj = tmp_path / "trajectory.pdb"
    traj.write_text("".join(lines), encoding="utf-8")

    freqs, _audit, frame_count, tail_count = _tail_pair_frequencies(
        traj,
        contact_cutoff_angstrom=7.0,
        min_separation=24,
        tail_fraction=1.0,
        audit_pairs=set(),
    )
    assert frame_count == 4
    assert tail_count == 4
    assert freqs[(1, 30)] == 0.5
    assert all(0.0 <= value <= 1.0 for value in freqs.values())


def test_adaptive_soft_guard_can_keep_dca_supported_low_chemical_pair() -> None:
    pair = (1, 30)
    payloads = [
        {"frequencies": {pair: 0.8}, "audit_reachability": {}}
        for _ in range(7)
    ]
    effective = {"strict": set(), "balanced": {pair}, "balanced_rescue": set(), "monitor": set(), "unknown": set()}
    # A-C has low chemical score, so hard threshold 0.5 rejects it.
    sequence = "A" + ("G" * 28) + "C" + ("G" * 30)

    hard = _evaluate_frequency_threshold(
        0.75,
        trajectory_payloads=payloads,
        row=object(),  # unused by the readout evaluator
        sequence=sequence,
        dca_scores={pair: 0.95},
        anchor_classes={pair: "balanced"},
        effective_lane_pairs=effective,
        balanced_dca_threshold=0.80,
        chemical_threshold=0.50,
        chemical_policy="hard_threshold",
        legacy_chemical_reference_threshold=0.50,
        vote_threshold=7,
        topology_mode="none",
        domain_boundaries_raw="",
        audit_pairs={pair},
    )
    soft = _evaluate_frequency_threshold(
        0.75,
        trajectory_payloads=payloads,
        row=object(),
        sequence=sequence,
        dca_scores={pair: 0.95},
        anchor_classes={pair: "balanced"},
        effective_lane_pairs=effective,
        balanced_dca_threshold=0.80,
        chemical_threshold=0.50,
        chemical_policy="adaptive_soft_guard",
        legacy_chemical_reference_threshold=0.50,
        vote_threshold=7,
        topology_mode="none",
        domain_boundaries_raw="",
        audit_pairs={pair},
    )

    assert hard["selected_pair_count"] == 0
    assert soft["selected_pair_count"] == 1
    assert soft["selected_balanced_core"] == [[1, 30]]
    assert soft["legacy_chemical_hard_gate_would_block_selected_count"] == 1
    assert soft["mean_frequency_by_selected_pair"]["1-30"] == 0.8
