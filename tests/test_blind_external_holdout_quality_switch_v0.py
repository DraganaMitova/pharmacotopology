from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.run_blind_external_holdout_battery_v0 import (  # noqa: E402
    SELF_CRITICAL_QUALITY_SWITCH_SELECTOR,
    _apply_self_critical_quality_switch,
)
from scripts.score_cached_contact_map_modes_v0 import (  # noqa: E402
    ANCHORED_SEQUENCE_COUPLING_BALANCE_SELECTOR,
    COHERENT_RELAXED_SEQUENCE_COUPLING_BALANCE_SELECTOR,
    CONFLICT_SHELL_DENSITY_RESCUE_SELECTOR,
    GLOBAL_CONTACT_MAP_COLLAPSE_SELECTOR,
    RELAXED_ANCHORED_SEQUENCE_COUPLING_BALANCE_SELECTOR,
    ROW_LOCAL_CRITICAL_MODE_SELECTOR,
    _sequence_coupling_expansion_decision,
    _sequence_feature_fallback_decision,
    _select_row_local_critical_mode,
)


def _row(
    *,
    target_id: str,
    quality: float,
    phase_f1: float,
    density_f1: float,
    selected_event_count: int = 0,
    phase_coverage_mode: str = "diagonal_half_radius",
    spine_f1: float = 0.0,
    conflict_f1: float = 0.0,
    shell_f1: float = 0.0,
    phase_conflict_f1: float = 0.0,
) -> dict[str, object]:
    return {
        "target_id": target_id,
        "selected_event_count": selected_event_count,
        "coupling_quality_mean_effseq_per_length": quality,
        "phase_coverage_scaffold_contact_count": 10,
        "phase_coverage_scaffold_exact_contact_precision": 0.25,
        "phase_coverage_scaffold_exact_long_range_contact_recall": phase_f1,
        "phase_coverage_scaffold_precision_recall_f1": phase_f1,
        "phase_coverage_scaffold_contact_map_perfect": False,
        "phase_coverage_scaffold_mode": phase_coverage_mode,
        "region_density_top_l_scaffold_contact_count": 20,
        "region_density_top_l_scaffold_exact_contact_precision": 0.5,
        "region_density_top_l_scaffold_exact_long_range_contact_recall": (
            density_f1
        ),
        "region_density_top_l_scaffold_precision_recall_f1": density_f1,
        "region_density_top_l_scaffold_contact_map_perfect": False,
        "phase_density_spine_scaffold_contact_count": 15,
        "phase_density_spine_scaffold_exact_contact_precision": 0.4,
        "phase_density_spine_scaffold_exact_long_range_contact_recall": (
            spine_f1
        ),
        "phase_density_spine_scaffold_precision_recall_f1": spine_f1,
        "phase_density_spine_scaffold_contact_map_perfect": False,
        "phase_density_conflict_consensus_scaffold_contact_count": 12,
        "phase_density_conflict_consensus_scaffold_exact_contact_precision": 0.45,
        "phase_density_conflict_consensus_scaffold_exact_long_range_contact_recall": (
            conflict_f1
        ),
        "phase_density_conflict_consensus_scaffold_precision_recall_f1": (
            conflict_f1
        ),
        "phase_density_conflict_consensus_scaffold_contact_map_perfect": False,
        "phase_density_conflict_shell_scaffold_contact_count": 18,
        "phase_density_conflict_shell_scaffold_exact_contact_precision": 0.35,
        "phase_density_conflict_shell_scaffold_exact_long_range_contact_recall": (
            shell_f1
        ),
        "phase_density_conflict_shell_scaffold_precision_recall_f1": shell_f1,
        "phase_density_conflict_shell_scaffold_contact_map_perfect": False,
        "phase_density_conflict_phase_confidence_scaffold_contact_count": 21,
        "phase_density_conflict_phase_confidence_scaffold_exact_contact_precision": 0.34,
        "phase_density_conflict_phase_confidence_scaffold_exact_long_range_contact_recall": (
            phase_conflict_f1
        ),
        "phase_density_conflict_phase_confidence_scaffold_precision_recall_f1": (
            phase_conflict_f1
        ),
        "phase_density_conflict_phase_confidence_scaffold_contact_map_perfect": (
            False
        ),
    }


def test_self_critical_quality_switch_uses_largest_non_native_quality_gap() -> None:
    rows = [
        _row(target_id="low_a", quality=5.2, phase_f1=0.24, density_f1=0.13),
        _row(target_id="low_b", quality=5.3, phase_f1=0.23, density_f1=0.18),
        _row(target_id="high_a", quality=12.7, phase_f1=0.09, density_f1=0.20),
        _row(target_id="high_b", quality=13.2, phase_f1=0.13, density_f1=0.35),
    ]

    audit = _apply_self_critical_quality_switch(rows)

    assert audit["selector_name"] == SELF_CRITICAL_QUALITY_SWITCH_SELECTOR
    assert audit["boundary_available"] is True
    assert audit["boundary_value"] == 9.0
    assert audit["selected_modes"] == (
        "phase_coverage",
        "region_density_top_l",
    )
    assert [row["self_critical_quality_switch_mode"] for row in rows] == [
        "phase_coverage",
        "phase_coverage",
        "region_density_top_l",
        "region_density_top_l",
    ]
    assert [
        row["self_critical_quality_switch_precision_recall_f1"]
        for row in rows
    ] == [0.24, 0.23, 0.20, 0.35]


def test_self_critical_quality_switch_stays_conservative_without_batch_boundary() -> None:
    rows = [
        _row(target_id="single", quality=13.2, phase_f1=0.13, density_f1=0.35)
    ]

    audit = _apply_self_critical_quality_switch(rows)

    assert audit["boundary_available"] is False
    assert rows[0]["self_critical_quality_switch_mode"] == "phase_coverage"
    assert rows[0]["self_critical_quality_switch_precision_recall_f1"] == 0.13


def test_self_critical_quality_switch_generates_spine_for_high_event_conflict() -> None:
    rows = [
        _row(
            target_id="low_quality",
            quality=5.2,
            selected_event_count=4,
            phase_f1=0.24,
            density_f1=0.13,
            spine_f1=0.21,
            conflict_f1=0.18,
            shell_f1=0.16,
        ),
        _row(
            target_id="high_quality_compact",
            quality=13.2,
            selected_event_count=3,
            phase_f1=0.13,
            density_f1=0.35,
            spine_f1=0.33,
            conflict_f1=0.31,
            shell_f1=0.29,
        ),
        _row(
            target_id="high_quality_conflict",
            quality=12.7,
            selected_event_count=30,
            phase_f1=0.09,
            density_f1=0.20,
            spine_f1=0.22,
            conflict_f1=0.24,
            shell_f1=0.26,
        ),
    ]

    audit = _apply_self_critical_quality_switch(rows)

    assert audit["event_boundary_available"] is True
    assert rows[0]["self_critical_quality_switch_mode"] == "phase_coverage"
    assert rows[1]["self_critical_quality_switch_mode"] == (
        "region_density_top_l"
    )
    assert rows[2]["self_critical_quality_switch_mode"] == (
        "phase_density_conflict_shell"
    )
    assert rows[2]["self_critical_quality_switch_precision_recall_f1"] == 0.26


def test_self_critical_quality_switch_uses_phase_confidence_for_square_high_event() -> None:
    rows = [
        _row(
            target_id="low_quality",
            quality=5.2,
            selected_event_count=4,
            phase_f1=0.24,
            density_f1=0.13,
            shell_f1=0.19,
            phase_conflict_f1=0.18,
        ),
        _row(
            target_id="high_quality_compact",
            quality=13.2,
            selected_event_count=3,
            phase_f1=0.13,
            density_f1=0.35,
            shell_f1=0.29,
            phase_conflict_f1=0.28,
        ),
        _row(
            target_id="high_quality_square_conflict",
            quality=12.7,
            selected_event_count=30,
            phase_coverage_mode="square_one_cell",
            phase_f1=0.09,
            density_f1=0.20,
            shell_f1=0.26,
            phase_conflict_f1=0.27,
        ),
    ]

    audit = _apply_self_critical_quality_switch(rows)

    assert audit["high_event_conflict_feature"] == "phase_coverage_scaffold_mode"
    assert rows[2]["self_critical_quality_switch_mode"] == (
        "phase_density_conflict_phase_confidence"
    )
    assert rows[2]["self_critical_quality_switch_precision_recall_f1"] == 0.27


def test_self_critical_quality_switch_uses_shell_for_low_quality_square_phase() -> None:
    rows = [
        _row(
            target_id="low_quality_diagonal",
            quality=5.2,
            selected_event_count=4,
            phase_coverage_mode="diagonal_half_radius",
            phase_f1=0.24,
            density_f1=0.13,
            shell_f1=0.19,
        ),
        _row(
            target_id="low_quality_square",
            quality=5.3,
            selected_event_count=8,
            phase_coverage_mode="square_one_cell",
            phase_f1=0.23,
            density_f1=0.18,
            shell_f1=0.25,
        ),
        _row(
            target_id="high_quality_compact",
            quality=13.2,
            selected_event_count=3,
            phase_coverage_mode="diagonal_half_radius",
            phase_f1=0.13,
            density_f1=0.35,
            shell_f1=0.29,
        ),
    ]

    audit = _apply_self_critical_quality_switch(rows)

    assert audit["low_quality_conflict_feature"] == "phase_coverage_scaffold_mode"
    assert rows[0]["self_critical_quality_switch_mode"] == "phase_coverage"
    assert rows[1]["self_critical_quality_switch_mode"] == (
        "phase_density_conflict_shell"
    )
    assert rows[1]["self_critical_quality_switch_precision_recall_f1"] == 0.25


def test_cached_row_local_critical_switch_routes_by_internal_phase_conflict() -> None:
    assert ROW_LOCAL_CRITICAL_MODE_SELECTOR == (
        "cached_row_local_critical_mode_switch_v0"
    )
    assert (
        _select_row_local_critical_mode(
            {
                "phase_mode": "point",
                "phase_radius": 1,
                "selected_event_count": 1,
                "direct_density_ratio": 2.5,
            }
        )
        == "phase_density_conflict_shell"
    )
    assert (
        _select_row_local_critical_mode(
            {
                "phase_mode": "square",
                "phase_radius": 1,
                "selected_event_count": 4,
                "direct_density_ratio": 2.5,
            }
        )
        == "region_density_top_l"
    )
    assert (
        _select_row_local_critical_mode(
            {
                "phase_mode": "square",
                "phase_radius": 2,
                "selected_event_count": 8,
                "direct_density_ratio": 4.8,
            }
        )
        == "phase_density_conflict_shell"
    )
    assert (
        _select_row_local_critical_mode(
            {
                "phase_mode": "diagonal",
                "phase_radius": 10,
                "selected_event_count": 5,
                "direct_density_ratio": 7.1,
            }
        )
        == "scaffold"
    )
    assert (
        _select_row_local_critical_mode(
            {
                "phase_mode": "diagonal",
                "phase_radius": 7,
                "selected_event_count": 3,
                "direct_density_ratio": 5.2,
            }
        )
        == "region_density_boundary"
    )


def test_anchored_sequence_coupling_balance_uses_compact_self_critical_gate() -> None:
    assert ANCHORED_SEQUENCE_COUPLING_BALANCE_SELECTOR == (
        "cached_anchored_sequence_coupling_balance_v0"
    )
    assert RELAXED_ANCHORED_SEQUENCE_COUPLING_BALANCE_SELECTOR == (
        "cached_relaxed_anchored_sequence_coupling_balance_v0"
    )
    assert COHERENT_RELAXED_SEQUENCE_COUPLING_BALANCE_SELECTOR == (
        "cached_coherent_relaxed_sequence_coupling_balance_v0"
    )
    assert CONFLICT_SHELL_DENSITY_RESCUE_SELECTOR == (
        "cached_conflict_shell_density_rescue_v0"
    )
    assert GLOBAL_CONTACT_MAP_COLLAPSE_SELECTOR == (
        "cached_global_contact_map_collapse_v0"
    )
    assert _sequence_coupling_expansion_decision(
        phase_mode="square",
        phase_radius=2,
        selected_event_count=8,
        anchor_count=342,
        addition_count=4,
    ) == (False, "square_phase_anchor_only", 21)
    assert _sequence_coupling_expansion_decision(
        phase_mode="point",
        phase_radius=1,
        selected_event_count=1,
        anchor_count=41,
        addition_count=6,
    ) == (False, "point_single_event_anchor_only", 8)
    assert _sequence_coupling_expansion_decision(
        phase_mode="point",
        phase_radius=1,
        selected_event_count=2,
        anchor_count=92,
        addition_count=45,
    ) == (True, "point_multi_event_sequence_balance", 11)
    assert _sequence_coupling_expansion_decision(
        phase_mode="diagonal",
        phase_radius=7,
        selected_event_count=3,
        anchor_count=16,
        addition_count=4,
    ) == (True, "diagonal_compact_sequence_balance", 11)
    assert _sequence_coupling_expansion_decision(
        phase_mode="diagonal",
        phase_radius=6,
        selected_event_count=3,
        anchor_count=76,
        addition_count=56,
    ) == (False, "diagonal_diffuse_or_singleton_anchor_only", 15)

    assert _sequence_feature_fallback_decision(
        anchor_mode="region_density_top_l",
        phase_mode="diagonal",
        compact_diagonal_limit=15,
        addition_count=2,
    ) == (True, "feature_compact_density_anchor_balance")
    assert _sequence_feature_fallback_decision(
        anchor_mode="phase_coverage",
        phase_mode="diagonal",
        compact_diagonal_limit=19,
        addition_count=1,
    ) == (False, "feature_fallback_requires_density_anchor")
    assert _sequence_feature_fallback_decision(
        anchor_mode="region_density_top_l",
        phase_mode="square",
        compact_diagonal_limit=21,
        addition_count=4,
    ) == (False, "feature_fallback_square_phase_anchor_only")
