from __future__ import annotations

from pathlib import Path

from src.pharmacotopology.folding_coupling_nucleus_selector import (
    COUPLING_ADAPTIVE_LOW_SIGNAL_FUTURE_MIN,
    build_coupling_nucleus_context,
    _adaptive_coupling_floor_profile,
    _passes_coupling_future_gate,
    _passes_low_signal_adaptive_rescue,
)
from src.pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from src.pharmacotopology.folding_real_coordinate_visual_benchmark import (
    load_real_coordinate_visual_rows,
)

ROOT = Path(__file__).resolve().parents[1]


def _context():
    rows = load_real_coordinate_visual_rows(
        ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
    )
    dataset = load_coupling_dataset(
        ROOT
        / "data"
        / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
    )
    return build_coupling_nucleus_context(rows=rows, coupling_dataset=dataset)


def _first_event(context, accession: str):
    row = next(row for row in context.rows if row.source_accession == accession)
    return next(event for event in context.competitive_events if event.row_id == row.row_id)


def test_adaptive_floor_uses_phase_complexity_and_msa_quality() -> None:
    context = _context()

    one_mbn = _adaptive_coupling_floor_profile(_first_event(context, "1MBN:A"), context)
    four_ake = _adaptive_coupling_floor_profile(_first_event(context, "4AKE:A"), context)
    one_tim = _adaptive_coupling_floor_profile(_first_event(context, "1TIM:A"), context)

    assert one_mbn.phase_mode in {"square", "diagonal"}
    assert one_mbn.sequence_complexity > 0.85
    assert one_mbn.coupling_depth_over_length > 8.0
    assert one_mbn.future_preservation_floor <= 0.36
    assert one_mbn.low_signal_rescue_enabled is True

    assert four_ake.phase_mode == "diagonal"
    assert four_ake.coupling_depth_over_length > 8.0
    assert four_ake.future_preservation_floor == 0.35
    assert four_ake.low_signal_rescue_enabled is True

    # 1TIM is the safety counterexample: good MSA, but its row-level future
    # signal ceiling stays below the low-signal rescue floor, so it is guarded.
    assert one_tim.row_future_preservation_ceiling < COUPLING_ADAPTIVE_LOW_SIGNAL_FUTURE_MIN
    assert one_tim.future_preservation_floor >= 0.42
    assert one_tim.adaptive_gate_enabled is False
    assert one_tim.low_signal_rescue_enabled is False


def test_low_signal_rescue_attacks_4ake_1mbn_without_opening_1tim() -> None:
    context = _context()
    counts: dict[str, tuple[int, int]] = {}
    for accession in ("1MBN:A", "4AKE:A", "1TIM:A"):
        row = next(row for row in context.rows if row.source_accession == accession)
        profile = _adaptive_coupling_floor_profile(
            next(event for event in context.competitive_events if event.row_id == row.row_id),
            context,
        )
        passed = [
            event
            for event in context.competitive_events
            if event.row_id == row.row_id and _passes_coupling_future_gate(event, context)
        ]
        rescued = [
            event
            for event in context.competitive_events
            if event.row_id == row.row_id
            and _passes_low_signal_adaptive_rescue(event, context, profile)
        ]
        counts[accession] = (len(passed), len(rescued))

    assert counts["1MBN:A"][0] > 0
    assert counts["1MBN:A"][1] > 0
    assert counts["4AKE:A"][0] > 0
    assert counts["4AKE:A"][1] > 0
    assert counts["1TIM:A"] == (0, 0)
