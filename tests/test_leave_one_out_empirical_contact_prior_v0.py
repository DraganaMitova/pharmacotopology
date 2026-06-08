from pathlib import Path

from pharmacotopology.folding_leave_one_out_empirical_contact_prior import (
    build_leave_one_out_empirical_contact_prior,
)
from pharmacotopology.folding_native_contact_eval import evaluate_contact_prediction
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows


def test_leave_one_out_empirical_prior_excludes_target_and_remains_native_free() -> None:
    rows = load_real_coordinate_visual_rows(Path("data/folding_real_coordinate_visual_8.locked.json"))
    target = next(row for row in rows if row.source_accession == "4AKE:A")
    training = tuple(row for row in rows if row.source_accession != "4AKE:A")

    packet = build_leave_one_out_empirical_contact_prior(
        target_row=target,
        training_rows=training,
        budget_fraction=2.06,
        max_degree=6,
    )
    metric = evaluate_contact_prediction(
        native_pairs=target.native_contact_pairs(),
        predicted_pairs=packet.selected_pairs,
    )

    assert "4AKE:A" not in packet.training_source_accessions
    assert packet.target_native_truth_used_before_selection is False
    assert packet.target_coordinate_truth_used_before_selection is False
    assert packet.alphafold_used_before_selection is False
    assert packet.msa_used_before_selection is False
    assert packet.raw_sequence_persisted is False
    # This learned prior is a real improvement over hand priors, but it is local-heavy.
    assert metric.native_contact_precision > 0.45
    assert metric.native_contact_recall > 0.35
    assert metric.long_range_contact_recall == 0.0
    assert not (metric.native_contact_precision >= 0.70 and metric.native_contact_recall >= 0.70)


def test_leave_one_out_prior_rejects_target_training_leakage() -> None:
    rows = load_real_coordinate_visual_rows(Path("data/folding_real_coordinate_visual_8.locked.json"))
    target = next(row for row in rows if row.source_accession == "4AKE:A")
    try:
        build_leave_one_out_empirical_contact_prior(target_row=target, training_rows=(target,))
    except ValueError as exc:
        assert "target" in str(exc)
    else:
        raise AssertionError("target leakage should be rejected")
