from pathlib import Path

from pharmacotopology.folding_independent_contact_evidence import (
    IndependentContactEvidencePair,
    candidate_region_evidence_from_events,
    evaluate_ensemble_contacts,
    native_coordinate_positive_control_evidence,
    parse_ca_pdb_points,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows


BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")


def _four_ake_row():
    rows = load_real_coordinate_visual_rows(BENCHMARK_FILE)
    return next(row for row in rows if row.source_accession == "4AKE:A")


def _manual_candidate_evidence(row, pair):
    return IndependentContactEvidencePair(
        row_id=row.row_id,
        source_accession=row.source_accession,
        source_id="candidate_region_pool",
        source_kind="candidate_region_sequence_closure_source_v0",
        source_family="candidate_region",
        i=pair[0],
        j=pair[1],
        confidence=0.5,
        sequence_separation=pair[1] - pair[0],
    )


def test_parse_ca_pdb_points_reads_chain_and_bfactor_confidence(tmp_path):
    pdb = tmp_path / "model.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1      11.111  12.222  13.333  1.00 95.00           C\n"
        "ATOM      2  CA  GLY A   2      14.111  15.222  16.333  1.00 80.00           C\n"
        "ATOM      3  CA  SER B   1      99.111  99.222  99.333  1.00 10.00           C\n",
        encoding="utf-8",
    )

    points = parse_ca_pdb_points(pdb, chain_id="A")

    assert len(points) == 2
    assert points[0].sequence_index == 1
    assert points[0].residue_number == 1
    assert points[0].confidence == 0.95
    assert points[1].confidence == 0.8


def test_ensemble_requires_independent_structure_source_by_default():
    row = _four_ake_row()
    native_long_pair = next(pair for pair in row.native_contact_pairs() if pair[1] - pair[0] >= 24)
    candidate = _manual_candidate_evidence(row, native_long_pair)

    report, decisions, _metric = evaluate_ensemble_contacts(
        row=row,
        evidence_sources=(candidate,),
        min_votes_required=2,
    )

    assert report.final_pair_count == 0
    assert report.benchmark_claim_allowed is False
    assert report.claim_rejection_reason == "missing_independent_structure_source"
    assert decisions[0].selection_reason == "missing_independent_structure_support"


def test_native_coordinate_positive_control_is_selected_but_rejected_as_leakage():
    row = _four_ake_row()
    native_long_pair = next(pair for pair in row.native_contact_pairs() if pair[1] - pair[0] >= 24)
    candidate = _manual_candidate_evidence(row, native_long_pair)
    native_positive_control = native_coordinate_positive_control_evidence(row)

    report, decisions, _metric = evaluate_ensemble_contacts(
        row=row,
        evidence_sources=(candidate, *native_positive_control),
        min_votes_required=2,
    )

    selected = [decision for decision in decisions if decision.selected]
    assert native_long_pair in {decision.pair() for decision in selected}
    assert report.final_pair_count == 1
    assert report.contact_precision == 1.0
    assert report.long_range_precision == 1.0
    assert report.benchmark_claim_allowed is False
    assert report.claim_rejection_reason == "independent_source_is_native_coordinate_leakage_positive_control"
    assert report.coordinate_truth_used_before_selection is True
    assert report.native_truth_used_before_selection is True
