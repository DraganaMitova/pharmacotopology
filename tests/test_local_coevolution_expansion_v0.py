from pathlib import Path

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_local_coevolution_expansion import build_local_coevolution_contacts, run_local_coevolution_expansion_packet
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

ROOT = Path(__file__).resolve().parents[1]
ROWS = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"


def _fixture():
    rows = load_real_coordinate_visual_rows(ROWS)
    constraints = load_coupling_dataset(COUPLINGS).constraints
    row = next(r for r in rows if r.source_accession == "4AKE:A")
    return rows, constraints, row


def test_local_coevolution_proxy_generates_contacts_without_native_leakage():
    _rows, constraints, row = _fixture()
    contacts, selected, scores = build_local_coevolution_contacts(row, constraints, top_anchor_count=50, window=5, threshold=0.5)
    assert contacts
    assert selected
    assert scores
    assert all(not c.native_truth_used_before_selection for c in contacts)
    assert all(not c.coordinate_truth_used_before_selection for c in contacts)
    assert any(c.channel == "local_coevolution_proxy_window" for c in contacts)


def test_packet_is_bounded_and_marks_proxy_mi_honestly():
    rows, constraints, _row = _fixture()
    packet = run_local_coevolution_expansion_packet(rows, constraints=constraints, evaluation_source_accessions=("4AKE:A",), top_anchor_count=50, window=5, threshold=0.5)
    assert packet.row_count == 1
    assert packet.raw_msa_available_for_true_local_mi is False
    assert packet.local_mi_channel_is_proxy_not_new_msa_calculation is True
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.rows[0].target_native_truth_attached_after_selection_for_evaluation is True


def test_local_coevolution_gate_does_not_auto_claim_solution():
    rows, constraints, _row = _fixture()
    packet = run_local_coevolution_expansion_packet(rows, constraints=constraints, evaluation_source_accessions=("4AKE:A",), top_anchor_count=50, window=5, threshold=0.5)
    assert packet.mean_native_contact_precision_after_audit >= 0.0
    assert packet.mean_long_range_contact_recall_after_audit >= 0.0
    assert packet.folding_problem_solved is False
