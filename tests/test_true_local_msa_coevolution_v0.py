from pathlib import Path

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows
from pharmacotopology.folding_true_local_msa_coevolution import build_true_local_msa_contacts, parse_msa, run_true_local_msa_packet

ROOT = Path(__file__).resolve().parents[1]
ROWS = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"


def _fixture():
    rows = load_real_coordinate_visual_rows(ROWS)
    constraints = load_coupling_dataset(COUPLINGS).constraints
    row = next(r for r in rows if r.source_accession == "4AKE:A")
    safe = [c for c in constraints if (c.row_id == row.row_id or c.source_accession == row.source_accession) and c.i < c.j and c.sequence_separation >= 30 and not c.coordinate_truth_used_to_build_constraint and not c.native_truth_used_before_coupling_selection and not c.structure_model_used]
    safe.sort(key=lambda c: (-(c.apc_corrected_score or c.confidence), c.i, c.j))
    anchor = next(c for c in safe if c.i + 1 < c.j + 1 and (c.j + 1) - (c.i + 1) >= 30 and c.i + 1 <= row.sequence_length and c.j + 1 <= row.sequence_length)
    return rows, constraints, row, anchor


def _write_synthetic_msa(tmp_path, row, anchor):
    i2 = anchor.i + 1
    j2 = anchor.j + 1
    records = [("query", row.sequence)]
    for idx in range(40):
        chars = list(row.sequence)
        # Keep anchor residues identical to the query so these sequences survive
        # the anchor-conditioned sub-MSA filter.
        chars[anchor.i - 1] = row.sequence[anchor.i - 1]
        chars[anchor.j - 1] = row.sequence[anchor.j - 1]
        if idx % 2 == 0:
            chars[i2 - 1] = "A"
            chars[j2 - 1] = "A"
        else:
            chars[i2 - 1] = "C"
            chars[j2 - 1] = "C"
        records.append(("synthetic_%02d" % idx, "".join(chars)))
    path = tmp_path / "4ake.synthetic_msa.afa"
    path.write_text("\n".join(">%s\n%s" % (name, seq) for name, seq in records), encoding="utf-8")
    return path


def test_parse_msa_maps_query_positions(tmp_path):
    _rows, _constraints, row, anchor = _fixture()
    path = _write_synthetic_msa(tmp_path, row, anchor)
    msa = parse_msa(path, row.sequence)
    assert msa.sequence_count == 41
    assert msa.query_to_alignment[1] == 0
    assert msa.query_sequence == row.sequence


def test_true_local_msa_channel_uses_real_mi_not_proxy(tmp_path):
    _rows, constraints, row, anchor = _fixture()
    path = _write_synthetic_msa(tmp_path, row, anchor)
    msa = parse_msa(path, row.sequence)
    contacts, selected, scores, resolved_threshold, candidate_local_pair_count = build_true_local_msa_contacts(
        row,
        constraints,
        msa,
        top_anchor_count=50,
        window=5,
        threshold=0.25,
        min_filtered_sequences=16,
    )
    assert selected
    assert scores
    assert resolved_threshold >= 0.0
    assert candidate_local_pair_count >= 0
    assert any(c.channel == "true_local_msa_mi_window" for c in contacts)
    assert all(c.true_msa_mi_used for c in contacts)
    assert all(not c.native_truth_used_before_selection for c in contacts)
    assert all(not c.coordinate_truth_used_before_selection for c in contacts)


def test_true_local_msa_packet_marks_raw_msa_available(tmp_path):
    rows, constraints, row, anchor = _fixture()
    path = _write_synthetic_msa(tmp_path, row, anchor)
    packet = run_true_local_msa_packet(rows, constraints, path, evaluation_source_accessions=("4AKE:A",), threshold=0.25, min_filtered_sequences=16)
    assert packet.row_count == 1
    assert packet.raw_msa_available_for_true_local_mi is True
    assert packet.rows[0].raw_msa_available_for_true_local_mi is True
    assert packet.universal_physical_law_claim_allowed is False
