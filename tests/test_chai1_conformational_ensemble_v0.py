from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_chai1_conformational_ensemble import build_conformational_ensemble_packet


def _write_ca_pdb(path: Path, *, d118_k136: float) -> None:
    lines = []
    serial = 1
    for i in range(1, 215):
        x = float(i) * 3.8
        y = 0.0
        z = 0.0
        if i == 118:
            x, y, z = 0.0, 0.0, 0.0
        if i == 136:
            x, y, z = d118_k136, 0.0, 0.0
        lines.append(
            f"ATOM  {serial:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 80.00           C\n"
        )
        serial += 1
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def test_conformational_ensemble_selects_open_proxy_short_d118_k136(tmp_path):
    closed = tmp_path / "closed_like.pdb"
    open_like = tmp_path / "open_like.pdb"
    _write_ca_pdb(closed, d118_k136=20.0)
    _write_ca_pdb(open_like, d118_k136=5.0)

    packet = build_conformational_ensemble_packet(
        structure_dir=tmp_path,
        source_id="chai1_single_sequence_conformational_ensemble",
        chain_id="A",
        minimum_ca_count=200,
    )

    assert packet.parsed_conformer_count == 2
    assert packet.selected_conformer_path == str(open_like)
    assert packet.native_truth_used_before_selection is False
    assert packet.msa_used_by_this_module is False
    assert packet.raw_sequence_exposed is False


def test_conformational_ensemble_rejects_alphafold_like_input(tmp_path):
    bad = tmp_path / "AF-P69441-F1-model_v4.pdb"
    _write_ca_pdb(bad, d118_k136=5.0)

    packet = build_conformational_ensemble_packet(
        structure_dir=tmp_path,
        source_id="chai1_but_path_contains_alphafold",
        chain_id="A",
        minimum_ca_count=200,
    )

    assert packet.alphafold_like_input_rejected is True
    assert packet.parsed_conformer_count == 0
    assert packet.contacts == ()
