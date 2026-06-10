from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
for _path in (REPO_ROOT / "scripts", REPO_ROOT / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from run_v15_4ake_dynamic_grammar_bridge_v0 import build_bridge  # noqa: E402
from run_v15_dynamic_separation_grammar_readout_v0 import (  # noqa: E402
    _global_status,
    normalize_1cll_dynamic,
    normalize_1ubq_dynamic,
    normalize_4ake_dynamic,
)


def test_4ake_bridge_does_not_synthesize_visual_only_positive(tmp_path: Path) -> None:
    visual = tmp_path / "visual"
    visual.mkdir()
    (visual / "contact_map_overlay.svg").write_text("<svg/>")
    audit = tmp_path / "audit_4ake_v13c.json"
    audit.write_text("[]")
    bridge = build_bridge(
        role_cert_paths=[],
        legacy_visual_dir=visual,
        legacy_cert_paths=[],
        input_material_paths=[audit],
    )
    row = bridge["protein_row"]
    assert row["positive_evidence_found"] is False
    assert row["claim_allowed"] is False
    assert "bridge_pending" in row["claim_lock_status"]
    assert row["fixed_residue_cutoff_used"] is False


def test_4ake_bridge_can_normalize_machine_readable_role_artifact(tmp_path: Path) -> None:
    role_cert = tmp_path / "role.json"
    role_cert.write_text(
        '{"selected_balanced_core": [[10, 30]], "noise_added": 0, "long_range_evidence_polluted": false}'
    )
    bridge = build_bridge(
        role_cert_paths=[role_cert],
        legacy_visual_dir=tmp_path / "none",
        legacy_cert_paths=[],
        input_material_paths=[],
    )
    row = bridge["protein_row"]
    assert row["positive_evidence_found"] is True
    assert row["selected_pairs"] == ["10-30"]
    assert row["claim_allowed"] is False
    assert row["dynamic_pair_roles"]["10-30"]["fixed_residue_cutoff_used"] is False


def test_v15_global_status_tracks_4ake_bridge_pending() -> None:
    rows = [
        normalize_4ake_dynamic({"protein_row": {"artifact_status": "present_legacy_visual_and_input_material_bridge_pending_machine_readable_role_artifact", "positive_evidence_found": False, "claim_lock_status": "bridge_pending_machine_readable_role_artifact"}}),
        normalize_1ubq_dynamic({"selected_frequency_band": {"selected_pair_count": 1, "selected_balanced_core": [[23, 48]]}}),
        normalize_1cll_dynamic({"selected_frequency_band": {"selected_pair_count": 1, "selected_C_domain_core": [[97, 133]]}}),
    ]
    assert _global_status(rows) == "dynamic_separation_grammar_positive_on_1UBQ_1CLL_4AKE_bridge_pending_claim_disabled"
