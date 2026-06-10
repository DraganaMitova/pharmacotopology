from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
for _path in (REPO_ROOT / "scripts", REPO_ROOT / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from run_v15_4ake_balanced_candidate_readout_v0 import (  # noqa: E402
    _domain_relation,
    _parse_grid,
    _sequence_context,
)
from run_v15_4ake_dynamic_grammar_bridge_v0 import build_bridge  # noqa: E402


def test_parse_frequency_grid_descending() -> None:
    grid = _parse_grid("0.98:0.50:0.01")
    assert grid[0] == 0.98
    assert grid[-1] == 0.5
    assert 0.76 in grid


def test_domain_relation_without_fixed_cutoff_for_4ake_domains() -> None:
    boundaries = ((1, 29), (30, 67), (68, 117), (118, 160), (161, 214))
    relation, role, evidence = _domain_relation((100, 150), boundaries)
    assert relation == "interdomain_D3_D4"
    assert role == "domain_hinge_or_interdomain_closure_candidate"
    assert evidence == "domain_hinge_or_interdomain_evidence"
    ctx = _sequence_context((100, 150), 214, relation)
    assert ctx["fixed_residue_cutoff_used"] is False
    assert ctx["separation_filter_applied"] is False


def test_4ake_balanced_candidate_cert_can_be_bridged(tmp_path: Path) -> None:
    role_cert = tmp_path / "v15_4ake_balanced_candidate_readout_certificate.json"
    role_cert.write_text(json.dumps({
        "selected_balanced_core": [[100, 150]],
        "selected_hinge_or_interdomain": ["100-150"],
        "support_by_selected_pair": {"100-150": 8},
        "noise_added": 0,
        "long_range_evidence_polluted": False,
        "classification_coverage_ratio": 1.0,
        "chemical_policy": "adaptive_soft_guard_report_only_not_hard_kill",
        "topology_policy": "hierarchical_domain_core_plus_interdomain_hinge_context",
    }))
    bridge = build_bridge(
        role_cert_paths=[role_cert],
        legacy_visual_dir=tmp_path / "none",
        legacy_cert_paths=[],
        input_material_paths=[],
    )
    row = bridge["protein_row"]
    assert row["positive_evidence_found"] is True
    assert "100-150" in row["selected_pairs"]
    assert row["fixed_residue_cutoff_used"] is False
    assert row["claim_allowed"] is False
