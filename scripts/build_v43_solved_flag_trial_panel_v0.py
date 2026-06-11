#!/usr/bin/env python3
from __future__ import annotations

"""Build the V43 solved-flag trial panel from the V42 sealed challenge panel."""

import argparse
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
V42_PANEL = REPO_ROOT / "data" / "de_novo_mechanism_language" / "V42" / "panel" / "panel_manifest.json"
DATA_ROOT = REPO_ROOT / "data" / "solved_flag_trial" / "V43"
PANEL_ROOT = DATA_ROOT / "panel"


CONTACT_SCORABLE_CLASSES = {
    "known_anchor",
    "membrane_channel_transporter",
    "soluble_compact_single_domain",
    "multistate_allosteric_metamorphic_binding",
}


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _contact_requirement(target: dict[str, Any]) -> dict[str, Any]:
    group = str(target.get("panel_group"))
    is_contact_scorable = group in CONTACT_SCORABLE_CLASSES and group != "disordered_or_partially_disordered"
    if group == "known_anchor" and target.get("target_name") == "alpha_synuclein_SNCA":
        is_contact_scorable = False
    if group == "weak_or_shallow_evolutionary_information":
        is_contact_scorable = False
    return {
        "contact_scorable": is_contact_scorable,
        "postseal_holdout_types": [
            "structure_or_region_constraint_scoring" if is_contact_scorable else "ensemble_or_function_scoring",
            "operator_region_overlap_scoring",
            "perturbation_effect_scoring",
        ],
        "contact_metric_required_for_candidate": is_contact_scorable,
    }


def build_panel() -> dict[str, Any]:
    v42_panel = _read_json(V42_PANEL, "V42 panel manifest")
    PANEL_ROOT.mkdir(parents=True, exist_ok=True)
    targets = []
    acquisition_entries = []
    for target in v42_panel.get("targets", []):
        enriched = {
            **target,
            "v43_panel_source": str(V42_PANEL),
            "structural_or_ensemble_holdout_status": "postseal_only",
            **_contact_requirement(target),
        }
        targets.append(enriched)
        _write_json(PANEL_ROOT / enriched["target_id"] / "prediction_input_manifest.json", enriched)
        acquisition_entries.append({
            "target_id": enriched["target_id"],
            "target_name": enriched["target_name"],
            "panel_group": enriched["panel_group"],
            "contact_scorable": enriched["contact_scorable"],
            "acquisition_status": "v42_panel_imported_for_v43_solved_flag_trial",
            "replacement_reason": None,
        })
    panel = {
        "kind": "V43_SOLVED_FLAG_TRIAL_PANEL_v0",
        "panel_target_count": len(targets),
        "source_panel": str(V42_PANEL),
        "panel_groups": {
            group: sum(1 for target in targets if target["panel_group"] == group)
            for group in sorted({target["panel_group"] for target in targets})
        },
        "contact_scorable_target_count": sum(1 for target in targets if target["contact_scorable"]),
        "targets": targets,
        "holdouts_created": False,
        "answer_key_available_to_prediction": False,
        "claim_allowed": False,
        "folding_problem_solved": False,
        "protein_folding_solved_candidate": False,
    }
    _write_json(PANEL_ROOT / "panel_manifest.json", panel)
    _write_json(PANEL_ROOT / "acquisition_log.json", {
        "kind": "V43_SOLVED_FLAG_TRIAL_PANEL_ACQUISITION_LOG_v0",
        "panel_target_count": len(targets),
        "replacement_count": 0,
        "entries": acquisition_entries,
    })
    return {
        "kind": "V43_PANEL_BUILD_v0",
        "panel_target_count": len(targets),
        "contact_scorable_target_count": panel["contact_scorable_target_count"],
        "panel_manifest": str(PANEL_ROOT / "panel_manifest.json"),
        "acquisition_log": str(PANEL_ROOT / "acquisition_log.json"),
        "claim_allowed": False,
        "folding_problem_solved": False,
        "protein_folding_solved_candidate": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V43 solved-flag trial panel.")
    parser.parse_args()
    print(json.dumps(build_panel(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
