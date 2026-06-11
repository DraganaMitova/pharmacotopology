#!/usr/bin/env python3
from __future__ import annotations

"""Build V38 blind mechanism generalization panel.

The panel contains three known positives from V36/V37 plus six near-decoys.
Decision-facing masked dossiers preserve mechanism evidence but remove target
names. The answer key is stored separately and is not needed for assignment.
"""

import argparse
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
V36_ROOT = REPO_ROOT / "data" / "external_evidence_dossiers"
PANEL_ROOT = REPO_ROOT / "data" / "blind_mechanism_panels" / "V38"


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object: {path}")
    return data


def _source(source_id: str, source_type: str, source_name: str, citation: str, statement: str) -> dict[str, Any]:
    return {
        "source_id": source_id,
        "source_type": source_type,
        "source_name": source_name,
        "source_url_or_citation": citation,
        "source_date_or_version": "accessed_2026-06-11",
        "evidence_statement": statement,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "claim_allowed": False,
    }


def _dossier(
    panel_id: str,
    target: str,
    group: str,
    expected_class: str,
    evidence_buckets: list[str],
    mechanism_evidence: list[str],
    grammar_rules: dict[str, Any],
    source_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "kind": "V38_PANEL_DOSSIER_v0",
        "panel_id": panel_id,
        "target": target,
        "group": group,
        "expected_mechanism_class_for_answer_key_only": expected_class,
        "evidence_buckets": evidence_buckets,
        "mechanism_evidence": mechanism_evidence,
        "grammar_rules": grammar_rules,
        "source_rows": source_rows,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
    }


def _mask_text(text: str, target: str) -> str:
    replacements = {
        "KcsA": "the masked protein",
        "XCL1": "the masked protein",
        "lymphotactin": "the masked protein",
        "alpha-synuclein": "the masked protein",
        "Alpha-synuclein": "The masked protein",
        "SNCA": "the masked protein",
        "bacteriorhodopsin": "the masked protein",
        "Bacteriorhodopsin": "The masked protein",
        "AQP1": "the masked protein",
        "aquaporin-1": "the masked protein",
        "CXCL8": "the masked protein",
        "interleukin-8": "the masked protein",
        "ubiquitin": "the masked protein",
        "lysozyme": "the masked protein",
        "myoglobin": "the masked protein",
    }
    out = text.replace(target, "the masked protein")
    for src, dst in replacements.items():
        out = out.replace(src, dst)
    return out


def _masked_dossier(dossier: dict[str, Any]) -> dict[str, Any]:
    target = str(dossier["target"])
    return {
        "kind": "V38_MASKED_DOSSIER_v0",
        "masked_id": dossier["panel_id"],
        "evidence_buckets": list(dossier["evidence_buckets"]),
        "mechanism_evidence": [_mask_text(str(item), target) for item in dossier["mechanism_evidence"]],
        "grammar_rules": dict(dossier["grammar_rules"]),
        "source_provenance": [
            {
                "source_id": row["source_id"],
                "source_type": row["source_type"],
                "source_url_or_citation": row["source_url_or_citation"],
                "coordinate_derived": row["coordinate_derived"],
                "internal_runtime_source": row["internal_runtime_source"],
                "native_metrics_used_for_selection": row["native_metrics_used_for_selection"],
                "coordinate_truth_used_before_selection": row["coordinate_truth_used_before_selection"],
                "claim_allowed": row["claim_allowed"],
            }
            for row in dossier["source_rows"]
        ],
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
    }


def _known_positive_dossiers() -> list[dict[str, Any]]:
    kcsa = _read_json(V36_ROOT / "KcsA" / "evidence_dossier.json")
    xcl1 = _read_json(V36_ROOT / "XCL1_lymphotactin" / "evidence_dossier.json")
    snca = _read_json(V36_ROOT / "alpha_synuclein_SNCA" / "evidence_dossier.json")
    return [
        _dossier(
            "TARGET_001",
            "KcsA",
            "known_positive",
            "membrane_pore_filter_oligomeric_ion_selectivity",
            list(kcsa["bucket_status"]["present_buckets"]),
            [
                "pH-gated potassium channel identity with K+ selectivity",
                "multi-pass membrane topology and pore/filter context",
                "TVGYG-style filter or signature context",
                "oligomeric/interface context present but no coordinate contacts",
            ],
            {"coordinate_contact_tables_allowed": False, "native_metrics_before_selection_allowed": False},
            [
                _source("V36_KCSA_1", "UniProt sequence/features/function annotations", "UniProtKB P0A334", "https://www.uniprot.org/uniprotkb/P0A334/entry", "pH-gated potassium channel identity and K+ preference."),
                _source("V36_KCSA_2", "external sequence conservation signatures", "Potassium-channel filter signature", "UniProtKB P0A334; Heginbotham et al. Biophys J. 1994;66:1061-1067", "Sequence-level filter/signature context."),
            ],
        ),
        _dossier(
            "TARGET_002",
            "XCL1_lymphotactin",
            "known_positive",
            "metamorphic_two_state_fold_switch",
            list(xcl1["bucket_status"]["present_buckets"]),
            [
                "metamorphic two-state protein with chemokine-like monomer state A",
                "beta-sandwich dimer state B with distinct biological function",
                "state-specific function evidence must remain separated",
                "mixed-state pooling is forbidden",
            ],
            {"mixed_state_pooling_allowed": False, "native_metrics_before_selection_allowed": False},
            [
                _source("V36_XCL1_1", "UniProt sequence/features/function annotations", "UniProtKB P47992", "https://www.uniprot.org/uniprotkb/P47992/entry", "C motif chemokine identity."),
                _source("V36_XCL1_2", "literature-derived state/function annotations", "XCL1 metamorphic literature", "Tuinstra RL et al. Proc Natl Acad Sci U S A. 2008;105:5057-5062", "Two unrelated native-state folds with state-specific functions."),
            ],
        ),
        _dossier(
            "TARGET_003",
            "alpha_synuclein_SNCA",
            "known_positive",
            "intrinsic_disorder_contextual_ensemble",
            list(snca["bucket_status"]["present_buckets"]),
            [
                "free-state intrinsic disorder evidence",
                "ensemble-not-single-fold rule",
                "disorder-to-order context when membrane-bound",
                "context-bound helix is not a folding-solved claim",
            ],
            {"single_native_fold_model_allowed": False, "ensemble_not_single_fold": True, "native_metrics_before_selection_allowed": False},
            [
                _source("V36_SNCA_1", "DisProt disorder annotations", "DisProt alpha-synuclein", "https://disprot.org/DP00070", "Free-state intrinsic disorder annotation."),
                _source("V36_SNCA_2", "literature-derived state/function annotations", "Alpha-synuclein membrane-bound helix literature", "Davidson WS et al. J Biol Chem. 1998;273:9443-9449", "Membrane-bound helical transition context."),
            ],
        ),
    ]


def _decoy_dossiers() -> list[dict[str, Any]]:
    return [
        _dossier(
            "TARGET_004",
            "bacteriorhodopsin",
            "membrane_decoy",
            "other_membrane_or_transport_context",
            ["membrane_topology_context", "transport_or_pump_context", "cofactor_context"],
            [
                "multi-pass membrane protein",
                "light-driven proton pump transport context",
                "retinal cofactor context",
                "no potassium filter or K+ ion-selectivity evidence",
            ],
            {"coordinate_contact_tables_allowed": False, "native_metrics_before_selection_allowed": False},
            [_source("DEC_BR_1", "UniProt sequence/features/function annotations", "UniProt bacteriorhodopsin", "UniProtKB P02945; reviewed entry", "Light-driven proton pump membrane protein.")],
        ),
        _dossier(
            "TARGET_005",
            "AQP1",
            "membrane_decoy",
            "other_membrane_or_transport_context",
            ["membrane_topology_context", "water_channel_context", "transport_context"],
            [
                "multi-pass membrane water channel",
                "transport context is water permeability",
                "no potassium selectivity or TVGYG-style filter evidence",
            ],
            {"coordinate_contact_tables_allowed": False, "native_metrics_before_selection_allowed": False},
            [_source("DEC_AQP1_1", "UniProt sequence/features/function annotations", "UniProtKB AQP1", "https://www.uniprot.org/uniprotkb/P29972/entry", "Aquaporin water-channel membrane transport annotation.")],
        ),
        _dossier(
            "TARGET_006",
            "CXCL8",
            "soluble_decoy",
            "soluble_single_or_contextual_fold_not_metamorphic",
            ["soluble_chemokine_context", "receptor_binding_context", "single_state_function_context"],
            [
                "soluble chemokine-like cytokine",
                "receptor-binding inflammatory function",
                "no two-native-state metamorphic switch evidence",
                "no mixed-state pooling rule needed",
            ],
            {"mixed_state_pooling_allowed": False, "native_metrics_before_selection_allowed": False},
            [_source("DEC_CXCL8_1", "UniProt sequence/features/function annotations", "UniProtKB CXCL8", "https://www.uniprot.org/uniprotkb/P10145/entry", "Interleukin-8 chemokine function annotation without metamorphic two-state evidence.")],
        ),
        _dossier(
            "TARGET_007",
            "ubiquitin",
            "soluble_decoy",
            "soluble_single_or_contextual_fold_not_metamorphic",
            ["soluble_single_fold_context", "compact_fold_context", "protein_modifier_context"],
            [
                "small soluble protein modifier",
                "compact single-fold context",
                "no two-state metamorphic evidence",
                "no intrinsic-disorder free-state ensemble evidence",
            ],
            {"single_native_fold_model_allowed": True, "native_metrics_before_selection_allowed": False},
            [_source("DEC_UB_1", "UniProt sequence/features/function annotations", "UniProtKB ubiquitin", "https://www.uniprot.org/uniprotkb/P0CG47/entry", "Ubiquitin modifier annotation as compact soluble contrast.")],
        ),
        _dossier(
            "TARGET_008",
            "lysozyme",
            "folded_decoy",
            "soluble_single_or_contextual_fold_not_metamorphic",
            ["soluble_enzyme_context", "compact_fold_context", "single_state_function_context"],
            [
                "soluble enzyme context",
                "compact folded protein contrast",
                "no free-state intrinsic disorder ensemble evidence",
                "no membrane-bound disorder-to-order rule",
            ],
            {"single_native_fold_model_allowed": True, "native_metrics_before_selection_allowed": False},
            [_source("DEC_LYZ_1", "UniProt sequence/features/function annotations", "UniProtKB lysozyme", "https://www.uniprot.org/uniprotkb/P00698/entry", "Lysozyme enzyme annotation as folded soluble contrast.")],
        ),
        _dossier(
            "TARGET_009",
            "myoglobin",
            "folded_decoy",
            "soluble_single_or_contextual_fold_not_metamorphic",
            ["soluble_heme_protein_context", "compact_fold_context", "single_state_function_context"],
            [
                "soluble heme oxygen-storage protein",
                "compact contextual fold contrast",
                "no intrinsic disorder/free-state ensemble evidence",
                "no metamorphic two-native-state evidence",
            ],
            {"single_native_fold_model_allowed": True, "native_metrics_before_selection_allowed": False},
            [_source("DEC_MYO_1", "UniProt sequence/features/function annotations", "UniProtKB myoglobin", "https://www.uniprot.org/uniprotkb/P02144/entry", "Myoglobin soluble heme-protein annotation as folded contrast.")],
        ),
    ]


def build_panel() -> dict[str, Any]:
    dossiers = _known_positive_dossiers() + _decoy_dossiers()
    answer_key = {
        row["panel_id"]: {
            "target": row["target"],
            "group": row["group"],
            "expected_mechanism_class": row["expected_mechanism_class_for_answer_key_only"],
            "known_positive": row["group"] == "known_positive",
        }
        for row in dossiers
    }
    masked = [_masked_dossier(row) for row in dossiers]
    return {
        "dossiers": dossiers,
        "masked_dossiers": masked,
        "answer_key": answer_key,
        "panel_manifest": {
            "kind": "V38_BLIND_MECHANISM_GENERALIZATION_PANEL_MANIFEST_v0",
            "panel_target_count": len(dossiers),
            "masked_panel_used": True,
            "answer_key_used_for_assignment": False,
            "known_positive_count": sum(1 for row in dossiers if row["group"] == "known_positive"),
            "decoy_count": sum(1 for row in dossiers if row["group"] != "known_positive"),
            "masked_ids": [row["panel_id"] for row in dossiers],
            "claim_allowed": False,
            "new_MD_allowed": False,
            "folding_problem_solved": False,
        },
        "acquisition_log": {
            "kind": "V38_ACQUISITION_LOG_v0",
            "acquisition_mode": "small_curated_external_noncoordinate_blind_panel",
            "source_rules": "No coordinate contacts, no internal runtime evidence, no predicted coordinates, no native metrics.",
            "replacements": [],
            "decoy_design": [
                "Two membrane decoys: bacteriorhodopsin, AQP1.",
                "Two chemokine/small soluble decoys: CXCL8, ubiquitin.",
                "Two folded non-IDP decoys: lysozyme, myoglobin.",
            ],
        },
    }


def write_panel(root: Path = PANEL_ROOT) -> dict[str, Any]:
    panel = build_panel()
    dossiers_dir = root / "dossiers"
    masked_dir = root / "masked_inputs"
    results_dir = root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    for dossier in panel["dossiers"]:
        _write_json(dossiers_dir / f"{dossier['panel_id']}_{dossier['target']}.json", dossier)
    for dossier in panel["masked_dossiers"]:
        _write_json(masked_dir / f"{dossier['masked_id']}.json", dossier)
    _write_json(root / "answer_key.json", panel["answer_key"])
    _write_json(root / "panel_manifest.json", panel["panel_manifest"])
    _write_json(root / "acquisition_log.json", panel["acquisition_log"])
    return {
        "kind": "V38_BLIND_MECHANISM_GENERALIZATION_PANEL_BUILD_v0",
        "panel_target_count": panel["panel_manifest"]["panel_target_count"],
        "known_positive_count": panel["panel_manifest"]["known_positive_count"],
        "decoy_count": panel["panel_manifest"]["decoy_count"],
        "masked_panel_used": True,
        "answer_key_used_for_assignment": False,
        "panel_root": str(root),
        "masked_input_dir": str(masked_dir),
        "answer_key": str(root / "answer_key.json"),
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V38 blind mechanism generalization panel.")
    parser.add_argument("--panel-root", type=Path, default=PANEL_ROOT)
    args = parser.parse_args()
    print(json.dumps(write_panel(args.panel_root), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
