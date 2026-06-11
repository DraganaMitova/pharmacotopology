#!/usr/bin/env python3
from __future__ import annotations

"""Build the V42 de novo universal mechanism-language challenge panel."""

import argparse
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "de_novo_mechanism_language" / "V42"
PANEL_ROOT = DATA_ROOT / "panel"


def _source(source_id: str, source_type: str, citation: str, statement: str, tags: list[str], regions: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "source_id": source_id,
        "source_type": source_type,
        "source_url_or_citation": citation,
        "source_date_or_version": "accessed_2026-06-11",
        "evidence_statement": statement,
        "feature_tags": tags,
        "region_hints": regions,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "answer_key_source": False,
        "holdout_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_prediction": False,
    }


def _region(name: str, span: str, evidence: str) -> dict[str, str]:
    return {"region_name": name, "span": span, "evidence": evidence}


TARGETS: list[dict[str, Any]] = [
    {
        "target_id": "TARGET_001",
        "target_name": "KcsA",
        "panel_group": "known_anchor",
        "sequence_status": "available",
        "prediction_sources": [
            _source("T001_UNIPROT", "UniProt feature/function text", "UniProtKB P0A334", "Potassium channel with membrane topology and K+ selectivity annotations.", ["membrane", "filter_signature", "ion_selectivity", "oligomer"], [_region("selectivity_filter_signature", "TVGYG motif", "family/function text"), _region("transmembrane_pore", "TM helices", "topology text")]),
            _source("T001_INTERPRO", "InterPro/Pfam family signatures", "InterPro IPR013099 / Pfam PF07885", "Potassium-channel family signature supports pore/filter grammar.", ["membrane", "filter_signature", "ion_selectivity"], [_region("potassium_channel_signature", "family signature", "InterPro/Pfam")]),
        ],
    },
    {
        "target_id": "TARGET_002",
        "target_name": "XCL1_lymphotactin",
        "panel_group": "known_anchor",
        "sequence_status": "available",
        "prediction_sources": [
            _source("T002_LITERATURE", "literature-derived state/function annotations", "Tuinstra RL et al. PNAS 2008; Alderson TR et al. ACS Chem Biol 2015", "XCL1 has separable Ltn10/Ltn40 state-function evidence.", ["metamorphic", "two_state", "state_separated", "function_switch"], [_region("state_A_chemokine_like", "N-terminal chemokine fold context", "state-function text"), _region("state_B_beta_sandwich_dimer", "alternate state context", "state-function text")]),
            _source("T002_REVIEW", "literature-derived state/function annotations", "Dishman AF and Volkman BF. Curr Opin Struct Biol. 2018", "Review-level evidence supports metamorphic fold switching.", ["metamorphic", "two_state", "state_separated"], [_region("state_switch_operator", "global state separation", "review text")]),
        ],
    },
    {
        "target_id": "TARGET_003",
        "target_name": "alpha_synuclein_SNCA",
        "panel_group": "known_anchor",
        "sequence_status": "available",
        "prediction_sources": [
            _source("T003_DISPROT", "DisProt disorder annotations", "DisProt DP00070", "Alpha-synuclein is intrinsically disordered in free state.", ["disorder", "idp", "ensemble"], [_region("acidic_C_terminal_tail", "96-140", "DisProt text"), _region("N_terminal_membrane_binding_region", "1-95", "state text")]),
            _source("T003_LITERATURE", "literature-derived state/function annotations", "Eliezer D et al. J Mol Biol. 2001; Davidson WS et al. JBC 1998", "Membrane context induces conditional helical ordering.", ["disorder", "context_bound_ordering", "membrane_binding"], [_region("membrane_induced_helix", "N-terminal amphipathic region", "literature text")]),
        ],
    },
    {
        "target_id": "TARGET_004",
        "target_name": "4AKE_adenylate_kinase",
        "panel_group": "known_anchor",
        "sequence_status": "available",
        "prediction_sources": [
            _source("T004_UNIPROT", "UniProt feature/function text", "UniProtKB P69441 adenylate kinase", "Soluble adenylate kinase has nucleotide-binding function and compact enzyme domain annotation.", ["compact", "soluble", "single_domain", "closure"], [_region("P_loop_nucleotide_binding", "Walker/P-loop region", "UniProt feature text"), _region("lid_closure_region", "mobile lid context", "function annotation")]),
            _source("T004_INTERPRO", "InterPro/Pfam family signatures", "InterPro adenylate kinase family signatures", "Family signature supports compact kinase core closure grammar.", ["compact", "soluble", "enzyme_core"], [_region("adenylate_kinase_core", "family core", "InterPro/Pfam")]),
        ],
    },
]


def _add(target_name: str, group: str, sources: list[dict[str, Any]], sequence_status: str = "available") -> None:
    TARGETS.append({
        "target_id": f"TARGET_{len(TARGETS) + 1:03d}",
        "target_name": target_name,
        "panel_group": group,
        "sequence_status": sequence_status,
        "prediction_sources": sources,
    })


_add("bacteriorhodopsin_BR", "membrane_channel_transporter", [
    _source("BR_UNIPROT", "UniProt feature/function text", "UniProtKB bacteriorhodopsin annotations", "Seven-transmembrane light-driven proton pump with retinal-binding context.", ["membrane", "transport", "cofactor", "multi_pass"], [_region("retinal_binding_operator", "lysine retinal region", "function text"), _region("seven_tm_bundle", "multi-pass membrane region", "topology text")]),
    _source("BR_INTERPRO", "InterPro/Pfam family signatures", "InterPro bacteriorhodopsin-like family", "Family signature supports membrane transport grammar.", ["membrane", "transport"], [_region("proton_pump_signature", "family signature", "InterPro")]),
])
_add("AQP1_aquaporin", "membrane_channel_transporter", [
    _source("AQP1_UNIPROT", "UniProt feature/function text", "UniProtKB P29972", "Aquaporin water channel with membrane topology and NPA motif annotations.", ["membrane", "channel", "signature_motif", "oligomer"], [_region("NPA_motif_pair", "NPA motifs", "feature text"), _region("membrane_channel_core", "six TM helices", "topology text")]),
])
_add("LacY_lactose_permease", "membrane_channel_transporter", [
    _source("LACY_UNIPROT", "UniProt feature/function text", "UniProtKB P02920", "Major facilitator transporter with multi-pass membrane topology.", ["membrane", "transport", "multi_pass", "alternating_access"], [_region("sugar_transport_operator", "MFS core", "function text"), _region("transmembrane_bundle", "multi-pass region", "topology text")]),
])
_add("CFTR_NBD1", "membrane_channel_transporter", [
    _source("CFTR_UNIPROT", "UniProt feature/function text", "UniProtKB P13569 CFTR", "ABC transporter/channel with nucleotide-binding domain and membrane-channel context.", ["membrane", "transport", "nucleotide_binding", "allosteric"], [_region("NBD1_operator", "NBD1", "feature text"), _region("channel_regulatory_context", "ABC/channel context", "function text")]),
])

for name, citation, region in [
    ("ubiquitin_1UBQ", "UniProtKB P0CG47 ubiquitin", "beta_grasp_core"),
    ("lysozyme_HEWL", "UniProtKB hen egg-white lysozyme", "catalytic_compact_core"),
    ("myoglobin_Mb", "UniProtKB myoglobin annotations", "heme_binding_globin_core"),
    ("barnase", "UniProtKB barnase/ribonuclease annotations", "ribonuclease_core"),
]:
    _add(name, "soluble_compact_single_domain", [
        _source(f"{name.upper()}_UNIPROT", "UniProt feature/function text", citation, "Soluble compact single-domain protein with stable function/family annotation.", ["compact", "soluble", "single_domain"], [_region(region, "domain core", "feature/function text")]),
    ])

for name, statement, tags, region in [
    ("SARS_CoV_2_ORF6", "Short viral accessory protein with sparse family depth and membrane-association annotations.", ["weak_evolutionary", "low_confidence", "membrane_association"], "accessory_region"),
    ("SARS_CoV_2_ORF8", "Rapidly evolving viral accessory protein with shallow evolutionary context.", ["weak_evolutionary", "low_confidence", "immune_context"], "accessory_fold_region"),
    ("p53_TAD", "Transactivation domain has low-complexity/disordered context and weak single-family structural signal.", ["weak_evolutionary", "low_confidence", "disorder", "folding_upon_binding"], "transactivation_domain"),
    ("Ebola_VP35_IID", "Viral interferon inhibitory domain has limited family depth in this local panel.", ["weak_evolutionary", "low_confidence", "contextual_binding"], "interferon_inhibitory_domain"),
]:
    _add(name, "weak_or_shallow_evolutionary_information", [
        _source(f"{name.upper()}_ANNOTATION", "UniProt/literature-derived function annotations", f"{name} feature/function annotations", statement, tags, [_region(region, "annotated region", "feature/function text")]),
    ])

for name, citation, region, extra_tags in [
    ("tau_K18", "DisProt / UniProt tau repeat-region annotations", "microtubule_binding_repeats", ["aggregation_context"]),
    ("FUS_low_complexity_domain", "DisProt / literature FUS low-complexity annotations", "low_complexity_prion_like_region", ["phase_separation"]),
    ("TDP43_C_terminal_LCD", "DisProt / UniProt TDP-43 C-terminal low-complexity annotations", "C_terminal_low_complexity_region", ["aggregation_context"]),
    ("hnRNPA1_LCD", "DisProt / literature hnRNPA1 low-complexity annotations", "low_complexity_RGG_region", ["phase_separation"]),
]:
    _add(name, "disordered_or_partially_disordered", [
        _source(f"{name.upper()}_DISORDER", "DisProt disorder annotations", citation, "Disorder/low-complexity annotation supports ensemble grammar.", ["disorder", "idp", "ensemble", *extra_tags], [_region(region, "annotated disordered region", "DisProt/UniProt text")]),
    ])

for name, citation, regions, tags in [
    ("RfaH_CTD", "RfaH metamorphic fold-switch literature annotations", ["alpha_state_operator", "beta_state_operator"], ["metamorphic", "two_state", "state_separated"]),
    ("Mad2", "Mad2 conformational switch literature annotations", ["open_state_region", "closed_state_region"], ["multistate", "conformational_switch", "state_separated"]),
    ("KaiB", "KaiB fold-switch clock protein literature annotations", ["ground_state_region", "fold_switched_region"], ["metamorphic", "two_state", "state_separated"]),
    ("calmodulin", "UniProt/literature calcium-binding allosteric switch annotations", ["EF_hand_lobes", "target_binding_switch"], ["multistate", "allosteric", "context_bound_ordering"]),
]:
    _add(name, "multistate_allosteric_metamorphic_binding", [
        _source(f"{name.upper()}_STATE", "literature-derived state/function annotations", citation, "State-separated or allosteric function annotations support switch grammar.", tags, [_region(region, "state/operator region", "state/function text") for region in regions]),
    ])


def _blocked_inputs(target: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"blocked_input": "PDB coordinates or mmCIF coordinates before sealing", "reason": "coordinates are post-hoc holdout only"},
        {"blocked_input": "AlphaFold/ESMFold/RoseTTAFold coordinates before sealing", "reason": "predicted coordinates cannot be prediction evidence"},
        {"blocked_input": "native contact maps or coordinate-derived contacts", "reason": "would leak structure truth"},
        {"blocked_input": "first_contact_clean_pharmacotopology_layer_run files as biological evidence", "reason": "runtime outputs are audit evidence only"},
        {"blocked_input": "V33/V34 coordinate-derived KcsA CSVs", "reason": "blocked coordinate-derived contacts"},
        {"blocked_input": "V42 holdout files before prediction sealing", "reason": "holdout must open only after sealed prediction hash"},
        {"blocked_input": "answer key or expected class label", "reason": "assignment must use prediction sources only"},
    ]


def build_panel() -> dict[str, Any]:
    PANEL_ROOT.mkdir(parents=True, exist_ok=True)
    targets = []
    acquisition_log = []
    for target in TARGETS:
        target = {**target, "blocked_prediction_inputs": _blocked_inputs(target)}
        targets.append(target)
        acquisition_log.append({
            "target_id": target["target_id"],
            "target_name": target["target_name"],
            "panel_group": target["panel_group"],
            "acquisition_status": "local_annotation_manifest_built",
            "replacement_reason": None,
            "prediction_source_count": len(target["prediction_sources"]),
        })
        target_dir = PANEL_ROOT / target["target_id"]
        target_dir.mkdir(parents=True, exist_ok=True)
        (target_dir / "prediction_input_manifest.json").write_text(json.dumps(target, indent=2, sort_keys=True), encoding="utf-8")
    panel = {
        "kind": "V42_DE_NOVO_UNIVERSAL_MECHANISM_CHALLENGE_PANEL_v0",
        "panel_target_count": len(targets),
        "panel_groups": {
            group: sum(1 for target in targets if target["panel_group"] == group)
            for group in sorted({target["panel_group"] for target in targets})
        },
        "targets": targets,
        "holdouts_created": False,
        "answer_key_available_to_prediction": False,
        "claim_allowed": False,
        "folding_problem_solved": False,
    }
    (PANEL_ROOT / "panel_manifest.json").write_text(json.dumps(panel, indent=2, sort_keys=True), encoding="utf-8")
    (PANEL_ROOT / "acquisition_log.json").write_text(json.dumps({
        "kind": "V42_PANEL_ACQUISITION_LOG_v0",
        "panel_target_count": len(targets),
        "replacement_count": 0,
        "entries": acquisition_log,
    }, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "kind": "V42_PANEL_BUILD_v0",
        "panel_target_count": len(targets),
        "panel_manifest": str(PANEL_ROOT / "panel_manifest.json"),
        "acquisition_log": str(PANEL_ROOT / "acquisition_log.json"),
        "targets": [target["target_id"] for target in targets],
        "claim_allowed": False,
        "folding_problem_solved": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V42 de novo mechanism challenge panel.")
    parser.parse_args()
    print(json.dumps(build_panel(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
