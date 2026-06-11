#!/usr/bin/env python3
from __future__ import annotations

"""Build V39 mechanism prediction packets and independent holdouts."""

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
PRED_ROOT = REPO_ROOT / "data" / "mechanism_predictions" / "V39"
HOLDOUT_ROOT = REPO_ROOT / "data" / "mechanism_holdouts" / "V39"
VALIDATION_ROOT = REPO_ROOT / "data" / "mechanism_validation" / "V39"

TARGET_ORDER = ["KcsA", "XCL1_lymphotactin", "alpha_synuclein_SNCA"]


def _source(source_id: str, source_type: str, source_name: str, citation: str, statement: str, bucket: str) -> dict[str, Any]:
    return {
        "source_id": source_id,
        "source_type": source_type,
        "source_name": source_name,
        "source_url_or_citation": citation,
        "source_date_or_version": "accessed_2026-06-11",
        "evidence_bucket": bucket,
        "evidence_statement": statement,
        "allowed_use": "V39_holdout_validation_or_prediction_source_no_coordinates_no_MD",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "answer_key_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "claim_allowed": False,
    }


def _prediction(pid: str, operator: str, text: str, holdout_bucket: str, falsifier: str) -> dict[str, Any]:
    return {
        "prediction_id": pid,
        "operator": operator,
        "prediction": text,
        "expected_holdout_bucket": holdout_bucket,
        "falsifiable_if": falsifier,
    }


TARGET_SPECS: dict[str, dict[str, Any]] = {
    "KcsA": {
        "mechanism_class": "membrane_pore_filter_oligomeric_ion_selectivity",
        "prediction_sources": [
            _source("KCSA_PRED_V36_DOSSIER", "V36 evidence dossier", "KcsA V36 evidence dossier", "data/external_evidence_dossiers/KcsA/evidence_dossier.json", "Non-coordinate dossier mechanism evidence.", "prediction_source"),
            _source("KCSA_PRED_V37_MAP", "V37 mechanism map", "KcsA V37 mechanism question map", "data/mechanism_maps/V37/KcsA_mechanism_question_map.json", "Mechanism-class map generated from non-coordinate evidence.", "prediction_source"),
        ],
        "predictions": [
            _prediction("KCSA_P1_FILTER", "filter_signature_required_for_K_selectivity", "A filter/signature operator should be required for K+ selectivity grammar.", "filter_or_signature_holdout", "independent holdout lacks potassium-channel filter/signature evidence"),
            _prediction("KCSA_P2_MEMBRANE", "membrane_topology_required", "The mechanism should require membrane/topology context rather than soluble-core grammar.", "membrane_topology_holdout", "independent holdout supports soluble-only context"),
            _prediction("KCSA_P3_OLIGOMER", "oligomer_context_without_whole_fold_claim", "Oligomer/interface context may matter, but must not become a whole-fold claim.", "oligomer_or_interface_holdout", "holdout contradicts oligomer/interface context or requires whole-fold claim"),
            _prediction("KCSA_P4_PERTURB_FILTER", "filter_ion_evidence_removal_breaks_hard_grammar", "Removing filter/ion-selectivity evidence should make KcsA hard grammar partial or invalid.", "perturbation_holdout", "filter/ion evidence can be removed without changing mechanism grammar"),
        ],
        "holdouts": [
            _source("KCSA_HOLDOUT_INTERPRO_K_CHANNEL", "InterPro family/domain signatures", "InterPro potassium-channel family signature", "InterPro IPR013099 / Pfam PF07885 potassium channel family annotations", "Independent family/signature evidence supports potassium-channel and filter grammar.", "filter_or_signature_holdout"),
            _source("KCSA_HOLDOUT_UNIPROT_FEATURE", "UniProt feature/function annotations", "UniProt KcsA feature/function holdout", "UniProtKB P0A334 feature/function annotations; not copied from prediction packet", "Independent annotation supports K+ preference and membrane/topology context.", "membrane_topology_holdout"),
            _source("KCSA_HOLDOUT_SUBUNIT", "literature-derived state/function annotations", "KcsA homotetramer/interface context", "UniProtKB P0A334 subunit annotation; Heginbotham et al. Biophys J. 1994;66:1061-1067", "Independent non-coordinate function/subunit evidence supports oligomeric channel context.", "oligomer_or_interface_holdout"),
            _source("KCSA_HOLDOUT_FILTER_PERTURBATION", "external sequence/family/conservation signatures", "Potassium-channel filter perturbation expectation", "Heginbotham L et al. Biophys J. 1994;66:1061-1067", "Filter signature/conservation evidence supports that removing filter/ion evidence invalidates hard K+ grammar.", "perturbation_holdout"),
        ],
    },
    "XCL1_lymphotactin": {
        "mechanism_class": "metamorphic_two_state_fold_switch",
        "prediction_sources": [
            _source("XCL1_PRED_V36_DOSSIER", "V36 evidence dossier", "XCL1 V36 evidence dossier", "data/external_evidence_dossiers/XCL1_lymphotactin/evidence_dossier.json", "Non-coordinate two-state dossier evidence.", "prediction_source"),
            _source("XCL1_PRED_V37_MAP", "V37 mechanism map", "XCL1 V37 mechanism question map", "data/mechanism_maps/V37/XCL1_lymphotactin_mechanism_question_map.json", "Mechanism-class map generated from non-coordinate evidence.", "prediction_source"),
        ],
        "predictions": [
            _prediction("XCL1_P1_NO_CONSENSUS", "no_single_consensus_fold", "XCL1 should not be represented as one consensus fold.", "metamorphic_two_state_holdout", "independent holdout supports only a single consensus fold"),
            _prediction("XCL1_P2_TWO_STATES", "two_state_specific_buckets_required", "There must be two state-specific evidence buckets.", "state_A_state_B_holdout", "one state explains all holdout evidence"),
            _prediction("XCL1_P3_DISTINCT_FUNCTION", "state_specific_function_labels_distinct", "State A and state B should have distinct function/context labels.", "state_function_holdout", "holdout merges state functions into one label"),
            _prediction("XCL1_P4_NO_POOLING", "mixed_state_pooling_forbidden", "Mixed-state pooling should be explicitly forbidden.", "pooling_rule_holdout", "mixed-state pooling is valid in holdout evidence"),
            _prediction("XCL1_P5_REMOVE_STATE", "one_state_removal_breaks_mechanism", "Removing one state should make the mechanism partial or invalid.", "perturbation_holdout", "one-state dossier remains high-confidence metamorphic grammar"),
        ],
        "holdouts": [
            _source("XCL1_HOLDOUT_ACS_ENGINEERING", "literature-derived state/function annotations", "Engineered metamorphic chemokine XCL1 holdout", "Alderson TR et al. ACS Chem Biol. 2015;10.1021/acschembio.5b00542", "Independent literature supports engineered XCL1 variants separating Ltn10/Ltn40 state functions.", "metamorphic_two_state_holdout"),
            _source("XCL1_HOLDOUT_STATE_A_B", "literature-derived state/function annotations", "XCL1 state A/state B function evidence", "Tuinstra RL et al. Proc Natl Acad Sci U S A. 2008;105:5057-5062", "Holdout state-function evidence supports chemokine-like monomer and beta-sandwich dimer contexts.", "state_A_state_B_holdout"),
            _source("XCL1_HOLDOUT_FUNCTIONS", "external literature state/function annotations", "XCL1 distinct state-function labels", "Dishman AF and Volkman BF. Curr Opin Struct Biol. 2018;50:90-98", "Independent review-level evidence supports distinct functions for XCL1 states.", "state_function_holdout"),
            _source("XCL1_HOLDOUT_POOLING", "non-coordinate experimental state evidence", "XCL1 no mixed pooling rule holdout", "Alderson TR et al. ACS Chem Biol. 2015;10.1021/acschembio.5b00542", "State-specific engineered variants support forbidding mixed-state pooling.", "pooling_rule_holdout"),
            _source("XCL1_HOLDOUT_PERTURB", "literature-derived state/function annotations", "XCL1 one-state removal perturbation holdout", "Tuinstra RL et al. Proc Natl Acad Sci U S A. 2008;105:5057-5062", "Two-state native interconversion evidence supports partial/invalid status if one state is removed.", "perturbation_holdout"),
        ],
    },
    "alpha_synuclein_SNCA": {
        "mechanism_class": "intrinsic_disorder_contextual_ensemble",
        "prediction_sources": [
            _source("SNCA_PRED_V36_DOSSIER", "V36 evidence dossier", "SNCA V36 evidence dossier", "data/external_evidence_dossiers/alpha_synuclein_SNCA/evidence_dossier.json", "Non-coordinate disorder/ensemble dossier evidence.", "prediction_source"),
            _source("SNCA_PRED_V37_MAP", "V37 mechanism map", "SNCA V37 mechanism question map", "data/mechanism_maps/V37/alpha_synuclein_SNCA_mechanism_question_map.json", "Mechanism-class map generated from non-coordinate evidence.", "prediction_source"),
        ],
        "predictions": [
            _prediction("SNCA_P1_DISORDER", "free_state_disorder_ensemble", "Free protein should be represented as disorder/ensemble, not a single native fold.", "disorder_holdout", "holdout supports compact single native fold for free protein"),
            _prediction("SNCA_P2_CONTEXT_BOUND", "context_bound_structure_only", "Context-bound structure may be allowed only as conditional state.", "context_bound_holdout", "holdout treats bound helix as solved free-state fold"),
            _prediction("SNCA_P3_DISORDER_TO_ORDER", "disorder_to_order_contextual", "Disorder-to-order or membrane-bound helix context should be contextual.", "disorder_to_order_holdout", "holdout lacks context dependence"),
            _prediction("SNCA_P4_BLOCK_SINGLE", "compact_single_fold_forbidden", "Forcing compact single-fold grammar should be blocked.", "single_fold_block_holdout", "single compact native fold grammar is supported"),
            _prediction("SNCA_P5_REMOVE_DISORDER", "disorder_removal_breaks_mechanism", "Removing disorder evidence should make the mechanism partial or invalid.", "perturbation_holdout", "removing disorder evidence leaves high-confidence IDP grammar"),
        ],
        "holdouts": [
            _source("SNCA_HOLDOUT_DISPROT_DP00070", "DisProt disorder annotations", "DisProt DP00070 alpha-synuclein holdout", "https://disprot.org/DP00070", "Independent disorder evidence supports free-state disorder/ensemble grammar.", "disorder_holdout"),
            _source("SNCA_HOLDOUT_MEMBRANE_HELIX", "literature-derived state/function annotations", "Alpha-synuclein membrane-bound helix holdout", "Eliezer D et al. J Mol Biol. 2001;307:1061-1073", "Independent NMR/literature state label supports membrane/SDS-bound helical context, not solved free fold.", "context_bound_holdout"),
            _source("SNCA_HOLDOUT_DAVIDSON", "non-coordinate experimental state evidence", "Alpha-synuclein disorder-to-order holdout", "Davidson WS et al. J Biol Chem. 1998;273:9443-9449", "Independent non-coordinate state evidence supports membrane-induced secondary structure.", "disorder_to_order_holdout"),
            _source("SNCA_HOLDOUT_NO_SINGLE", "DisProt disorder annotations", "Alpha-synuclein no single native fold rule holdout", "DisProt DP00070 / UniProtKB P37840", "Disorder annotation supports blocking compact single-fold grammar.", "single_fold_block_holdout"),
            _source("SNCA_HOLDOUT_PERTURB", "external literature state/function annotations", "SNCA disorder removal perturbation holdout", "Uversky VN et al. Proteins. 2001;45:515-521", "Independent disorder literature supports partial/invalid status when disorder evidence is removed.", "perturbation_holdout"),
        ],
    },
}


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})


def _prediction_packet(target: str) -> dict[str, Any]:
    spec = TARGET_SPECS[target]
    return {
        "kind": "V39_MECHANISM_PREDICTION_PACKET_v0",
        "target": target,
        "mechanism_class": spec["mechanism_class"],
        "prediction_source_evidence": spec["prediction_sources"],
        "falsifiable_predictions": spec["predictions"],
        "answer_key_used_for_prediction": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
    }


def _holdout_manifest(target: str) -> dict[str, Any]:
    spec = TARGET_SPECS[target]
    prediction_source_ids = {row["source_id"] for row in spec["prediction_sources"]}
    return {
        "kind": "V39_MECHANISM_HOLDOUT_MANIFEST_v0",
        "target": target,
        "mechanism_class": spec["mechanism_class"],
        "holdout_validation_evidence": [
            {**row, "independent_from_prediction_sources": row["source_id"] not in prediction_source_ids}
            for row in spec["holdouts"]
        ],
        "prediction_source_ids_excluded": sorted(prediction_source_ids),
        "answer_key_used_for_prediction": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
    }


def write_target(target: str) -> dict[str, str]:
    packet = _prediction_packet(target)
    holdout = _holdout_manifest(target)
    pred_dir = PRED_ROOT / target
    holdout_dir = HOLDOUT_ROOT / target
    validation_dir = VALIDATION_ROOT / target
    validation_dir.mkdir(parents=True, exist_ok=True)
    _write_json(pred_dir / "prediction_packet.json", packet)
    _write_json(holdout_dir / "holdout_manifest.json", holdout)
    _write_csv(
        pred_dir / "prediction_table.csv",
        packet["falsifiable_predictions"],
        ["prediction_id", "operator", "prediction", "expected_holdout_bucket", "falsifiable_if"],
    )
    _write_csv(
        holdout_dir / "holdout_evidence_table.csv",
        holdout["holdout_validation_evidence"],
        ["source_id", "source_type", "source_name", "source_url_or_citation", "source_date_or_version", "evidence_bucket", "evidence_statement", "independent_from_prediction_sources", "coordinate_derived", "internal_runtime_source", "answer_key_source"],
    )
    return {
        "prediction_packet": str(pred_dir / "prediction_packet.json"),
        "prediction_table": str(pred_dir / "prediction_table.csv"),
        "holdout_manifest": str(holdout_dir / "holdout_manifest.json"),
        "holdout_evidence_table": str(holdout_dir / "holdout_evidence_table.csv"),
        "validation_dir": str(validation_dir),
    }


def build_all() -> dict[str, Any]:
    artifacts = {target: write_target(target) for target in TARGET_ORDER}
    return {
        "kind": "V39_MECHANISM_PREDICTION_HOLDOUT_BUILD_v0",
        "target_count": len(TARGET_ORDER),
        "prediction_packet_count": len(TARGET_ORDER),
        "holdout_source_count": sum(len(TARGET_SPECS[target]["holdouts"]) for target in TARGET_ORDER),
        "targets": TARGET_ORDER,
        "artifacts": artifacts,
        "answer_key_used_for_prediction": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V39 mechanism prediction holdouts.")
    parser.parse_args()
    print(json.dumps(build_all(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
