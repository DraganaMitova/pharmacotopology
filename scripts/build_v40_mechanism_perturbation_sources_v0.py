#!/usr/bin/env python3
from __future__ import annotations

"""Build V40 mechanism perturbation pressure source manifests and packets."""

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "mechanism_perturbations" / "V40"
SOURCE_ROOT = DATA_ROOT / "sources"
PRED_ROOT = DATA_ROOT / "predictions"
VALIDATION_ROOT = DATA_ROOT / "validation"

TARGET_ORDER = ["KcsA", "XCL1_lymphotactin", "alpha_synuclein_SNCA"]


def _source(
    source_id: str,
    source_role: str,
    source_type: str,
    source_name: str,
    citation: str,
    statement: str,
    bucket: str,
) -> dict[str, Any]:
    return {
        "source_id": source_id,
        "source_role": source_role,
        "source_type": source_type,
        "source_name": source_name,
        "source_url_or_citation": citation,
        "source_date_or_version": "accessed_2026-06-11",
        "evidence_bucket": bucket,
        "evidence_statement": statement,
        "allowed_use": "V40_mechanism_perturbation_pressure_no_coordinates_no_MD",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "answer_key_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "claim_allowed": False,
    }


def _prediction(
    perturbation_id: str,
    operator: str,
    bucket: str,
    pressure: str,
    behavior: str,
    falsifier: str,
) -> dict[str, Any]:
    return {
        "perturbation_id": perturbation_id,
        "operator": operator,
        "perturbation_bucket": bucket,
        "expected_validation_bucket": bucket,
        "predicted_pressure": pressure,
        "expected_behavior": behavior,
        "falsifiable_if": falsifier,
    }


TARGET_SPECS: dict[str, dict[str, Any]] = {
    "KcsA": {
        "mechanism_class": "membrane_pore_filter_oligomeric_ion_selectivity",
        "scientist_question": "Can the grammar identify that perturbations affecting filter/ion-selectivity evidence should matter more than generic channel annotation?",
        "prediction_sources": [
            _source("KCSA_V40_V36_DOSSIER", "prediction_seed", "V36 evidence dossier", "KcsA V36 non-coordinate dossier", "data/external_evidence_dossiers/KcsA/evidence_dossier.json", "Seed evidence supplies non-coordinate KcsA mechanism context.", "prediction_seed"),
            _source("KCSA_V40_V37_MAP", "prediction_seed", "V37 mechanism map", "KcsA V37 mechanism map", "data/mechanism_maps/V37/KcsA_mechanism_question_map.json", "Mechanism map separates filter, membrane, and oligomer context.", "prediction_seed"),
            _source("KCSA_V40_V39_PACKET", "prediction_seed", "V39 prediction packet", "KcsA V39 falsifiable prediction packet", "data/mechanism_predictions/V39/KcsA/prediction_packet.json", "V39 packet identifies filter and selectivity as causal grammar operators.", "prediction_seed"),
            _source("KCSA_V40_V39_VALIDATION", "prediction_seed", "V39 validation result", "KcsA V39 validation result", "data/mechanism_validation/V39/KcsA/validation_result.json", "V39 validation supported KcsA operator buckets without coordinate evidence.", "prediction_seed"),
        ],
        "perturbations": [
            _prediction("KCSA_V40_P1_FILTER_SIGNATURE", "remove_filter_signature", "filter_signature_perturbation", "break", "Removing or relabeling the filter signature must break hard KcsA grammar.", "generic channel annotation alone preserves hard KcsA grammar"),
            _prediction("KCSA_V40_P2_ION_SELECTIVITY", "remove_ion_selectivity_label", "ion_selectivity_perturbation", "break", "Removing K+ selectivity evidence must break or invalidate the hard selectivity grammar.", "hard KcsA grammar remains validated without ion-selectivity evidence"),
            _prediction("KCSA_V40_P3_MEMBRANE_CONTEXT", "remove_membrane_topology_context", "membrane_topology_context_perturbation", "weaken", "Generic soluble/channel annotation alone must not preserve membrane-pore grammar.", "soluble-only annotation preserves hard KcsA grammar"),
            _prediction("KCSA_V40_P4_OLIGOMER_CONTEXT", "remove_oligomer_interface_context", "oligomer_interface_context_perturbation", "weaken", "Oligomer/interface context may weaken mechanism support but must not become a whole-fold prediction.", "oligomer context is converted into a whole-fold claim"),
        ],
        "validation_holdouts": [
            _source("KCSA_V40_H1_FILTER_SIGNATURE", "validation_holdout", "InterPro/Pfam family signatures", "Potassium-channel filter family signature", "InterPro IPR013099 / Pfam PF07885 potassium channel annotations", "Independent family signatures support filter/signature pressure as causal for hard KcsA grammar.", "filter_signature_perturbation"),
            _source("KCSA_V40_H2_ION_SELECTIVITY", "validation_holdout", "literature-derived mutation/function annotations", "KcsA selectivity filter perturbation evidence", "Heginbotham L et al. Biophys J. 1994;66:1061-1067", "Non-coordinate selectivity evidence supports breakage when K+ selectivity evidence is removed.", "ion_selectivity_perturbation"),
            _source("KCSA_V40_H3_MEMBRANE_CONTEXT", "validation_holdout", "UniProt feature/function annotations", "KcsA membrane topology annotations", "UniProtKB P0A334 feature and function annotations", "Independent membrane annotation supports topology context pressure without coordinate contacts.", "membrane_topology_context_perturbation"),
            _source("KCSA_V40_H4_OLIGOMER_CONTEXT", "validation_holdout", "literature-derived state/function annotations", "KcsA subunit/interface context", "UniProtKB P0A334 subunit annotation; Heginbotham L et al. Biophys J. 1994;66:1061-1067", "Independent subunit/function evidence supports oligomer context while keeping whole-fold claims disabled.", "oligomer_interface_context_perturbation"),
        ],
    },
    "XCL1_lymphotactin": {
        "mechanism_class": "metamorphic_two_state_fold_switch",
        "scientist_question": "Can the grammar identify that perturbing state separation or deleting one state destroys the mechanism?",
        "prediction_sources": [
            _source("XCL1_V40_V36_DOSSIER", "prediction_seed", "V36 evidence dossier", "XCL1 V36 non-coordinate dossier", "data/external_evidence_dossiers/XCL1_lymphotactin/evidence_dossier.json", "Seed evidence supplies two-state XCL1 context.", "prediction_seed"),
            _source("XCL1_V40_V37_MAP", "prediction_seed", "V37 mechanism map", "XCL1 V37 mechanism map", "data/mechanism_maps/V37/XCL1_lymphotactin_mechanism_question_map.json", "Mechanism map separates state-specific operators.", "prediction_seed"),
            _source("XCL1_V40_V39_PACKET", "prediction_seed", "V39 prediction packet", "XCL1 V39 falsifiable prediction packet", "data/mechanism_predictions/V39/XCL1_lymphotactin/prediction_packet.json", "V39 packet identifies two-state separation and no-pooling rules.", "prediction_seed"),
            _source("XCL1_V40_V39_VALIDATION", "prediction_seed", "V39 validation result", "XCL1 V39 validation result", "data/mechanism_validation/V39/XCL1_lymphotactin/validation_result.json", "V39 validation supported XCL1 two-state operator buckets without coordinate evidence.", "prediction_seed"),
        ],
        "perturbations": [
            _prediction("XCL1_V40_P1_STATE_A_LOSS", "remove_state_A_bucket", "state_A_loss_or_weakening", "break", "Removing state A makes the metamorphic mechanism partial or invalid.", "one state explains all XCL1 behavior"),
            _prediction("XCL1_V40_P2_STATE_B_LOSS", "remove_state_B_bucket", "state_B_loss_or_weakening", "break", "Removing state B makes the metamorphic mechanism partial or invalid.", "one state explains all XCL1 behavior"),
            _prediction("XCL1_V40_P3_MIXED_POOLING", "force_mixed_state_pooling", "mixed_state_pooling_error", "block", "Mixed-state pooling must be blocked.", "pooled state evidence remains high-confidence"),
            _prediction("XCL1_V40_P4_FUNCTION_DECOUPLING", "merge_state_function_labels", "state_function_decoupling", "weaken", "State-specific function labels must remain separated.", "state functions can be merged without weakening the mechanism"),
            _prediction("XCL1_V40_P5_BALANCE_SHIFT", "shift_two_state_balance", "two_state_balance_shift", "shift", "Perturbations that shift the balance between states should change mechanism support.", "balance shifts leave the two-state mechanism unchanged"),
        ],
        "validation_holdouts": [
            _source("XCL1_V40_H1_STATE_A", "validation_holdout", "literature-derived state/function annotations", "XCL1 Ltn10 state function evidence", "Tuinstra RL et al. Proc Natl Acad Sci U S A. 2008;105:5057-5062", "Independent state-function evidence supports partial/invalid status when state A is removed.", "state_A_loss_or_weakening"),
            _source("XCL1_V40_H2_STATE_B", "validation_holdout", "literature-derived state/function annotations", "XCL1 Ltn40 state function evidence", "Tuinstra RL et al. Proc Natl Acad Sci U S A. 2008;105:5057-5062", "Independent state-function evidence supports partial/invalid status when state B is removed.", "state_B_loss_or_weakening"),
            _source("XCL1_V40_H3_POOLING", "validation_holdout", "non-coordinate experimental perturbation evidence", "Engineered XCL1 state separation evidence", "Alderson TR et al. ACS Chem Biol. 2015;10.1021/acschembio.5b00542", "State-specific engineered variants support blocking mixed-state pooling.", "mixed_state_pooling_error"),
            _source("XCL1_V40_H4_FUNCTION", "validation_holdout", "external literature state/function annotations", "XCL1 distinct state-function labels", "Dishman AF and Volkman BF. Curr Opin Struct Biol. 2018;50:90-98", "Independent review-level evidence supports separate state-function labels.", "state_function_decoupling"),
            _source("XCL1_V40_H5_BALANCE", "validation_holdout", "literature-derived state/function annotations", "XCL1 two-state balance pressure", "Alderson TR et al. ACS Chem Biol. 2015;10.1021/acschembio.5b00542", "Independent engineering literature supports balance-shift pressure between XCL1 states.", "two_state_balance_shift"),
        ],
    },
    "alpha_synuclein_SNCA": {
        "mechanism_class": "intrinsic_disorder_contextual_ensemble",
        "scientist_question": "Can the grammar identify that forcing SNCA into one compact fold is wrong, and that removing disorder evidence weakens the mechanism?",
        "prediction_sources": [
            _source("SNCA_V40_V36_DOSSIER", "prediction_seed", "V36 evidence dossier", "SNCA V36 non-coordinate dossier", "data/external_evidence_dossiers/alpha_synuclein_SNCA/evidence_dossier.json", "Seed evidence supplies disorder/ensemble context.", "prediction_seed"),
            _source("SNCA_V40_V37_MAP", "prediction_seed", "V37 mechanism map", "SNCA V37 mechanism map", "data/mechanism_maps/V37/alpha_synuclein_SNCA_mechanism_question_map.json", "Mechanism map separates free disorder, context-bound helix, and aggregation context.", "prediction_seed"),
            _source("SNCA_V40_V39_PACKET", "prediction_seed", "V39 prediction packet", "SNCA V39 falsifiable prediction packet", "data/mechanism_predictions/V39/alpha_synuclein_SNCA/prediction_packet.json", "V39 packet identifies disorder and context-bound operators.", "prediction_seed"),
            _source("SNCA_V40_V39_VALIDATION", "prediction_seed", "V39 validation result", "SNCA V39 validation result", "data/mechanism_validation/V39/alpha_synuclein_SNCA/validation_result.json", "V39 validation supported SNCA disorder/ensemble operator buckets without coordinate evidence.", "prediction_seed"),
        ],
        "perturbations": [
            _prediction("SNCA_V40_P1_DISORDER_LOSS", "remove_disorder_evidence", "disorder_evidence_loss", "weaken", "Removing disorder evidence makes the mechanism partial or invalid.", "free SNCA remains high-confidence IDP grammar without disorder evidence"),
            _prediction("SNCA_V40_P2_SINGLE_FOLD", "force_compact_single_fold", "ensemble_to_single_fold_forcing", "block", "Forcing compact single-fold grammar must be blocked.", "compact single-fold grammar validates free SNCA"),
            _prediction("SNCA_V40_P3_BOUND_HELIX", "overpromote_context_bound_helix", "context_bound_helix_overpromotion", "block", "Membrane-bound helix evidence can be conditional only, not a solved free fold.", "bound helix is promoted to solved free-state fold"),
            _prediction("SNCA_V40_P4_AGGREGATION_NATIVE", "overpromote_aggregation_context", "aggregation_context_overpromotion", "block", "Aggregation or amyloid context must not become a native-fold claim.", "aggregation context is treated as native-fold evidence"),
            _prediction("SNCA_V40_P5_CONTEXT_SHIFT", "shift_disorder_to_order_context", "disorder_to_order_context_shift", "shift", "Disorder-to-order shifts should remain context-dependent.", "context dependence is ignored"),
        ],
        "validation_holdouts": [
            _source("SNCA_V40_H1_DISORDER", "validation_holdout", "DisProt disorder annotations", "DisProt DP00070 alpha-synuclein disorder annotation", "https://disprot.org/DP00070", "Independent disorder annotation supports weakening when disorder evidence is removed.", "disorder_evidence_loss"),
            _source("SNCA_V40_H2_SINGLE_FOLD", "validation_holdout", "DisProt and UniProt annotations", "Alpha-synuclein no compact native fold rule", "DisProt DP00070 / UniProtKB P37840", "Independent disorder annotations support blocking compact single-fold grammar.", "ensemble_to_single_fold_forcing"),
            _source("SNCA_V40_H3_BOUND_HELIX", "validation_holdout", "literature-derived state/function annotations", "Alpha-synuclein membrane-bound helix context", "Eliezer D et al. J Mol Biol. 2001;307:1061-1073", "Independent state evidence supports conditional bound-helix context only.", "context_bound_helix_overpromotion"),
            _source("SNCA_V40_H4_AGGREGATION", "validation_holdout", "external literature state/function annotations", "Alpha-synuclein aggregation context boundary", "Uversky VN et al. Proteins. 2001;45:515-521", "Independent disorder/aggregation literature supports blocking amyloid context as a native-fold claim.", "aggregation_context_overpromotion"),
            _source("SNCA_V40_H5_CONTEXT_SHIFT", "validation_holdout", "non-coordinate experimental perturbation evidence", "Alpha-synuclein disorder-to-order membrane context", "Davidson WS et al. J Biol Chem. 1998;273:9443-9449", "Independent non-coordinate state evidence supports context-dependent disorder-to-order shift.", "disorder_to_order_context_shift"),
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


def _source_manifest(target: str) -> dict[str, Any]:
    spec = TARGET_SPECS[target]
    prediction_ids = {row["source_id"] for row in spec["prediction_sources"]}
    holdouts = [
        {**row, "independent_from_prediction_sources": row["source_id"] not in prediction_ids}
        for row in spec["validation_holdouts"]
    ]
    return {
        "kind": "V40_MECHANISM_PERTURBATION_SOURCE_MANIFEST_v0",
        "target": target,
        "mechanism_class": spec["mechanism_class"],
        "scientist_question": spec["scientist_question"],
        "prediction_seed_sources": spec["prediction_sources"],
        "validation_holdout_sources": holdouts,
        "prediction_source_ids_excluded_from_holdouts": sorted(prediction_ids),
        "answer_key_used_for_prediction": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
    }


def _prediction_packet(target: str) -> dict[str, Any]:
    spec = TARGET_SPECS[target]
    return {
        "kind": "V40_MECHANISM_PERTURBATION_PREDICTION_PACKET_v0",
        "target": target,
        "mechanism_class": spec["mechanism_class"],
        "scientist_question": spec["scientist_question"],
        "perturbation_predictions": spec["perturbations"],
        "answer_key_used_for_prediction": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
    }


def write_target(target: str) -> dict[str, str]:
    manifest = _source_manifest(target)
    packet = _prediction_packet(target)
    source_dir = SOURCE_ROOT / target
    pred_dir = PRED_ROOT / target
    validation_dir = VALIDATION_ROOT / target
    validation_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = source_dir / "perturbation_source_manifest.json"
    packet_path = pred_dir / "perturbation_prediction_packet.json"
    table_path = pred_dir / "perturbation_table.csv"
    _write_json(manifest_path, manifest)
    _write_json(packet_path, packet)
    _write_csv(
        table_path,
        packet["perturbation_predictions"],
        [
            "perturbation_id",
            "operator",
            "perturbation_bucket",
            "expected_validation_bucket",
            "predicted_pressure",
            "expected_behavior",
            "falsifiable_if",
        ],
    )
    return {
        "perturbation_source_manifest": str(manifest_path),
        "perturbation_prediction_packet": str(packet_path),
        "perturbation_table": str(table_path),
        "validation_dir": str(validation_dir),
    }


def build_all() -> dict[str, Any]:
    artifacts = {target: write_target(target) for target in TARGET_ORDER}
    return {
        "kind": "V40_MECHANISM_PERTURBATION_SOURCE_BUILD_v0",
        "target_count": len(TARGET_ORDER),
        "perturbation_packet_count": len(TARGET_ORDER),
        "perturbation_count_by_target": {
            target: len(TARGET_SPECS[target]["perturbations"]) for target in TARGET_ORDER
        },
        "validation_holdout_count_by_target": {
            target: len(TARGET_SPECS[target]["validation_holdouts"]) for target in TARGET_ORDER
        },
        "targets": TARGET_ORDER,
        "artifacts": artifacts,
        "answer_key_used_for_prediction": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V40 mechanism perturbation pressure inputs.")
    parser.parse_args()
    print(json.dumps(build_all(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
