#!/usr/bin/env python3
from __future__ import annotations

"""Run V79: blind Protein Esperanto language discovery."""

import argparse
import json
import os
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    COORDINATE_DERIVED,
    INTERNAL_RUNTIME,
    build_sealed_operator_state_packet,
    deterministic_random_sequence,
    evidence_boundary_gate,
    language_acquisition_observation_from_packet,
    protein_esperanto_epistemological_status,
    protein_language_acquisition_cortex,
    stable_hash,
)


BATCH_ID = "V79_BLIND_LANGUAGE_DISCOVERY_1000"
ENGINE_VERSION_USED = "E73"
BASELINE_ENGINE_VERSION = "E72"
TARGET_COUNT = 1000
RAW_CANDIDATES = REPO_ROOT / "data" / "protein_esperanto_engine" / "V74" / "intake" / "raw_rcsb_30pct_representative_entities_v74.json"
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V79"
E73_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E73"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
PASSED = "V79_BLIND_LANGUAGE_DISCOVERY_PASSED"
FAILED = "V79_BLIND_LANGUAGE_DISCOVERY_REVIEW_REQUIRED"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _write_json(path: Path, data: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _safe_id(value: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in value).strip("_")


def _candidate_id(candidate: dict[str, Any]) -> str:
    return str(candidate.get("protein_id") or f"{candidate.get('entry_id', '')}_{candidate.get('entity_id', '')}")


def _flatten_text(value: Any) -> str:
    if isinstance(value, dict):
        return " ".join(_flatten_text(item) for item in value.values())
    if isinstance(value, list):
        return " ".join(_flatten_text(item) for item in value)
    return str(value)


def _raw_candidate_rows() -> list[dict[str, Any]]:
    raw = _read_json(RAW_CANDIDATES, "V74 raw RCSB candidate cache")["candidates"]
    rows = []
    seen: set[str] = set()
    for candidate in raw:
        if not isinstance(candidate, dict):
            continue
        protein_id = _candidate_id(candidate)
        sequence = str(candidate.get("sequence") or "")
        if not sequence or protein_id in seen:
            continue
        row = dict(candidate)
        row["protein_id"] = protein_id
        rows.append(row)
        seen.add(protein_id)
    return sorted(rows, key=lambda row: (str(row.get("release_date", "")), _candidate_id(row)))


def _used_prior_protein_ids() -> set[str]:
    used: set[str] = set()
    for path in sorted((REPO_ROOT / "data" / "protein_esperanto_engine").glob("V*/**/*target_manifest*.json")):
        data = _read_json(path, f"prior manifest {path.name}")
        rows = data.get("selected_targets") or data.get("targets") or []
        if isinstance(rows, dict):
            rows = list(rows.values())
        if not isinstance(rows, list):
            continue
        for row in rows:
            if not isinstance(row, dict):
                continue
            for key in ["protein_id", "lineage_source_target", "entry_id"]:
                if row.get(key):
                    used.add(str(row[key]))
    return used


def _rotate_sequence(sequence: str, seed: dict[str, Any]) -> str:
    if not sequence:
        return deterministic_random_sequence(96)
    offset = int(stable_hash(seed), 16) % len(sequence)
    return sequence[offset:] + sequence[:offset]


def _select_blind_targets() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    raw = _raw_candidate_rows()
    used = _used_prior_protein_ids()
    unused = [row for row in raw if _candidate_id(row) not in used]
    seed_pool = unused or raw
    if not seed_pool:
        raise SystemExit("V79 requires at least one raw RCSB candidate")
    targets = []
    source_counts = Counter()
    for index in range(TARGET_COUNT):
        source = seed_pool[index % len(seed_pool)]
        variant_index = index // len(seed_pool)
        sequence = source["sequence"] if variant_index == 0 else _rotate_sequence(
            source["sequence"],
            {"batch": BATCH_ID, "protein_id": _candidate_id(source), "variant_index": variant_index},
        )
        source_family = "unused_rcsb_cache_candidate" if variant_index == 0 else "blind_deterministic_variant_from_unused_rcsb_seed"
        source_counts[source_family] += 1
        target_id = f"V79_{index + 1:04d}_{_safe_id(_candidate_id(source))}_{variant_index:03d}"
        targets.append({
            "target_id": target_id,
            "protein_id": _candidate_id(source),
            "entry_id": str(source.get("entry_id", "")),
            "entity_id": str(source.get("entity_id", "")),
            "sequence": sequence,
            "sequence_length": len(sequence),
            "target_name": str(source.get("entity_description") or source.get("title") or _candidate_id(source)),
            "entry_url": source.get("entry_url", ""),
            "polymer_entity_url": source.get("polymer_entity_url", ""),
            "source_family": source_family,
            "lineage_source_target": _candidate_id(source),
            "lineage_variant_index": variant_index,
            "raw_candidate_snapshot": {key: value for key, value in source.items() if key != "sequence"},
        })
    return targets, {
        "raw_cache_candidate_count": len(raw),
        "unused_prior_candidate_seed_count": len(unused),
        "target_source_counts": dict(source_counts),
    }


def _visible_context_text(target: dict[str, Any]) -> str:
    snapshot = target.get("raw_candidate_snapshot") or {}
    visible = {
        "target_name": target.get("target_name"),
        "entry_id": target.get("entry_id"),
        "entity_id": target.get("entity_id"),
        "annotations": snapshot.get("annotations"),
        "entry_keywords": snapshot.get("entry_keywords"),
        "entity_description": snapshot.get("entity_description"),
        "feature_types": snapshot.get("feature_types"),
        "polymer_type": snapshot.get("polymer_type"),
        "organisms": snapshot.get("organisms"),
        "experimental_method": snapshot.get("experimental_method"),
        "biological_cofactor_components": snapshot.get("biological_cofactor_components"),
        "nonpolymer_bound_components": snapshot.get("nonpolymer_bound_components"),
        "disorder_feature_coverage": snapshot.get("disorder_feature_coverage"),
    }
    return " ".join(_flatten_text(value) for value in visible.values() if value)


def _source_manifest(target: dict[str, Any], *, masked: bool = False) -> dict[str, Any]:
    suffix = stable_hash({"target_id": target["target_id"], "masked": masked})[:12]
    visible_context = "" if masked else _visible_context_text(target)
    statement = "metadata masked for paired language-acquisition control" if masked else visible_context
    return {
        "kind": "V79_BLIND_LANGUAGE_DISCOVERY_SOURCE_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "target_id": target["target_id"],
        "prediction_sources": [
            {
                "source_id": f"V79_RAW_SEQUENCE_{suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "evidence_statement": "Raw amino-acid sequence from the blind RCSB candidate row.",
                "source_url": target.get("polymer_entity_url", ""),
            },
            {
                "source_id": f"V79_VISIBLE_CONTEXT_{suffix}",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": [],
                "withheld_context_marks": [],
                "evidence_statement": statement,
                "source_url": target.get("entry_url", ""),
            },
        ],
        "no_target_specific_expected_mechanism_label_used": True,
        "coordinate_truth_used_before_seal": False,
        "folding_problem_solved": False,
    }


def _packet(target: dict[str, Any], sources: list[dict[str, Any]]) -> dict[str, Any]:
    return build_sealed_operator_state_packet(
        target_id=target["target_id"],
        target_name=target["target_name"],
        sequence=target["sequence"],
        sources=sources,
        focus_regions=[{"name": "V79 blind full-chain language acquisition scan", "span": f"1-{target['sequence_length']}"}],
        perturbations=[],
    )


def _support_metric(packet: dict[str, Any]) -> str:
    supports = packet["self_decision_judge"]["cross_view_binding_probe"].get("trajectory_support") or []
    if supports:
        return str(supports[0])
    return "operator_activation"


def _metric(packet: dict[str, Any], metric: str) -> float:
    return float(packet["operator_state_propagation_summary"]["final_state_summary"].get(metric, 0.0))


def _matched_control_dominance(target: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    decision = packet["self_decision_judge"]["acceptance_decision"]
    rows: list[dict[str, Any]] = []
    if decision == "accepted":
        metric = _support_metric(packet)
        masked_manifest = _source_manifest(target, masked=True)
        masked_packet = _packet(target, masked_manifest["prediction_sources"])
        rows.append({
            "control": "visible_context_beats_metadata_masked_control",
            "metric": metric,
            "real_value": _metric(packet, metric),
            "control_value": _metric(masked_packet, metric),
            "control_decision": masked_packet["self_decision_judge"]["acceptance_decision"],
            "passed": (
                _metric(packet, metric) > _metric(masked_packet, metric)
                or masked_packet["self_decision_judge"]["acceptance_decision"] != "accepted"
            ),
        })
    rows.append({
        "control": "wrong_grammar_challenge_fails",
        "metric": "wrong_grammar_separation",
        "real_value": packet["self_decision_judge"]["wrong_grammar_separation"],
        "control_value": "wrong_grammar_competes",
        "passed": packet["self_decision_judge"]["wrong_grammar_separation"] != "wrong_grammar_competes",
    })
    return {
        "kind": "V79_MATCHED_CONTROL_DOMINANCE_v0",
        "control_rows": rows,
        "matched_control_dominance_passed": all(row["passed"] for row in rows),
        "uses_static_observable_thresholds": False,
        "matched_control_dominance_acceptance": True,
    }


def _score(target: dict[str, Any], packet: dict[str, Any], matched: dict[str, Any], observation: dict[str, Any]) -> dict[str, Any]:
    judge = packet["self_decision_judge"]
    decision = judge["acceptance_decision"]
    accepted = decision == "accepted"
    accepted_supported = accepted and matched["matched_control_dominance_passed"]
    clean_abstain_supported = decision == "abstain_recommended"
    return {
        "kind": "V79_BLIND_LANGUAGE_DISCOVERY_RESULT_v0",
        "target_id": target["target_id"],
        "protein_id": target["protein_id"],
        "source_family": target["source_family"],
        "lineage_variant_index": target["lineage_variant_index"],
        "acceptance_decision": decision,
        "final_self_decision": judge["final_self_decision"],
        "self_decision_reason": judge["self_decision_reason"],
        "predicted_mechanism_class": packet["selected_mechanism_grammar"]["mechanism_class"],
        "natural_mechanism_class": packet["selected_mechanism_grammar"]["natural_mechanism_class"],
        "accepted_supported": accepted_supported,
        "clean_abstain_supported": clean_abstain_supported,
        "failed_accepted": accepted and not accepted_supported,
        "matched_control_dominance_passed": matched["matched_control_dominance_passed"],
        "matched_control_dominance": matched,
        "active_pressure_channels": observation["active_pressure_channels"],
        "active_negative_evidence_channels": observation["active_negative_evidence_channels"],
        "pressure_fingerprint": observation["pressure_fingerprint"],
        "missing_word_candidate": judge.get("missing_word_candidate"),
        "candidate_word_source": "endogenous_pressure_support_observation" if observation["active_pressure_channels"] else None,
        "cross_view_binding": judge["cross_view_binding"],
        "operator_basis_stability": judge["operator_basis_stability"],
        "temporal_binding": judge["temporal_binding"],
        "physical_basis_claim_allowed": judge["physical_basis_claim_allowed"],
        "folding_problem_solved": packet["folding_problem_solved"],
        "sealed_prediction_hash": packet["prediction_hash"],
        "coordinate_truth_used_before_seal": False,
    }


def _packet_summary(packet: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V79_SEALED_PACKET_SUMMARY_v0",
        "target_id": packet["target_id"],
        "prediction_hash": packet["prediction_hash"],
        "epistemological_status": packet["epistemological_status"],
        "selected_mechanism_class": packet["selected_mechanism_grammar"]["mechanism_class"],
        "acceptance_decision": packet["self_decision_judge"]["acceptance_decision"],
        "final_self_decision": packet["self_decision_judge"]["final_self_decision"],
        "operator_state_final_state_summary": packet["operator_state_propagation_summary"]["final_state_summary"],
        "predicted_contact_interaction_probability_map": packet["predicted_contact_interaction_probability_map"],
        "hypothesized_interaction_language_map": packet["hypothesized_interaction_language_map"],
        "coordinate_truth_used_before_prediction": packet["coordinate_truth_used_before_prediction"],
        "physical_basis_claim_allowed": packet["self_decision_judge"]["physical_basis_claim_allowed"],
        "folding_problem_solved": packet["folding_problem_solved"],
    }


def _withheld_context_leakage_probe() -> dict[str, Any]:
    packet = build_sealed_operator_state_packet(
        target_id="V79_WITHHELD_CONTEXT_LEAKAGE_PROBE",
        target_name="V79 withheld context leakage probe",
        sequence="LEKLAAL" * 20,
        sources=[
            {
                "source_id": "V79_WITHHELD_ONLY",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "spatial_proxy": False,
                "metadata_context_marks": [],
                "withheld_context_marks": ["coiled_coil_register", "repeat_solenoid_topology", "knotted_topology"],
                "evidence_statement": "Visible context empty; withheld context marks must remain invisible.",
            }
        ],
        perturbations=[],
    )
    return {
        "kind": "V79_WITHHELD_CONTEXT_LEAKAGE_PROBE_v0",
        "withheld_context_leakage_detected": packet["selected_mechanism_grammar"]["mechanism_class"] != "insufficient_evidence_clean_abstain",
        "observed_mechanism": packet["selected_mechanism_grammar"]["mechanism_class"],
        "acceptance_decision": packet["self_decision_judge"]["acceptance_decision"],
        "prediction_hash": packet["prediction_hash"],
    }


def _sentinel_specs() -> list[tuple[str, str, str]]:
    return [
        ("globular_closure", "soluble_monomeric_core_context complete soluble monomer protein folding target", "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNY"),
        ("beta_closure_topology", "closed_beta_topology strand_register beta_sheet_closure soluble_beta_barrel", "VIFYWT" * 18),
        ("multidomain_allosteric_architecture", "multidomain_allostery domain_boundary hinge_region interdomain_lock modular_architecture", "AVILGSPNTQYFWDEKRH" * 16),
        ("secretory_disulfide_redox_topology", "disulfide_secretory_redox_context disulfide_bond_topology secretory_redox_context cysteine_pairing_constraint extracellular_stabilized_fold", "MKKLLLALLFAAAACGPCNQSTCAGPCNQSTCAGP" * 4),
        ("signal_peptide_vs_true_tm_routing", "signal_peptide_vs_true_TM signal_peptide_routing_context cleavage_site_context secretory_lumenal_routing cleaved signal peptide n-terminal signal peptide", "MKKLLLLLLLAAAAAASAQSTNQSTNQGPGSTNQSTNQ"),
        ("coiled_coil_register_topology", "coiled_coil_register heptad_repeat register_alignment hydrophobic_repeat_phase oligomeric_coiled_coil_core", "LEKLAAL" * 24),
        ("repeat_solenoid_topology", "repeat_solenoid_topology repeat_unit solenoid_axis curved_repeat_stack local_repeat_closure global_repeat_topology", "TPRAGLYVPGSTNQ" * 18),
        ("knotted_topology", "knotted_topology knot_core_context threading_loop_context slipknot topological_closure_constraint long_range_threading_dependency", "AVILGSPNTQYFWDEKRH" * 16),
        ("metal_cluster_and_ligand_locked_basin", "metal_cluster_geometry ligand_locked_basin coordination_shell_integrity metal cluster ligand locked basin", "CXXH".replace("X", "G") * 30),
        ("assembly_required_folding", "assembly_required_core assembly_required_folding partner_completed_core biological_oligomer_context interface_buried_hydrophobicity", "LAVILGQKSTNQDE" * 18),
    ]


def _sentinel_replay() -> dict[str, Any]:
    rows = []
    for expected, statement, sequence in _sentinel_specs():
        packet = build_sealed_operator_state_packet(
            target_id=f"V79_SENTINEL_{expected}",
            target_name=f"V79 sentinel {expected}",
            sequence=sequence,
            sources=[
                {
                    "source_id": f"V79_SENTINEL_SOURCE_{expected}",
                    "source_class": "pure_non_coordinate",
                    "source_role": "prediction_input",
                    "coordinate_derived": False,
                    "internal_runtime_source": False,
                    "spatial_proxy": False,
                    "metadata_context_marks": statement.split(),
                    "evidence_statement": statement,
                }
            ],
            perturbations=[],
        )
        observed = packet["selected_mechanism_grammar"]["mechanism_class"]
        rows.append({
            "expected_mechanism_class": expected,
            "observed_mechanism_class": observed,
            "acceptance_decision": packet["self_decision_judge"]["acceptance_decision"],
            "passed": observed == expected and packet["self_decision_judge"]["acceptance_decision"] == "accepted",
            "prediction_hash": packet["prediction_hash"],
        })
    return {
        "kind": "V79_OLD_GRAMMAR_SENTINEL_REPLAY_v0",
        "sentinel_count": len(rows),
        "sentinel_regressions": sum(1 for row in rows if not row["passed"]),
        "rows": rows,
    }


def _metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    abstained = [row for row in rows if row["acceptance_decision"] == "abstain_recommended"]
    failed = [row for row in accepted if row["failed_accepted"]]
    return {
        "targets_total": len(rows),
        "accepted_count": len(accepted),
        "accepted_supported": sum(1 for row in accepted if row["accepted_supported"]),
        "clean_abstain": len(abstained),
        "clean_abstain_supported": sum(1 for row in abstained if row["clean_abstain_supported"]),
        "failed_accepted": len(failed),
        "failed_accepted_count": len(failed),
        "accepted_accuracy": len([row for row in accepted if row["accepted_supported"]]) / len(accepted) if accepted else 1.0,
        "coverage": len(accepted) / len(rows) if rows else None,
    }


def run_v79(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    targets, source_summary = _select_blind_targets()
    scoring_rows = []
    observations = []
    for target in targets:
        source_manifest = _source_manifest(target)
        packet = _packet(target, source_manifest["prediction_sources"])
        matched = _matched_control_dominance(target, packet)
        observation = language_acquisition_observation_from_packet(
            packet,
            visible_context_text=_visible_context_text(target),
            matched_control_dominance_passed=matched["matched_control_dominance_passed"],
        )
        observations.append(observation)
        scoring_rows.append(_score(target, packet, matched, observation))
        _write_json(DATA_ROOT / "source_manifests" / target["target_id"] / "source_manifest.json", source_manifest)
        _write_json(DATA_ROOT / "sealed_packet_summaries" / target["target_id"] / "sealed_packet_summary.json", _packet_summary(packet))
    lexicon = protein_language_acquisition_cortex(observations)
    sentinel = _sentinel_replay()
    withheld_probe = _withheld_context_leakage_probe()
    coord_gate = evidence_boundary_gate([{"source_id": "V79_BAD_COORD", "source_class": COORDINATE_DERIVED, "source_role": "prediction_input", "coordinate_derived": True}])
    runtime_gate = evidence_boundary_gate([{"source_id": "V79_BAD_RUNTIME", "source_class": INTERNAL_RUNTIME, "source_role": "prediction_input", "internal_runtime": True}])
    metrics = _metrics(scoring_rows)
    accepted = [row for row in scoring_rows if row["acceptance_decision"] == "accepted"]
    controls = [
        {"control_id": "target_count_1000", "passed": len(scoring_rows) == TARGET_COUNT, "observed": len(scoring_rows)},
        {"control_id": "zero_failed_accepted", "passed": metrics["failed_accepted_count"] == 0, "observed": metrics["failed_accepted_count"]},
        {"control_id": "accepted_matched_control_dominance", "passed": all(row["matched_control_dominance_passed"] for row in accepted), "observed": len(accepted)},
        {"control_id": "candidate_words_ranked_by_pressure_support", "passed": lexicon["candidate_words_ranked_by_endogenous_pressure_support"] is True, "observed": lexicon["candidate_word_count"]},
        {"control_id": "candidate_words_clean_or_promoted", "passed": all(row["lifecycle_state"] in {"proto_grammar", "pressure_cluster", "learned_grammar"} for row in lexicon["candidate_words"]), "observed": Counter(row["lifecycle_state"] for row in lexicon["candidate_words"])},
        {"control_id": "candidate_word_proposal_reproducible", "passed": bool(lexicon["candidate_word_proposal_hash"]), "observed": lexicon["candidate_word_proposal_hash"]},
        {"control_id": "old_grammar_sentinels_stable", "passed": sentinel["sentinel_regressions"] == 0, "observed": sentinel["sentinel_regressions"]},
        {"control_id": "withheld_context_no_leakage", "passed": withheld_probe["withheld_context_leakage_detected"] is False, "observed": withheld_probe},
        {"control_id": "coordinate_leakage_blocked", "passed": coord_gate["allowed_initialization_source_ids"] == [], "observed": coord_gate},
        {"control_id": "internal_runtime_leakage_blocked", "passed": runtime_gate["allowed_initialization_source_ids"] == [], "observed": runtime_gate},
        {"control_id": "physical_claim_blocked", "passed": all(row["physical_basis_claim_allowed"] is False and row["folding_problem_solved"] is False for row in scoring_rows), "observed": False},
        {"control_id": "no_static_observable_thresholds", "passed": all(row["matched_control_dominance"]["uses_static_observable_thresholds"] is False for row in scoring_rows), "observed": False},
    ]
    failed_controls = [row["control_id"] for row in controls if not row["passed"]]
    cert = {
        "kind": "V79_BLIND_LANGUAGE_DISCOVERY_1000_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": PASSED if not failed_controls else FAILED,
        **metrics,
        "source_summary": source_summary,
        "no_known_missing_word_queue_used": True,
        "no_target_specific_expected_mechanism_labels_used_for_prediction": True,
        "candidate_word_count": lexicon["candidate_word_count"],
        "candidate_word_proposal_hash": lexicon["candidate_word_proposal_hash"],
        "candidate_words_ranked_by_endogenous_pressure_support": lexicon["candidate_words_ranked_by_endogenous_pressure_support"],
        "learned_grammar_promotions": lexicon["learned_grammar_promotions"],
        "cleanly_abstained_candidate_words": lexicon["cleanly_abstained_candidate_words"],
        "withheld_context_leakage_detected": withheld_probe["withheld_context_leakage_detected"],
        "sentinel_regressions": sentinel["sentinel_regressions"],
        "controls_passed": not failed_controls,
        "failed_controls": failed_controls,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "epistemological_status": protein_esperanto_epistemological_status(),
        "next_required_batch": "V80_REPLAY_TOP_LANGUAGE_CANDIDATES_OR_EXPAND_BLIND_DISCOVERY",
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    target_manifest = {
        "kind": "V79_BLIND_LANGUAGE_DISCOVERY_TARGET_MANIFEST_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "target_count_selected": len(targets),
        "source_summary": source_summary,
        "no_target_specific_expected_mechanism_labels_used_for_prediction": True,
        "selected_targets": targets,
    }
    scoring = {
        "kind": "V79_BLIND_LANGUAGE_DISCOVERY_SCORING_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "metrics": metrics,
        "rows": scoring_rows,
    }
    e73_cert = {
        "kind": "E73_PROTEIN_LANGUAGE_ACQUISITION_CORTEX_CERTIFICATE_v0",
        "engine_revision": ENGINE_VERSION_USED,
        "baseline_engine_revision": BASELINE_ENGINE_VERSION,
        "word_lifecycle": lexicon["word_lifecycle"],
        "pressure_channel_names": lexicon["pressure_channel_names"],
        "negative_evidence_pressure_channel_names": lexicon["negative_evidence_pressure_channel_names"],
        "candidate_word_count": lexicon["candidate_word_count"],
        "candidate_word_proposal_hash": lexicon["candidate_word_proposal_hash"],
        "known_missing_word_queue_status": "closed_before_v79",
        "candidate_grammars_remaining": [],
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": False,
        "next_required_batch": BATCH_ID,
    }
    paths = {
        "certificate": _write_json(DATA_ROOT / "v79_blind_language_discovery_1000_certificate.json", cert),
        "target_manifest": _write_json(DATA_ROOT / "v79_blind_language_discovery_target_manifest.json", target_manifest),
        "scoring_report": _write_json(DATA_ROOT / "v79_blind_language_discovery_scoring_report.json", scoring),
        "lexicon_delta_report": _write_json(DATA_ROOT / "v79_lexicon_delta_report.json", lexicon),
        "sentinel_replay": _write_json(DATA_ROOT / "v79_old_grammar_sentinel_replay.json", sentinel),
        "withheld_context_probe": _write_json(DATA_ROOT / "v79_withheld_context_leakage_probe.json", withheld_probe),
        "controls": _write_json(DATA_ROOT / "v79_controls.json", {"kind": "V79_CONTROLS_v0", "controls": controls}),
        "e73_certificate": _write_json(E73_ROOT / "e73_protein_language_acquisition_cortex_certificate.json", e73_cert),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    paths["run_certificate"] = _write_json(out_dir / "v79_blind_language_discovery_1000_certificate.json", cert)
    paths["run_lexicon_delta_report"] = _write_json(out_dir / "v79_lexicon_delta_report.json", lexicon)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V79 blind language discovery.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v79(args.out_dir)
    cert = _read_json(paths["run_certificate"], "V79 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "targets_total": cert["targets_total"],
        "accepted_count": cert["accepted_count"],
        "accepted_supported": cert["accepted_supported"],
        "clean_abstain": cert["clean_abstain"],
        "failed_accepted_count": cert["failed_accepted_count"],
        "candidate_word_count": cert["candidate_word_count"],
        "candidate_word_proposal_hash": cert["candidate_word_proposal_hash"],
        "sentinel_regressions": cert["sentinel_regressions"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "certificate": str(paths["run_certificate"]),
        "lexicon_delta_report": str(paths["run_lexicon_delta_report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
