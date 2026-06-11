#!/usr/bin/env python3
from __future__ import annotations

"""Build V47 RfaH-CTD metamorphic fold-switch source manifests."""

import argparse
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "live_unsolved_targets" / "V47" / "RfaH_CTD"
SOURCE_ROOT = DATA_ROOT / "sources"
PREDICTION_ROOT = DATA_ROOT / "prediction_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts_postseal"
VALIDATION_ROOT = DATA_ROOT / "validation"
REPORT_ROOT = DATA_ROOT / "reports"

RFAH_UNIPROT_ACCESSION = "P0AFW0"
RFAH_UNIPROT_ID = "RFAH_ECOLI"
RFAH_FULL_SEQUENCE = (
    "MQSWYLLYCKRGQLQRAQEHLERQAVNCLAPMITLEKIVRGKRTAVSEPLFPNYLFVEFDPEVIHTTTINATRGVSHFVRFGASPAIVPSAVIHQLSVYK"
    "PKDIVDPATPYPGDKVIITEGAFEGFQAIFTEPDGEARSMLLLNLINKEIKHSVKNTEFRKL"
)
RFAH_CTD_START = 101
RFAH_CTD_END = 162
RFAH_CTD_SEQUENCE = RFAH_FULL_SEQUENCE[RFAH_CTD_START - 1:RFAH_CTD_END]


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _composition(sequence: str, offset: int = 1) -> dict[str, Any]:
    counts = Counter(sequence)
    length = len(sequence)
    hydrophobic = "AILMFWVY"
    charged = "DEKRH"
    polar = "NQST"
    return {
        "length": length,
        "counts": {aa: counts[aa] for aa in sorted(counts)},
        "hydrophobic_fraction": sum(counts[aa] for aa in hydrophobic) / length,
        "charged_fraction": sum(counts[aa] for aa in charged) / length,
        "polar_fraction": sum(counts[aa] for aa in polar) / length,
        "proline_positions": [{"position": offset + index, "residue": aa} for index, aa in enumerate(sequence) if aa == "P"],
        "glycine_positions": [{"position": offset + index, "residue": aa} for index, aa in enumerate(sequence) if aa == "G"],
        "aromatic_positions": [{"position": offset + index, "residue": aa} for index, aa in enumerate(sequence) if aa in "FWY"],
    }


def _domain_architecture() -> list[dict[str, Any]]:
    return [
        {
            "domain": "RfaH_NTD_NGN_like_transcription_domain",
            "start": 1,
            "end": 100,
            "annotation_basis": "UniProt sequence/function/domain text; no coordinate records used",
        },
        {
            "domain": "RfaH_CTD_metamorphic_fold_switch_domain",
            "start": RFAH_CTD_START,
            "end": RFAH_CTD_END,
            "annotation_basis": "RfaH fold-switch literature state labels; no PDB/mmCIF coordinates used",
        },
    ]


def _uniprot_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "UNIPROT_P0AFW0_SEQUENCE_FUNCTION_NONCOORDINATE",
        "source_type": "UniProtKB reviewed sequence, function, interaction, and GO annotation with coordinate records excluded",
        "source_url_or_citation": "https://rest.uniprot.org/uniprotkb/P0AFW0.json",
        "accession": RFAH_UNIPROT_ACCESSION,
        "uniprot_id": RFAH_UNIPROT_ID,
        "retrieved_or_checked": "2026-06-11",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_prediction": False,
        "answer_key_source": False,
        "sequence_region_scope": {
            "protein": "Escherichia coli K-12 transcription antitermination protein RfaH",
            "full_length": len(RFAH_FULL_SEQUENCE),
            "ctd_start": RFAH_CTD_START,
            "ctd_end": RFAH_CTD_END,
            "ctd_sequence": RFAH_CTD_SEQUENCE,
            "full_sequence": RFAH_FULL_SEQUENCE,
        },
        "noncoordinate_annotations_used": [
            "RfaH is a transcription antitermination protein recruited to RNAP complexes by the ops element.",
            "RfaH suppresses pausing and termination in specialized long operons.",
            "RfaH interacts with nontemplate DNA and RNAP.",
            "UniProt GO records include translation activator activity and positive regulation of translation.",
            "The N-terminal domain cavity is described as buried by the C-terminal domain and unmasked in activation context.",
        ],
        "features_excluded_before_prediction": [
            "UniProt PDB cross-references",
            "UniProt secondary-structure features with PDB evidence",
            "PDB/mmCIF coordinate files",
            "coordinate-derived contacts",
            "native contact maps",
            "AlphaFoldDB coordinates",
        ],
    }


def _fold_switch_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "RFAH_CTD_ALPHA_BETA_SWITCH_CELL_2012_NONCOORDINATE",
        "source_type": "literature-derived state/function labels with structural coordinates excluded",
        "source_url_or_citation": "Burmann et al. Cell 2012, PMID 22817892, PMCID PMC3430373, DOI 10.1016/j.cell.2012.05.042",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_prediction": False,
        "answer_key_source": False,
        "state_labels_used": [
            "RfaH-CTD alpha-helical hairpin state",
            "RfaH-CTD beta-barrel or beta-roll state",
            "transcription-factor state",
            "translation-factor state",
            "context-dependent refolding",
        ],
        "coordinate_material_blocked": [
            "PDB IDs and atomic coordinates from the paper are not opened before sealing",
            "no native contacts or coordinate-derived distances are imported",
            "state labels are used as mechanism annotations, not as coordinate truth",
        ],
    }


def _partner_context_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "RFAH_NTD_PARTNER_CONTEXT_FOLD_SWITCH_2021_NONCOORDINATE",
        "source_type": "literature-derived NTD/partner-context mechanism annotation without coordinates",
        "source_url_or_citation": "Galaz-Davison et al. PLoS Comput Biol 2021, PMID 34478435, PMCID PMC8454952, DOI 10.1371/journal.pcbi.1008882",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_prediction": False,
        "answer_key_source": False,
        "context_labels_used": [
            "NTD-bound autoinhibition",
            "CTD release",
            "NTD active role in fold switching",
            "partner/context-dependent refolding",
        ],
        "coordinate_material_blocked": [
            "molecular models, PDB structures, and transition coordinates are excluded before sealing",
            "only text-level mechanism labels are available to prediction",
        ],
    }


def _ctd_sequence_features_source() -> dict[str, Any]:
    return {
        "source_id": "RFAH_CTD_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS",
        "source_type": "derived non-coordinate CTD sequence and domain-scope features",
        "source_url_or_citation": "derived from UniProt P0AFW0 sequence residues 101-162",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_prediction": False,
        "answer_key_source": False,
        "derived_from_source_ids": ["UNIPROT_P0AFW0_SEQUENCE_FUNCTION_NONCOORDINATE"],
        "domain_architecture": _domain_architecture(),
        "composition": {
            "full_length": _composition(RFAH_FULL_SEQUENCE),
            "ctd_101_162": _composition(RFAH_CTD_SEQUENCE, RFAH_CTD_START),
        },
    }


def _problem_context_sources() -> list[dict[str, Any]]:
    return [
        {
            "source_id": "RFAH_ALPHA_BETA_SWITCH_CELL_2012_HOLDOUT",
            "source_type": "post-seal alpha/beta fold-switch literature holdout",
            "source_url_or_citation": "https://europepmc.org/article/MED/22817892",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": [
                "RfaH-CTD is described as switching from an alpha-helical hairpin to a beta-barrel-like translation-factor state.",
                "The holdout is opened only after the sealed packet hash exists.",
            ],
        },
        {
            "source_id": "RFAH_NTD_ACTIVE_ROLE_2021_HOLDOUT",
            "source_type": "post-seal NTD/context fold-switch mechanism holdout",
            "source_url_or_citation": "https://europepmc.org/article/MED/34478435",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["The N-terminal domain is treated as an active contributor to the RfaH fold-switch mechanism."],
        },
        {
            "source_id": "RFAH_INTERDOMAIN_RESIDUES_2025_HOLDOUT",
            "source_type": "post-seal interdomain residue and transition-mechanism holdout",
            "source_url_or_citation": "https://europepmc.org/article/MED/40522227",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["Recent work still tests interdomain residues controlling the RfaH fold switch."],
        },
        {
            "source_id": "RFAH_TRANSCRIPTION_TRANSLATION_COUPLING_HOLDOUT",
            "source_type": "post-seal functional context holdout",
            "source_url_or_citation": "UniProt P0AFW0 GO/function evidence including PMID 22817892 and transcription/translation coupling reviews",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["RfaH is linked to transcription antitermination and translation activation/coupling contexts."],
        },
    ]


def _blocked_prediction_inputs() -> list[dict[str, str]]:
    return [
        {"blocked_input": "PDB coordinates or mmCIF coordinates before sealing", "reason": "coordinates would leak state-specific structural truth"},
        {"blocked_input": "coordinate-derived contacts, distances, or native contact maps", "reason": "V47 tests mechanism grammar, not contact lookup"},
        {"blocked_input": "AlphaFold/ESMFold/RoseTTAFold coordinates before sealing", "reason": "predicted coordinate models are blocked as evidence"},
        {"blocked_input": "UniProt PDB-derived secondary-structure features before sealing", "reason": "helix/beta records with PDB evidence are coordinate-adjacent and excluded"},
        {"blocked_input": "first_contact_clean_pharmacotopology_layer_run files as biological evidence", "reason": "runtime artifacts are audit outputs only"},
        {"blocked_input": "validation holdouts before prediction sealing", "reason": "holdouts open only after sealed prediction hash exists"},
        {"blocked_input": "answer key or pass/fail status during prediction", "reason": "classification must use source-separated allowed evidence"},
        {"blocked_input": "target name only", "reason": "RfaH label alone cannot create a full fold-switch packet"},
        {"blocked_input": "generic transcription-factor annotation only", "reason": "generic function cannot explain CTD alpha/beta metamorphic switching"},
        {"blocked_input": "single consensus fold forcing", "reason": "V47 must reject one-state overpromotion"},
        {"blocked_input": "swapped XCL1/KaiB/Mad2 fold-switch evidence", "reason": "non-RfaH metamorphic examples cannot validate RfaH-specific predictions"},
    ]


def build_sources() -> dict[str, Any]:
    for root in [SOURCE_ROOT, PREDICTION_ROOT, HOLDOUT_ROOT, VALIDATION_ROOT, REPORT_ROOT]:
        root.mkdir(parents=True, exist_ok=True)
    prediction_sources = [
        _uniprot_prediction_source(),
        _fold_switch_prediction_source(),
        _partner_context_prediction_source(),
        _ctd_sequence_features_source(),
    ]
    context_sources = _problem_context_sources()
    blocked_inputs = _blocked_prediction_inputs()
    source_manifest = {
        "kind": "V47_RFAH_CTD_METAMORPHIC_SOURCE_MANIFEST_v0",
        "target_id": "V47_RFAH_CTD",
        "target": "Escherichia coli RfaH C-terminal domain",
        "uniprot_accession": RFAH_UNIPROT_ACCESSION,
        "uniprot_id": RFAH_UNIPROT_ID,
        "sequence_region_scope": f"RfaH {RFAH_UNIPROT_ACCESSION} CTD residues {RFAH_CTD_START}-{RFAH_CTD_END}",
        "prediction_input_source_count": len(prediction_sources),
        "prediction_sources": prediction_sources,
        "problem_context_sources_not_opened_for_prediction": context_sources,
        "blocked_prediction_inputs": blocked_inputs,
        "prediction_inputs_separated_from_validation_holdouts": True,
        "holdouts_created": False,
        "answer_key_available_to_prediction": False,
        "coordinate_sources_available_before_prediction": False,
        "internal_runtime_sources_available_to_prediction": False,
        "folding_problem_solved": False,
        "live_fold_switch_solution_packet": False,
    }
    _write_json(SOURCE_ROOT / "source_manifest.json", source_manifest)
    _write_json(SOURCE_ROOT / "uniprot_p0afw0_prediction_source.json", prediction_sources[0])
    _write_json(SOURCE_ROOT / "rfah_ctd_alpha_beta_switch_prediction_source.json", prediction_sources[1])
    _write_json(SOURCE_ROOT / "rfah_ntd_partner_context_prediction_source.json", prediction_sources[2])
    _write_json(SOURCE_ROOT / "rfah_ctd_sequence_pattern_features.json", prediction_sources[3])
    _write_json(SOURCE_ROOT / "problem_context_sources_not_opened_for_prediction.json", {"sources": context_sources})
    _write_json(SOURCE_ROOT / "blocked_prediction_inputs_manifest.json", {
        "kind": "V47_BLOCKED_PREDICTION_INPUTS_MANIFEST_v0",
        "blocked_prediction_inputs": blocked_inputs,
    })
    _write_json(SOURCE_ROOT / "acquisition_log.json", {
        "kind": "V47_RFAH_CTD_SOURCE_ACQUISITION_LOG_v0",
        "target_id": "V47_RFAH_CTD",
        "external_noncoordinate_source_count": len(prediction_sources) + len(context_sources),
        "prediction_input_source_count": len(prediction_sources),
        "holdout_source_count_before_prediction": 0,
        "coordinate_derived_source_count_before_prediction": 0,
        "internal_runtime_source_count_for_prediction": 0,
        "entries": [
            {"source_id": source["source_id"], "allowed_use": source["allowed_use"]}
            for source in [*prediction_sources, *context_sources]
        ],
    })
    ctd_comp = _composition(RFAH_CTD_SEQUENCE, RFAH_CTD_START)
    return {
        "kind": "V47_RFAH_CTD_SOURCE_BUILD_v0",
        "target_id": "V47_RFAH_CTD",
        "target": "Escherichia coli RfaH C-terminal domain",
        "sequence_region_scope": f"RfaH {RFAH_UNIPROT_ACCESSION} CTD residues {RFAH_CTD_START}-{RFAH_CTD_END}",
        "sequence_length": len(RFAH_FULL_SEQUENCE),
        "ctd_length": len(RFAH_CTD_SEQUENCE),
        "ctd_sequence": RFAH_CTD_SEQUENCE,
        "domain_count": len(_domain_architecture()),
        "ctd_hydrophobic_fraction": ctd_comp["hydrophobic_fraction"],
        "ctd_charged_fraction": ctd_comp["charged_fraction"],
        "prediction_input_source_count": len(prediction_sources),
        "holdouts_created": False,
        "source_manifest": str(SOURCE_ROOT / "source_manifest.json"),
        "blocked_inputs_manifest": str(SOURCE_ROOT / "blocked_prediction_inputs_manifest.json"),
        "folding_problem_solved": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V47 RfaH-CTD metamorphic fold-switch source manifests.")
    parser.parse_args()
    print(json.dumps(build_sources(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
