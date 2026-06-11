#!/usr/bin/env python3
from __future__ import annotations

"""Build V48 SARS-CoV-2 ORF6 viral accessory source manifests."""

import argparse
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "live_unsolved_targets" / "V48" / "SARS2_ORF6"
SOURCE_ROOT = DATA_ROOT / "sources"
PREDICTION_ROOT = DATA_ROOT / "prediction_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts_postseal"
VALIDATION_ROOT = DATA_ROOT / "validation"
REPORT_ROOT = DATA_ROOT / "reports"

ORF6_UNIPROT_ACCESSION = "P0DTC6"
ORF6_UNIPROT_ID = "NS6_SARS2"
ORF6_FULL_SEQUENCE = "MFHLVDFQVTIAEILLIIMRTFKVSIWNLDYIINLIIKNLSKSLTENKYSQLDEEQPMEID"
ORF6_CTERM_START = 38
ORF6_CTERM_END = 61
ORF6_CTERM_SEQUENCE = ORF6_FULL_SEQUENCE[ORF6_CTERM_START - 1:ORF6_CTERM_END]
ORF6_MINIMAL_CTERM_MOTIF_START = 50
ORF6_MINIMAL_CTERM_MOTIF_END = 61
ORF6_MINIMAL_CTERM_MOTIF = ORF6_FULL_SEQUENCE[ORF6_MINIMAL_CTERM_MOTIF_START - 1:ORF6_MINIMAL_CTERM_MOTIF_END]


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _composition(sequence: str, offset: int = 1) -> dict[str, Any]:
    counts = Counter(sequence)
    length = len(sequence)
    hydrophobic = "AILMFWVY"
    charged = "DEKRH"
    acidic = "DE"
    polar = "NQST"
    return {
        "length": length,
        "counts": {aa: counts[aa] for aa in sorted(counts)},
        "hydrophobic_fraction": sum(counts[aa] for aa in hydrophobic) / length,
        "charged_fraction": sum(counts[aa] for aa in charged) / length,
        "acidic_fraction": sum(counts[aa] for aa in acidic) / length,
        "polar_fraction": sum(counts[aa] for aa in polar) / length,
        "methionine_positions": [{"position": offset + index, "residue": aa} for index, aa in enumerate(sequence) if aa == "M"],
        "acidic_positions": [{"position": offset + index, "residue": aa} for index, aa in enumerate(sequence) if aa in acidic],
        "aromatic_positions": [{"position": offset + index, "residue": aa} for index, aa in enumerate(sequence) if aa in "FWY"],
    }


def _region_windows() -> list[dict[str, Any]]:
    windows = []
    for start, end, label in [
        (1, 17, "N_terminal_hydrophobic_region"),
        (18, 24, "Golgi_localization_region"),
        (25, 37, "central_linker_region"),
        (38, 61, "C_terminal_host_hijacking_region"),
        (50, 61, "minimal_C_terminal_RAE1_NUP98_motif_region"),
    ]:
        sub = ORF6_FULL_SEQUENCE[start - 1:end]
        comp = _composition(sub, start)
        windows.append({
            "label": label,
            "span": f"{start}-{end}",
            "sequence": sub,
            "length": len(sub),
            "hydrophobic_fraction": comp["hydrophobic_fraction"],
            "charged_fraction": comp["charged_fraction"],
            "acidic_fraction": comp["acidic_fraction"],
            "methionine_positions": comp["methionine_positions"],
            "acidic_positions": comp["acidic_positions"],
            "sequence_pattern_note": "non-coordinate ORF6 short-region operator window",
        })
    return windows


def _uniprot_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "UNIPROT_P0DTC6_SEQUENCE_FUNCTION_NONCOORDINATE",
        "source_type": "UniProtKB reviewed sequence, function, interaction, localization, and mutagenesis text with coordinate records excluded",
        "source_url_or_citation": "https://rest.uniprot.org/uniprotkb/P0DTC6.json",
        "accession": ORF6_UNIPROT_ACCESSION,
        "uniprot_id": ORF6_UNIPROT_ID,
        "retrieved_or_checked": "2026-06-11",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_prediction": False,
        "answer_key_source": False,
        "sequence_region_scope": {
            "protein": "SARS-CoV-2 ORF6 / accessory protein 6 / ns6",
            "full_length": len(ORF6_FULL_SEQUENCE),
            "sequence": ORF6_FULL_SEQUENCE,
            "c_terminal_region": f"{ORF6_CTERM_START}-{ORF6_CTERM_END}",
            "c_terminal_sequence": ORF6_CTERM_SEQUENCE,
        },
        "noncoordinate_annotations_used": [
            "ORF6 disrupts bidirectional nucleocytoplasmic transport.",
            "ORF6 interacts via its C-terminus with host RAE1 in the NUP98-RAE1 complex.",
            "ORF6 prevents STAT1 nuclear translocation and antagonizes interferon signaling.",
            "ORF6 localizes to host ER/Golgi membrane contexts and virus-induced vesicular structures.",
            "C-terminal mutagenesis affects NUP98-RAE1 interaction, nuclear transport blockade, and IFN antagonistic function.",
        ],
        "features_excluded_before_prediction": [
            "UniProt PDB cross-references",
            "PDB/mmCIF coordinate files",
            "coordinate-derived contacts",
            "native contact maps",
            "AlphaFold/ESMFold/RoseTTAFold coordinates",
        ],
    }


def _host_hijacking_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "ORF6_RAE1_NUP98_HOST_HIJACKING_TEXT_NONCOORDINATE",
        "source_type": "literature-derived host-interaction and nuclear-transport annotations without coordinates",
        "source_url_or_citation": "Miorin et al. PNAS 2020 PMID 33097660; Addetia et al. mBio 2021 PMID 33849972; Kato et al. BBRC 2021 PMID 33360543",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_prediction": False,
        "answer_key_source": False,
        "labels_used": [
            "C-terminal ORF6 host interaction",
            "RAE1/NUP98 nuclear pore complex context",
            "nuclear import/export disruption",
            "STAT/IRF/interferon antagonism as functional consequence",
            "short-region host-hijacking interface",
        ],
        "coordinate_material_blocked": [
            "PDB IDs and peptide-complex coordinates are not opened before sealing",
            "no contacts, distances, or structural templates are imported",
        ],
    }


def _localization_context_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "ORF6_LOCALIZATION_MEMBRANE_CONTEXT_TEXT_NONCOORDINATE",
        "source_type": "literature-derived localization and membrane-context annotation without coordinate models",
        "source_url_or_citation": "Wong et al. J Cell Sci 2022 PMID 35187564 plus UniProt P0DTC6 localization text",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_prediction": False,
        "answer_key_source": False,
        "labels_used": [
            "host ER/Golgi membrane localization context",
            "localization can be decoupled from IFN antagonism",
            "localization is context, not membrane-channel grammar",
        ],
    }


def _short_region_features_source() -> dict[str, Any]:
    return {
        "source_id": "ORF6_SHORT_REGION_SEQUENCE_FEATURES_DERIVED_FROM_ALLOWED_INPUTS",
        "source_type": "derived non-coordinate sequence, C-terminal motif, and composition features",
        "source_url_or_citation": "derived from UniProt P0DTC6 sequence and allowed text annotations",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_prediction": False,
        "answer_key_source": False,
        "derived_from_source_ids": ["UNIPROT_P0DTC6_SEQUENCE_FUNCTION_NONCOORDINATE"],
        "region_windows": _region_windows(),
        "composition": {
            "full_length": _composition(ORF6_FULL_SEQUENCE),
            "C_terminal_38_61": _composition(ORF6_CTERM_SEQUENCE, ORF6_CTERM_START),
            "minimal_C_terminal_50_61": _composition(ORF6_MINIMAL_CTERM_MOTIF, ORF6_MINIMAL_CTERM_MOTIF_START),
        },
    }


def _problem_context_sources() -> list[dict[str, Any]]:
    return [
        {
            "source_id": "ORF6_UNIPROT_P0DTC6_HOLDOUT_INDEX",
            "source_type": "post-seal UniProt feature/function/mutagenesis holdout index",
            "source_url_or_citation": "https://www.uniprot.org/uniprotkb/P0DTC6/entry",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["Full UniProt C-terminal mutagenesis and interaction records are checked only after sealing."],
        },
        {
            "source_id": "ORF6_MIORIN_NUP98_IFN_HOLDOUT",
            "source_type": "post-seal NUP98/IFN antagonism holdout",
            "source_url_or_citation": "https://europepmc.org/article/MED/33097660",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["ORF6 hijacks Nup98 to block STAT nuclear import and antagonize interferon signaling."],
        },
        {
            "source_id": "ORF6_ADDETIA_RAE1_NUP98_TRANSPORT_HOLDOUT",
            "source_type": "post-seal RAE1/NUP98 transport and C-terminal perturbation holdout",
            "source_url_or_citation": "https://europepmc.org/article/MED/33849972",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["C-terminal ORF6 regions and M58 perturbation affect NUP98-RAE1 interaction and nuclear transport blockade."],
        },
        {
            "source_id": "ORF6_KATO_RAE1_NUP98_DISLOCATION_HOLDOUT",
            "source_type": "post-seal nuclear pore dislocation holdout",
            "source_url_or_citation": "https://europepmc.org/article/MED/33360543",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["ORF6 overexpression dislocates RAE1 and NUP98 from the nuclear pore complex."],
        },
        {
            "source_id": "ORF6_WONG_LOCALIZATION_IFN_DECOUPLING_HOLDOUT",
            "source_type": "post-seal localization/context holdout",
            "source_url_or_citation": "https://europepmc.org/article/MED/35187564",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["ORF6 localization and interferon antagonism can be separated experimentally."],
        },
    ]


def _blocked_prediction_inputs() -> list[dict[str, str]]:
    return [
        {"blocked_input": "PDB coordinates or mmCIF coordinates before sealing", "reason": "coordinates would leak ORF6-host interface truth"},
        {"blocked_input": "PDB-derived contacts or native contact maps", "reason": "V48 tests short-region mechanism grammar, not contact lookup"},
        {"blocked_input": "AlphaFold/ESMFold/RoseTTAFold coordinates before sealing", "reason": "predicted coordinate models are blocked as evidence"},
        {"blocked_input": "ORF6 peptide-complex coordinate templates before sealing", "reason": "host-interface coordinates are validation/structure truth"},
        {"blocked_input": "internal runtime reports as biological evidence", "reason": "runtime artifacts are audit outputs only"},
        {"blocked_input": "validation holdouts before prediction sealing", "reason": "holdouts open only after sealed prediction hash exists"},
        {"blocked_input": "answer key or pass/fail status during prediction", "reason": "classification must use source-separated allowed evidence"},
        {"blocked_input": "target name only", "reason": "ORF6 label alone cannot create a host-hijacking solution packet"},
        {"blocked_input": "generic viral accessory annotation only", "reason": "generic viral identity cannot explain C-terminal RAE1/NUP98 host hijacking"},
        {"blocked_input": "compact stable globular fold forcing", "reason": "V48 must not overpromote ORF6 into a normal single-fold target"},
        {"blocked_input": "swapped ORF8/ORF3a/NSP evidence", "reason": "other viral proteins cannot validate ORF6-specific C-terminal predictions"},
    ]


def build_sources() -> dict[str, Any]:
    for root in [SOURCE_ROOT, PREDICTION_ROOT, HOLDOUT_ROOT, VALIDATION_ROOT, REPORT_ROOT]:
        root.mkdir(parents=True, exist_ok=True)
    prediction_sources = [
        _uniprot_prediction_source(),
        _host_hijacking_prediction_source(),
        _localization_context_prediction_source(),
        _short_region_features_source(),
    ]
    context_sources = _problem_context_sources()
    blocked_inputs = _blocked_prediction_inputs()
    source_manifest = {
        "kind": "V48_SARS2_ORF6_SOURCE_MANIFEST_v0",
        "target_id": "V48_SARS2_ORF6",
        "target": "SARS-CoV-2 ORF6 / accessory protein 6 / ns6",
        "uniprot_accession": ORF6_UNIPROT_ACCESSION,
        "uniprot_id": ORF6_UNIPROT_ID,
        "sequence_region_scope": f"ORF6 {ORF6_UNIPROT_ACCESSION} full length 1-61 with C-terminal host-interface focus {ORF6_CTERM_START}-{ORF6_CTERM_END}",
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
        "live_viral_host_hijacking_solution_packet": False,
    }
    _write_json(SOURCE_ROOT / "source_manifest.json", source_manifest)
    _write_json(SOURCE_ROOT / "uniprot_p0dtc6_prediction_source.json", prediction_sources[0])
    _write_json(SOURCE_ROOT / "orf6_host_hijacking_prediction_source.json", prediction_sources[1])
    _write_json(SOURCE_ROOT / "orf6_localization_context_prediction_source.json", prediction_sources[2])
    _write_json(SOURCE_ROOT / "orf6_short_region_sequence_features.json", prediction_sources[3])
    _write_json(SOURCE_ROOT / "problem_context_sources_not_opened_for_prediction.json", {"sources": context_sources})
    _write_json(SOURCE_ROOT / "blocked_prediction_inputs_manifest.json", {
        "kind": "V48_BLOCKED_PREDICTION_INPUTS_MANIFEST_v0",
        "blocked_prediction_inputs": blocked_inputs,
    })
    _write_json(SOURCE_ROOT / "acquisition_log.json", {
        "kind": "V48_SARS2_ORF6_SOURCE_ACQUISITION_LOG_v0",
        "target_id": "V48_SARS2_ORF6",
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
    full_comp = _composition(ORF6_FULL_SEQUENCE)
    cterm_comp = _composition(ORF6_CTERM_SEQUENCE, ORF6_CTERM_START)
    return {
        "kind": "V48_SARS2_ORF6_SOURCE_BUILD_v0",
        "target_id": "V48_SARS2_ORF6",
        "target": "SARS-CoV-2 ORF6 / accessory protein 6 / ns6",
        "sequence_region_scope": f"ORF6 {ORF6_UNIPROT_ACCESSION} full length 1-61 with C-terminal host-interface focus {ORF6_CTERM_START}-{ORF6_CTERM_END}",
        "sequence_length": len(ORF6_FULL_SEQUENCE),
        "c_terminal_region": f"{ORF6_CTERM_START}-{ORF6_CTERM_END}",
        "c_terminal_sequence": ORF6_CTERM_SEQUENCE,
        "minimal_c_terminal_motif": ORF6_MINIMAL_CTERM_MOTIF,
        "hydrophobic_fraction": full_comp["hydrophobic_fraction"],
        "c_terminal_acidic_fraction": cterm_comp["acidic_fraction"],
        "region_window_count": len(_region_windows()),
        "prediction_input_source_count": len(prediction_sources),
        "holdouts_created": False,
        "source_manifest": str(SOURCE_ROOT / "source_manifest.json"),
        "blocked_inputs_manifest": str(SOURCE_ROOT / "blocked_prediction_inputs_manifest.json"),
        "folding_problem_solved": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V48 SARS-CoV-2 ORF6 viral accessory source manifests.")
    parser.parse_args()
    print(json.dumps(build_sources(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
