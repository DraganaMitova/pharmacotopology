#!/usr/bin/env python3
from __future__ import annotations

"""Build V46 CFTR F508del membrane multidomain source manifests."""

import argparse
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "live_unsolved_targets" / "V46" / "CFTR_F508del"
SOURCE_ROOT = DATA_ROOT / "sources"
PREDICTION_ROOT = DATA_ROOT / "prediction_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts_postseal"
VALIDATION_ROOT = DATA_ROOT / "validation"
REPORT_ROOT = DATA_ROOT / "reports"

CFTR_UNIPROT_ACCESSION = "P13569"
CFTR_FULL_SEQUENCE = (
    "MQRSPLEKASVVSKLFFSWTRPILRKGYRQRLELSDIYQIPSVDSADNLSEKLEREWDRELASKKNPKLINALRRCFFWRFMFYGIFLYLGEVTKAVQPL"
    "LLGRIIASYDPDNKEERSIAIYLGIGLCLLFIVRTLLLHPAIFGLHHIGMQMRIAMFSLIYKKTLKLSSRVLDKISIGQLVSLLSNNLNKFDEGLALAHF"
    "VWIAPLQVALLMGLIWELLQASAFCGLGFLIVLALFQAGLGRMMMKYRDQRAGKISERLVITSEMIENIQSVKAYCWEEAMEKMIENLRQTELKLTRKAA"
    "YVRYFNSSAFFFSGFFVVFLSVLPYALIKGIILRKIFTTISFCIVLRMAVTRQFPWAVQTWYDSLGAINKIQDFLQKQEYKTLEYNLTTTEVVMENVTAF"
    "WEEGFGELFEKAKQNNNNRKTSNGDDSLFFSNFSLLGTPVLKDINFKIERGQLLAVAGSTGAGKTSLLMVIMGELEPSEGKIKHSGRISFCSQFSWIMPG"
    "TIKENIIFGVSYDEYRYRSVIKACQLEEDISKFAEKDNIVLGEGGITLSGGQRARISLARAVYKDADLYLLDSPFGYLDVLTEKEIFESCVCKLMANKTR"
    "ILVTSKMEHLKKADKILILHEGSSYFYGTFSELQNLQPDFSSKLMGCDSFDQFSAERRNSILTETLHRFSLEGDAPVSWTETKKQSFKQTGEFGEKRKNS"
    "ILNPINSIRKFSIVQKTPLQMNGIEEDSDEPLERRLSLVPDSEQGEAILPRISVISTGPTLQARRRQSVLNLMTHSVNQGQNIHRKTTASTRKVSLAPQA"
    "NLTELDIYSRRLSQETGLEISEEINEEDLKECFFDDMESIPAVTTWNTYLRYITVHKSLIFVLIWCLVIFLAEVAASLVVLWLLGNTPLQDKGNSTHSRN"
    "NSYAVIITSTSSYYVFYIYVGVADTLLAMGFFRGLPLVHTLITVSKILHHKMLHSVLQAPMSTLNTLKAGGILNRFSKDIAILDDLLPLTIFDFIQLLLI"
    "VIGAIAVVAVLQPYIFVATVPVIVAFIMLRAYFLQTSQQLKQLESEGRSPIFTHLVTSLKGLWTLRAFGRQPYFETLFHKALNLHTANWFLYLSTLRWFQ"
    "MRIEMIFVIFFIAVTFISILTTGEGEGRVGIILTLAMNIMSTLQWAVNSSIDVDSLMRSVSRVFKFIDMPTEGKPTKSTKPYKNGQLSKVMIIENSHVKK"
    "DDIWPSGGQMTVKDLTAKYTEGGNAILENISFSISPGQRVGLLGRTGSGKSTLLSAFLRLLNTEGEIQIDGVSWDSITLQQWRKAFGVIPQKVFIFSGTF"
    "RKNLDPYEQWSDQEIWKVADEVGLRSVIEQFPGKLDFVLVDGGCVLSHGHKQLMCLARSVLSKAKILLLDEPSAHLDPVTYQIIRRTLKQAFADCTVILC"
    "EHRIEAMLECQQFLVIEENKVRQYDSIQKLLNERSLFRQAISPSDRVKLFPHRNSSKCKSKPQIAALKEETEEEVQDTRL"
)


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _composition(sequence: str) -> dict[str, Any]:
    counts = Counter(sequence)
    length = len(sequence)
    hydrophobic = "AILMFWVY"
    charged = "DEKRH"
    return {
        "length": length,
        "counts": {aa: counts[aa] for aa in sorted(counts)},
        "hydrophobic_fraction": sum(counts[aa] for aa in hydrophobic) / length,
        "charged_fraction": sum(counts[aa] for aa in charged) / length,
    }


def _domain_architecture() -> list[dict[str, Any]]:
    return [
        {"domain": "MSD1", "feature_type": "ABC transmembrane type-1 1", "start": 81, "end": 365},
        {"domain": "NBD1", "feature_type": "ABC transporter 1", "start": 423, "end": 646},
        {"domain": "R_region", "feature_type": "Disordered regulatory region", "start": 654, "end": 831},
        {"domain": "MSD2", "feature_type": "ABC transmembrane type-1 2", "start": 859, "end": 1155},
        {"domain": "NBD2", "feature_type": "ABC transporter 2", "start": 1210, "end": 1443},
    ]


def _transmembrane_segments() -> list[dict[str, Any]]:
    return [
        {"segment": 1, "start": 78, "end": 98},
        {"segment": 2, "start": 123, "end": 146},
        {"segment": 3, "start": 196, "end": 216},
        {"segment": 4, "start": 223, "end": 243},
        {"segment": 5, "start": 299, "end": 319},
        {"segment": 6, "start": 340, "end": 358},
        {"segment": 7, "start": 859, "end": 879},
        {"segment": 8, "start": 919, "end": 939},
        {"segment": 9, "start": 991, "end": 1011},
        {"segment": 10, "start": 1014, "end": 1034},
        {"segment": 11, "start": 1096, "end": 1116},
        {"segment": 12, "start": 1131, "end": 1151},
    ]


def _uniprot_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "UNIPROT_P13569_SEQUENCE_FEATURES_NONCOORDINATE",
        "source_type": "UniProtKB sequence, topology, domain, function, and non-coordinate feature annotation",
        "source_url_or_citation": "https://rest.uniprot.org/uniprotkb/P13569.json",
        "accession": CFTR_UNIPROT_ACCESSION,
        "retrieved_or_checked": "2026-06-11",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "answer_key_source": False,
        "sequence_region_scope": {
            "protein": "human CFTR / cystic fibrosis transmembrane conductance regulator",
            "region": "full-length CFTR with F508del focus in NBD1",
            "full_length": len(CFTR_FULL_SEQUENCE),
            "f508_position": 508,
            "f508_context": CFTR_FULL_SEQUENCE[487:528],
            "sequence": CFTR_FULL_SEQUENCE,
        },
        "domain_architecture_used": _domain_architecture(),
        "transmembrane_segments_used": _transmembrane_segments(),
        "function_annotations_used": [
            "epithelial chloride channel / ABC transporter family context",
            "ATP-dependent channel regulation and R-region phosphorylation context",
            "cell-surface trafficking and membrane localization context",
        ],
        "features_excluded_before_prediction": [
            "PDB/mmCIF coordinate records",
            "coordinate-derived contacts",
            "native contact maps",
            "AlphaFold/ESMFold/RoseTTAFold coordinates",
        ],
    }


def _f508del_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "UNIPROT_P13569_F508DEL_VARIANT_NONCOORDINATE",
        "source_type": "UniProtKB disease-variant annotation",
        "source_url_or_citation": "https://rest.uniprot.org/uniprotkb/P13569.json",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "answer_key_source": False,
        "variant": {
            "label": "F508del / deletion of phenylalanine 508",
            "position": 508,
            "wild_type_residue": "F",
            "domain_context": "NBD1 / ABC transporter 1 domain",
            "noncoordinate_annotation": [
                "impairs protein folding and stability",
                "causes local changes to a surface mediating interactions between domains",
                "impairs maturation and trafficking to the cell membrane",
                "primes degradation through a proteasome-linked pathway",
                "decreases frequency of channel opening in vitro",
            ],
        },
    }


def _domain_pattern_source() -> dict[str, Any]:
    nbd1 = CFTR_FULL_SEQUENCE[422:646]
    msd1 = CFTR_FULL_SEQUENCE[80:365]
    r_region = CFTR_FULL_SEQUENCE[653:831]
    return {
        "source_id": "CFTR_DOMAIN_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS",
        "source_type": "derived non-coordinate domain/topology composition features",
        "source_url_or_citation": "derived from UniProt P13569 sequence and domain annotations",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "answer_key_source": False,
        "derived_from_source_ids": ["UNIPROT_P13569_SEQUENCE_FEATURES_NONCOORDINATE"],
        "domain_architecture": _domain_architecture(),
        "transmembrane_segment_count": len(_transmembrane_segments()),
        "composition": {
            "full_length": _composition(CFTR_FULL_SEQUENCE),
            "MSD1_81_365": _composition(msd1),
            "NBD1_423_646": _composition(nbd1),
            "R_region_654_831": _composition(r_region),
        },
    }


def _problem_context_sources() -> list[dict[str, Any]]:
    return [
        {
            "source_id": "CFTR_F508DEL_CORRECTION_REQUIREMENTS_PMC3266553",
            "source_type": "post-seal review/experimental literature context",
            "source_url_or_citation": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3266553/",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": [
                "Efficient correction of F508del-CFTR requires more than a local mutation fix.",
                "NBD1 stability and domain/interface effects are part of the rescue grammar.",
            ],
        },
        {
            "source_id": "CFTR_F508DEL_NBD1_STABILITY_LITERATURE_HOLDOUT",
            "source_type": "post-seal literature context",
            "source_url_or_citation": "PubMed literature on F508del NBD1 folding/stability defects",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["F508del lowers NBD1 folding competence/stability and causes maturation defects."],
        },
        {
            "source_id": "CFTR_F508DEL_CORRECTOR_RESCUE_CONTEXT_HOLDOUT",
            "source_type": "post-seal corrector/rescue literature context",
            "source_url_or_citation": "CFTR corrector/rescue literature; VX-809/lumacaftor and combination corrector logic",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["Corrector action maps to multidomain folding, interface, and proteostasis rescue rather than one-site repair."],
        },
    ]


def _blocked_prediction_inputs() -> list[dict[str, str]]:
    return [
        {"blocked_input": "PDB coordinates or mmCIF coordinates before sealing", "reason": "coordinates are validation/structure truth and not prediction evidence"},
        {"blocked_input": "PDB-derived contacts or native contact maps", "reason": "would leak structure/interface truth into prediction"},
        {"blocked_input": "AlphaFold/ESMFold/RoseTTAFold coordinates before sealing", "reason": "predicted coordinate models are blocked as evidence"},
        {"blocked_input": "CFTR coordinate models as prediction evidence", "reason": "V46 tests mechanism grammar, not coordinate lookup"},
        {"blocked_input": "internal runtime reports as biological evidence", "reason": "runtime artifacts are audit outputs only"},
        {"blocked_input": "validation holdouts before prediction sealing", "reason": "holdouts open only after sealed prediction hash exists"},
        {"blocked_input": "answer key/class labels during prediction", "reason": "assignment must use source-separated allowed evidence"},
        {"blocked_input": "generic channel annotation only", "reason": "channel identity alone cannot explain F508del folding rescue"},
        {"blocked_input": "single local F508 deletion explanation only", "reason": "V46 must model NBD1 plus interface plus proteostasis/maturation grammar"},
    ]


def build_sources() -> dict[str, Any]:
    for root in [SOURCE_ROOT, PREDICTION_ROOT, HOLDOUT_ROOT, VALIDATION_ROOT, REPORT_ROOT]:
        root.mkdir(parents=True, exist_ok=True)
    prediction_sources = [
        _uniprot_prediction_source(),
        _f508del_prediction_source(),
        _domain_pattern_source(),
    ]
    context_sources = _problem_context_sources()
    blocked_inputs = _blocked_prediction_inputs()
    source_manifest = {
        "kind": "V46_CFTR_F508DEL_SOURCE_MANIFEST_v0",
        "target_id": "V46_CFTR_F508DEL",
        "target": "human CFTR F508del membrane multidomain folding defect",
        "uniprot_accession": CFTR_UNIPROT_ACCESSION,
        "sequence_region_scope": "full-length CFTR P13569 with F508del focus in NBD1",
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
        "live_membrane_solution_packet": False,
    }
    _write_json(SOURCE_ROOT / "source_manifest.json", source_manifest)
    _write_json(SOURCE_ROOT / "uniprot_p13569_prediction_source.json", prediction_sources[0])
    _write_json(SOURCE_ROOT / "f508del_variant_prediction_source.json", prediction_sources[1])
    _write_json(SOURCE_ROOT / "cftr_domain_pattern_features.json", prediction_sources[2])
    _write_json(SOURCE_ROOT / "problem_context_sources_not_opened_for_prediction.json", {"sources": context_sources})
    _write_json(SOURCE_ROOT / "blocked_prediction_inputs_manifest.json", {
        "kind": "V46_BLOCKED_PREDICTION_INPUTS_MANIFEST_v0",
        "blocked_prediction_inputs": blocked_inputs,
    })
    _write_json(SOURCE_ROOT / "acquisition_log.json", {
        "kind": "V46_CFTR_F508DEL_SOURCE_ACQUISITION_LOG_v0",
        "target_id": "V46_CFTR_F508DEL",
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
    return {
        "kind": "V46_CFTR_F508DEL_SOURCE_BUILD_v0",
        "target_id": "V46_CFTR_F508DEL",
        "target": "human CFTR F508del membrane multidomain folding defect",
        "sequence_region_scope": "full-length CFTR P13569 with F508del focus in NBD1",
        "sequence_length": len(CFTR_FULL_SEQUENCE),
        "f508_residue": CFTR_FULL_SEQUENCE[507],
        "domain_count": len(_domain_architecture()),
        "transmembrane_segment_count": len(_transmembrane_segments()),
        "prediction_input_source_count": len(prediction_sources),
        "holdouts_created": False,
        "source_manifest": str(SOURCE_ROOT / "source_manifest.json"),
        "blocked_inputs_manifest": str(SOURCE_ROOT / "blocked_prediction_inputs_manifest.json"),
        "folding_problem_solved": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V46 CFTR F508del membrane multidomain source manifests.")
    parser.parse_args()
    print(json.dumps(build_sources(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
