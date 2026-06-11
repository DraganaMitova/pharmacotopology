#!/usr/bin/env python3
from __future__ import annotations

"""Build V44 FUS-LC live unsolved target source manifests."""

import argparse
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "live_unsolved_targets" / "V44" / "FUS_LC"
SOURCE_ROOT = DATA_ROOT / "sources"
PREDICTION_ROOT = DATA_ROOT / "prediction_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts_postseal"
VALIDATION_ROOT = DATA_ROOT / "validation"
REPORT_ROOT = DATA_ROOT / "reports"

FUS_UNIPROT_ACCESSION = "P35637"
FUS_DISPROT_ID = "DP01102"
FUS_FULL_SEQUENCE = (
    "MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQ"
    "SQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSS"
    "SYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQD"
    "QSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYE"
    "PRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSDNNTIFVQGLGENVTIES"
    "VADYFKQIGIIKTNKKTGQPMINLYTDRETGKLKGEATVSFDDPPSAKAAIDWFDGKEFS"
    "GNPIKVSFATRRADFNRGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGGQ"
    "QRAGDWKCPNPTCENMNFSWRNECNQCKAPKPDGPGGGPGGSHMGGNYGDDRRGGRGGYD"
    "RGGYRGRGGDRGGFRGGRGGGDRGGFGPGKMDSRGEHRQDRRERPY"
)
FUS_LC_SEQUENCE = FUS_FULL_SEQUENCE[:214]


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _positions(sequence: str, residue: str) -> list[int]:
    return [index + 1 for index, aa in enumerate(sequence) if aa == residue]


def _composition(sequence: str) -> dict[str, Any]:
    counts = Counter(sequence)
    length = len(sequence)
    return {
        "length": length,
        "counts": {aa: counts[aa] for aa in sorted(counts)},
        "fractions": {aa: counts[aa] / length for aa in sorted(counts)},
        "qgsy_count": sum(counts[aa] for aa in "QGSY"),
        "qgsy_fraction": sum(counts[aa] for aa in "QGSY") / length,
        "tyrosine_positions": _positions(sequence, "Y"),
        "serine_positions": _positions(sequence, "S"),
    }


def _sequence_windows(sequence: str) -> list[dict[str, Any]]:
    windows = []
    for start, end in [(1, 75), (76, 112), (113, 164), (165, 214), (39, 95), (110, 150)]:
        sub = sequence[start - 1:end]
        comp = _composition(sub)
        windows.append({
            "span": f"{start}-{end}",
            "length": len(sub),
            "qgsy_fraction": comp["qgsy_fraction"],
            "tyrosine_count": comp["counts"].get("Y", 0),
            "serine_count": comp["counts"].get("S", 0),
            "glutamine_count": comp["counts"].get("Q", 0),
            "glycine_count": comp["counts"].get("G", 0),
            "sequence_pattern_note": "sequence-derived sticker/spacer and low-complexity window",
        })
    return windows


def _uniprot_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "UNIPROT_P35637_SEQUENCE_FEATURES_NONCOORDINATE",
        "source_type": "UniProtKB sequence and non-coordinate feature annotation",
        "source_url_or_citation": "https://rest.uniprot.org/uniprotkb/P35637.json",
        "accession": FUS_UNIPROT_ACCESSION,
        "retrieved_or_checked": "2026-06-11",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "answer_key_source": False,
        "features_used": [
            {"type": "Region", "description": "Disordered", "start": 1, "end": 286},
            {"type": "Compositional bias", "description": "Polar residues", "start": 1, "end": 14},
            {"type": "Compositional bias", "description": "Low complexity", "start": 17, "end": 75},
            {"type": "Compositional bias", "description": "Low complexity", "start": 83, "end": 164},
            {"type": "Compositional bias", "description": "Gly residues", "start": 165, "end": 177},
            {"type": "Compositional bias", "description": "Gly residues", "start": 186, "end": 209},
            {"type": "Modified residue", "description": "Phosphoserine", "start": 26, "end": 26},
            {"type": "Modified residue", "description": "Phosphoserine", "start": 30, "end": 30},
            {"type": "Modified residue", "description": "Phosphoserine; by ATM", "start": 42, "end": 42},
        ],
        "features_excluded_before_prediction": [
            "PDB/mmCIF coordinate records",
            "coordinate-derived contacts",
            "native contact maps",
            "UniProt beta strand/turn records when coordinate-derived",
            "AlphaFold/ESMFold/RoseTTAFold coordinates",
        ],
        "sequence_region_scope": {
            "protein": "human FUS / RNA-binding protein FUS",
            "region": "FUS-LC residues 1-214",
            "full_length": len(FUS_FULL_SEQUENCE),
            "region_start": 1,
            "region_end": 214,
            "sequence": FUS_LC_SEQUENCE,
        },
    }


def _disprot_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "DISPROT_DP01102_DISORDER_ONLY_NONCOORDINATE",
        "source_type": "DisProt disorder annotation subset",
        "source_url_or_citation": "https://disprot.org/api/P35637",
        "accession": FUS_UNIPROT_ACCESSION,
        "disprot_id": FUS_DISPROT_ID,
        "retrieved_or_checked": "2026-06-11",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "answer_key_source": False,
        "regions_used": [
            {
                "region_id": "DP01102r004",
                "term": "disorder",
                "namespace": "Structural state",
                "start": 1,
                "end": 507,
                "reference_id": "29677513",
            },
            {
                "region_id": "DP01102r014",
                "term": "disorder",
                "namespace": "Structural state",
                "start": 1,
                "end": 163,
                "reference_id": "26455390",
            },
            {
                "region_id": "DP01102r029",
                "term": "disorder",
                "namespace": "Structural state",
                "start": 98,
                "end": 214,
                "reference_id": "28942918",
            },
        ],
        "regions_excluded_until_holdout": [
            "amyloid fibril formation regions",
            "RNA binding regions",
            "molecular condensate scaffold activity regions",
            "protein-partner phase-separation inhibition regions",
        ],
    }


def _pattern_prediction_source() -> dict[str, Any]:
    composition = _composition(FUS_LC_SEQUENCE)
    return {
        "source_id": "FUS_LC_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS",
        "source_type": "derived sequence-composition and sticker-spacer features",
        "source_url_or_citation": "derived from UniProt P35637 residues 1-214 sequence",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "answer_key_source": False,
        "derived_from_source_ids": ["UNIPROT_P35637_SEQUENCE_FEATURES_NONCOORDINATE"],
        "composition": composition,
        "sequence_windows": _sequence_windows(FUS_LC_SEQUENCE),
    }


def _problem_context_sources() -> list[dict[str, Any]]:
    return [
        {
            "source_id": "FUS_LC_REVIEW_PROBLEM_CONTEXT_ARXIV_2303_04215",
            "source_type": "review/problem-context literature",
            "source_url_or_citation": "https://arxiv.org/abs/2303.04215",
            "allowed_use": "problem_context_and_postseal_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": [
                "FUS-LC residues 1-214 are condition-sensitive IDP/condensate material.",
                "Gels, fibrils, and glass-like behavior vary by experimental condition.",
                "Phosphorylation effects are site-specific and alter state stability.",
            ],
        },
        {
            "source_id": "DISPROT_DP01102_FULL_REGION_HOLDOUT_INDEX",
            "source_type": "DisProt post-seal region/function holdout index",
            "source_url_or_citation": "https://disprot.org/api/P35637",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": [
                "Full DisProt region/function annotations are kept out of prediction packet construction.",
                "Only disorder-only records are available before sealing.",
            ],
        },
    ]


def _blocked_prediction_inputs() -> list[dict[str, str]]:
    return [
        {"blocked_input": "PDB coordinates or mmCIF coordinates before sealing", "reason": "coordinates are validation/structure truth and not prediction evidence"},
        {"blocked_input": "coordinate-derived contacts or native contact maps", "reason": "would leak structure truth into the live solution packet"},
        {"blocked_input": "AlphaFold/ESMFold/RoseTTAFold coordinates before sealing", "reason": "predicted coordinate models are blocked as evidence for this live mechanism packet"},
        {"blocked_input": "first_contact_clean_pharmacotopology_layer_run files as biological evidence", "reason": "runtime artifacts are audit outputs only"},
        {"blocked_input": "validation holdouts before prediction sealing", "reason": "holdouts open only after sealed prediction hash exists"},
        {"blocked_input": "answer key/class labels during prediction", "reason": "class assignment must use source-separated allowed evidence"},
        {"blocked_input": "generic target-name-only assignment", "reason": "FUS name alone cannot create a solution packet"},
    ]


def build_sources() -> dict[str, Any]:
    for root in [SOURCE_ROOT, PREDICTION_ROOT, HOLDOUT_ROOT, VALIDATION_ROOT, REPORT_ROOT]:
        root.mkdir(parents=True, exist_ok=True)
    prediction_sources = [
        _uniprot_prediction_source(),
        _disprot_prediction_source(),
        _pattern_prediction_source(),
    ]
    context_sources = _problem_context_sources()
    blocked_inputs = _blocked_prediction_inputs()
    source_manifest = {
        "kind": "V44_FUS_LC_LIVE_UNSOLVED_SOURCE_MANIFEST_v0",
        "target_id": "V44_FUS_LC",
        "target": "human FUS low-complexity domain",
        "uniprot_accession": FUS_UNIPROT_ACCESSION,
        "disprot_id": FUS_DISPROT_ID,
        "sequence_region_scope": "FUS P35637 residues 1-214",
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
        "live_unsolved_target_solution_packet": False,
    }
    _write_json(SOURCE_ROOT / "source_manifest.json", source_manifest)
    _write_json(SOURCE_ROOT / "uniprot_p35637_prediction_source.json", prediction_sources[0])
    _write_json(SOURCE_ROOT / "disprot_dp01102_disorder_prediction_source.json", prediction_sources[1])
    _write_json(SOURCE_ROOT / "fus_lc_sequence_pattern_features.json", prediction_sources[2])
    _write_json(SOURCE_ROOT / "problem_context_sources_not_opened_for_prediction.json", {"sources": context_sources})
    _write_json(SOURCE_ROOT / "blocked_prediction_inputs_manifest.json", {
        "kind": "V44_BLOCKED_PREDICTION_INPUTS_MANIFEST_v0",
        "blocked_prediction_inputs": blocked_inputs,
    })
    _write_json(SOURCE_ROOT / "acquisition_log.json", {
        "kind": "V44_FUS_LC_SOURCE_ACQUISITION_LOG_v0",
        "target_id": "V44_FUS_LC",
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
    composition = _composition(FUS_LC_SEQUENCE)
    return {
        "kind": "V44_FUS_LC_SOURCE_BUILD_v0",
        "target_id": "V44_FUS_LC",
        "target": "human FUS low-complexity domain",
        "sequence_region_scope": "FUS P35637 residues 1-214",
        "sequence_length": len(FUS_LC_SEQUENCE),
        "qgsy_fraction": composition["qgsy_fraction"],
        "tyrosine_count": composition["counts"].get("Y", 0),
        "serine_count": composition["counts"].get("S", 0),
        "prediction_input_source_count": len(prediction_sources),
        "holdouts_created": False,
        "source_manifest": str(SOURCE_ROOT / "source_manifest.json"),
        "blocked_inputs_manifest": str(SOURCE_ROOT / "blocked_prediction_inputs_manifest.json"),
        "folding_problem_solved": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V44 FUS-LC live unsolved target sources.")
    parser.parse_args()
    print(json.dumps(build_sources(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
