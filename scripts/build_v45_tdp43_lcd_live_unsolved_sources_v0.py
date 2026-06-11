#!/usr/bin/env python3
from __future__ import annotations

"""Build V45 TDP-43 LCD live unsolved target source manifests."""

import argparse
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "live_unsolved_targets" / "V45" / "TDP43_LCD"
SOURCE_ROOT = DATA_ROOT / "sources"
PREDICTION_ROOT = DATA_ROOT / "prediction_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts_postseal"
VALIDATION_ROOT = DATA_ROOT / "validation"
REPORT_ROOT = DATA_ROOT / "reports"

TDP43_UNIPROT_ACCESSION = "Q13148"
TDP43_DISPROT_ID = "DP01108"
TDP43_FULL_SEQUENCE = (
    "MSEYIRVTEDENDEPIEIPSEDDGTVLLSTVTAQFPGACGLRYRNPVSQCMRGVRLVEGI"
    "LHAPDAGWGNLVYVVNYPKDNKRKMDETDASSAVKVKRAVQKTSDLIVLGLPWKTTEQDL"
    "KEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRWCDCKLPNS"
    "KQSQDEPLRSRKVFVGRCTEDMTEDELREFFSQYGDVMDVFIPKPFRAFAFVTFADDQIA"
    "QSLCGEDLIIKGISVHISNAEPKHNSNRQLERSGRFGGNPGGFGNQGGFGNSRGGGAGLG"
    "NNQGSNMGGGMNFGAFSINPAMMAAAQAALQSSWGMMGMLASQQNQSGPSGNNQNQGNMQ"
    "REPNQAFGSGNNSYSGSNSGAAIGWGSASNAGSGSGFNGGFGSSMDSKSSGWGM"
)
TDP43_LCD_SEQUENCE = TDP43_FULL_SEQUENCE[273:414]


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _positions(sequence: str, residues: str, offset: int) -> list[dict[str, Any]]:
    return [
        {"position": offset + index, "residue": aa}
        for index, aa in enumerate(sequence)
        if aa in residues
    ]


def _composition(sequence: str, offset: int = 274) -> dict[str, Any]:
    counts = Counter(sequence)
    length = len(sequence)
    low_complexity_residues = "GNQSYFM"
    return {
        "length": length,
        "counts": {aa: counts[aa] for aa in sorted(counts)},
        "fractions": {aa: counts[aa] / length for aa in sorted(counts)},
        "gnqsyfm_count": sum(counts[aa] for aa in low_complexity_residues),
        "gnqsyfm_fraction": sum(counts[aa] for aa in low_complexity_residues) / length,
        "glycine_fraction": counts["G"] / length,
        "qn_fraction": (counts["Q"] + counts["N"]) / length,
        "aromatic_positions": _positions(sequence, "FYW", offset),
        "methionine_positions": _positions(sequence, "M", offset),
        "serine_positions": _positions(sequence, "S", offset),
    }


def _sequence_windows(sequence: str) -> list[dict[str, Any]]:
    windows = []
    for start, end in [(274, 303), (304, 343), (321, 343), (344, 366), (367, 414), (403, 414)]:
        sub = TDP43_FULL_SEQUENCE[start - 1:end]
        comp = _composition(sub, start)
        windows.append({
            "span": f"{start}-{end}",
            "length": len(sub),
            "glycine_fraction": comp["glycine_fraction"],
            "qn_fraction": comp["qn_fraction"],
            "aromatic_count": len(comp["aromatic_positions"]),
            "methionine_count": len(comp["methionine_positions"]),
            "serine_count": len(comp["serine_positions"]),
            "sequence_pattern_note": "sequence-derived TDP-43 CTD/LCD operator window",
        })
    return windows


def _uniprot_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "UNIPROT_Q13148_SEQUENCE_FEATURES_NONCOORDINATE",
        "source_type": "UniProtKB sequence and non-coordinate feature annotation",
        "source_url_or_citation": "https://rest.uniprot.org/uniprotkb/Q13148.json",
        "accession": TDP43_UNIPROT_ACCESSION,
        "retrieved_or_checked": "2026-06-11",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "answer_key_source": False,
        "features_used": [
            {"type": "Region", "description": "Disordered", "start": 261, "end": 303},
            {"type": "Region", "description": "Disordered", "start": 341, "end": 373},
            {"type": "Compositional bias", "description": "Gly residues", "start": 275, "end": 303},
            {"type": "Compositional bias", "description": "Low complexity", "start": 342, "end": 358},
            {"type": "Region", "description": "Interaction with UBQLN2", "start": 216, "end": 414},
            {"type": "Modified residue", "description": "Phosphoserine", "start": 292, "end": 292},
            {"type": "Modified residue", "description": "Omega-N-methylarginine", "start": 293, "end": 293},
        ],
        "features_excluded_before_prediction": [
            "PDB/mmCIF coordinate records",
            "coordinate-derived contacts",
            "native contact maps",
            "UniProt beta strand/helix records when coordinate-derived",
            "AlphaFold/ESMFold/RoseTTAFold coordinates",
        ],
        "sequence_region_scope": {
            "protein": "human TDP-43 / TAR DNA-binding protein 43",
            "region": "TDP-43 LCD/CTD residues 274-414",
            "full_length": len(TDP43_FULL_SEQUENCE),
            "region_start": 274,
            "region_end": 414,
            "sequence": TDP43_LCD_SEQUENCE,
        },
    }


def _disprot_prediction_source() -> dict[str, Any]:
    return {
        "source_id": "DISPROT_DP01108_DISORDER_ONLY_NONCOORDINATE",
        "source_type": "DisProt disorder annotation subset",
        "source_url_or_citation": "https://disprot.org/api/Q13148",
        "accession": TDP43_UNIPROT_ACCESSION,
        "disprot_id": TDP43_DISPROT_ID,
        "retrieved_or_checked": "2026-06-11",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "answer_key_source": False,
        "regions_used": [
            {"region_id": "DP01108r001", "term": "disorder", "namespace": "Structural state", "start": 266, "end": 414, "reference_id": "29511089"},
            {"region_id": "DP01108r005", "term": "disorder", "namespace": "Structural state", "start": 267, "end": 414, "reference_id": "27545621"},
            {"region_id": "DP01108r013", "term": "disorder", "namespace": "Structural state", "start": 263, "end": 414, "reference_id": "26735904"},
        ],
        "regions_excluded_until_holdout": [
            "molecular condensate scaffold activity",
            "amyloid fibril formation",
            "disorder-to-order transition",
            "DNA/RNA/nucleic-acid binding context",
            "alpha-helical LLPS segment literature",
            "aromatic-residue mutational evidence",
        ],
    }


def _pattern_prediction_source() -> dict[str, Any]:
    composition = _composition(TDP43_LCD_SEQUENCE)
    return {
        "source_id": "TDP43_LCD_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS",
        "source_type": "derived sequence-composition and sparse-sticker features",
        "source_url_or_citation": "derived from UniProt Q13148 residues 274-414 sequence",
        "allowed_use": "prediction_input_before_sealing",
        "holdout_source": False,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "answer_key_source": False,
        "derived_from_source_ids": ["UNIPROT_Q13148_SEQUENCE_FEATURES_NONCOORDINATE"],
        "composition": composition,
        "sequence_windows": _sequence_windows(TDP43_LCD_SEQUENCE),
    }


def _problem_context_sources() -> list[dict[str, Any]]:
    return [
        {
            "source_id": "DISPROT_DP01108_FULL_REGION_HOLDOUT_INDEX",
            "source_type": "DisProt post-seal region/function holdout index",
            "source_url_or_citation": "https://disprot.org/api/Q13148",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": [
                "Full DisProt region/function annotations are kept out of prediction packet construction.",
                "Only disorder-only records are available before sealing.",
            ],
        },
        {
            "source_id": "TDP43_AROMATIC_LLPS_PUBMED_29511089",
            "source_type": "post-seal literature context",
            "source_url_or_citation": "https://pubmed.ncbi.nlm.nih.gov/29511089/",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["TDP-43 CTD LLPS is mediated by a small number of aromatic residues."],
        },
        {
            "source_id": "TDP43_ALPHA_HELIX_LLPS_PUBMED_27545621",
            "source_type": "post-seal literature context",
            "source_url_or_citation": "https://pubmed.ncbi.nlm.nih.gov/27545621/",
            "allowed_use": "postseal_holdout_validation_only",
            "holdout_source": True,
            "coordinate_derived": False,
            "internal_runtime_source": False,
            "notes": ["ALS mutations disrupt phase separation mediated by an alpha-helical structure in the low-complexity CTD."],
        },
    ]


def _blocked_prediction_inputs() -> list[dict[str, str]]:
    return [
        {"blocked_input": "PDB coordinates or mmCIF coordinates before sealing", "reason": "coordinates are validation/structure truth and not prediction evidence"},
        {"blocked_input": "coordinate-derived contacts or native contact maps", "reason": "would leak structure truth into the live solution packet"},
        {"blocked_input": "AlphaFold/ESMFold/RoseTTAFold coordinates before sealing", "reason": "predicted coordinate models are blocked as evidence"},
        {"blocked_input": "first_contact_clean_pharmacotopology_layer_run files as biological evidence", "reason": "runtime artifacts are audit outputs only"},
        {"blocked_input": "validation holdouts before prediction sealing", "reason": "holdouts open only after sealed prediction hash exists"},
        {"blocked_input": "answer key/class labels during prediction", "reason": "class assignment must use source-separated allowed evidence"},
        {"blocked_input": "FUS-LC answer transfer", "reason": "TDP-43 LCD must be solved as a distinct sparse-aromatic/prion-like grammar"},
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
        "kind": "V45_TDP43_LCD_LIVE_UNSOLVED_SOURCE_MANIFEST_v0",
        "target_id": "V45_TDP43_LCD",
        "target": "human TDP-43 low-complexity C-terminal domain",
        "uniprot_accession": TDP43_UNIPROT_ACCESSION,
        "disprot_id": TDP43_DISPROT_ID,
        "sequence_region_scope": "TDP-43 Q13148 residues 274-414",
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
    _write_json(SOURCE_ROOT / "uniprot_q13148_prediction_source.json", prediction_sources[0])
    _write_json(SOURCE_ROOT / "disprot_dp01108_disorder_prediction_source.json", prediction_sources[1])
    _write_json(SOURCE_ROOT / "tdp43_lcd_sequence_pattern_features.json", prediction_sources[2])
    _write_json(SOURCE_ROOT / "problem_context_sources_not_opened_for_prediction.json", {"sources": context_sources})
    _write_json(SOURCE_ROOT / "blocked_prediction_inputs_manifest.json", {
        "kind": "V45_BLOCKED_PREDICTION_INPUTS_MANIFEST_v0",
        "blocked_prediction_inputs": blocked_inputs,
    })
    _write_json(SOURCE_ROOT / "acquisition_log.json", {
        "kind": "V45_TDP43_LCD_SOURCE_ACQUISITION_LOG_v0",
        "target_id": "V45_TDP43_LCD",
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
    composition = _composition(TDP43_LCD_SEQUENCE)
    return {
        "kind": "V45_TDP43_LCD_SOURCE_BUILD_v0",
        "target_id": "V45_TDP43_LCD",
        "target": "human TDP-43 low-complexity C-terminal domain",
        "sequence_region_scope": "TDP-43 Q13148 residues 274-414",
        "sequence_length": len(TDP43_LCD_SEQUENCE),
        "gnqsyfm_fraction": composition["gnqsyfm_fraction"],
        "glycine_fraction": composition["glycine_fraction"],
        "aromatic_count": len(composition["aromatic_positions"]),
        "methionine_count": len(composition["methionine_positions"]),
        "serine_count": len(composition["serine_positions"]),
        "prediction_input_source_count": len(prediction_sources),
        "holdouts_created": False,
        "source_manifest": str(SOURCE_ROOT / "source_manifest.json"),
        "blocked_inputs_manifest": str(SOURCE_ROOT / "blocked_prediction_inputs_manifest.json"),
        "folding_problem_solved": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V45 TDP-43 LCD live unsolved target sources.")
    parser.parse_args()
    print(json.dumps(build_sources(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
