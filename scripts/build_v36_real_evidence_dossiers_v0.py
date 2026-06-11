#!/usr/bin/env python3
from __future__ import annotations

"""Build V36 real external evidence dossiers.

V36 is acquisition/dossier work, not another prediction gate. It records small,
auditable, non-coordinate evidence packages for three hard protein classes:
KcsA membrane-channel grammar, XCL1 metamorphic state-switch grammar, and
alpha-synuclein disorder/ensemble grammar.
"""

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "external_evidence_dossiers"

TARGET_ORDER = ["KcsA", "XCL1_lymphotactin", "alpha_synuclein_SNCA"]


def _row(
    target: str,
    source_id: str,
    evidence_bucket: str,
    source_type: str,
    source_name: str,
    source_url_or_citation: str,
    source_date_or_version: str,
    evidence_statement: str,
    grammar_role: str,
) -> dict[str, Any]:
    return {
        "target": target,
        "source_id": source_id,
        "evidence_bucket": evidence_bucket,
        "source_type": source_type,
        "source_name": source_name,
        "source_url_or_citation": source_url_or_citation,
        "source_date_or_version": source_date_or_version,
        "source_boundary": "external_noncoordinate_scientific_annotation_or_state_evidence",
        "allowed_use": "V36_claim_disabled_target_operator_grammar_dossier_only",
        "provenance_notes": (
            "External annotation/literature/state evidence only; not a PDB contact table, "
            "not an internal runtime report, not native-metric selection, no MD."
        ),
        "evidence_statement": evidence_statement,
        "grammar_role": grammar_role,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "claim_allowed": False,
    }


def _kcsa_rows() -> list[dict[str, Any]]:
    target = "KcsA"
    return [
        _row(
            target,
            "KCSA_UNIPROT_P0A334_IDENTITY",
            "sequence_or_family_identity",
            "UniProt sequence/features/function annotations",
            "UniProtKB P0A334 pH-gated potassium channel KcsA",
            "https://www.uniprot.org/uniprotkb/P0A334/entry",
            "UniProtKB reviewed P0A334; last annotation update 2026-01-28; accessed 2026-06-11",
            "UniProt identifies KcsA as a pH-gated potassium channel from Streptomyces lividans and places it in the potassium channel family.",
            "KcsA target identity and potassium-channel family class.",
        ),
        _row(
            target,
            "KCSA_UNIPROT_P0A334_K_SELECTIVITY",
            "ion_selectivity_context",
            "UniProt sequence/features/function annotations",
            "UniProtKB P0A334 K+ selectivity/function annotation",
            "https://rest.uniprot.org/uniprotkb/P0A334.json",
            "UniProtKB reviewed P0A334; last annotation update 2026-01-28; accessed 2026-06-11",
            "UniProt function annotation reports pH-gated potassium-channel activity and monovalent cation preference K(+) > Rb(+) > NH4(+) >> Na(+) > Li(+).",
            "K+ selectivity grammar; not a coordinate-derived ion-contact table.",
        ),
        _row(
            target,
            "KCSA_UNIPROT_P0A334_MEMBRANE",
            "membrane_topology_context",
            "UniProt sequence/features/function annotations",
            "UniProtKB P0A334 membrane topology annotation",
            "https://rest.uniprot.org/uniprotkb/P0A334.json",
            "UniProtKB reviewed P0A334; last annotation update 2026-01-28; accessed 2026-06-11",
            "UniProt annotates KcsA as a cell-membrane multi-pass membrane protein with transmembrane helices and notes residues 62-79 are situated in the membrane and important for channel structure/properties.",
            "Membrane/transmembrane operator context before any native-coordinate selection.",
        ),
        _row(
            target,
            "KCSA_SEQUENCE_TVGYG_FILTER_SIGNATURE",
            "filter_or_signature_context",
            "external sequence conservation signatures",
            "KcsA canonical potassium-channel TVGYG selectivity-filter signature",
            "UniProtKB P0A334 canonical sequence; Heginbotham L et al. Biophys J. 1994;66:1061-1067; PubMed:8038378",
            "accessed_2026-06-11",
            "The KcsA sequence contains the canonical potassium-channel TVGYG filter signature, used here only as sequence/family grammar and not as a coordinate-contact constraint.",
            "Filter/signature grammar for K+-channel pore context.",
        ),
    ]


def _xcl1_rows() -> list[dict[str, Any]]:
    target = "XCL1_lymphotactin"
    return [
        _row(
            target,
            "XCL1_UNIPROT_P47992_IDENTITY",
            "metamorphic_two_state_context",
            "UniProt sequence/features/function annotations",
            "UniProtKB P47992 lymphotactin / XCL1",
            "https://www.uniprot.org/uniprotkb/P47992/entry",
            "UniProtKB reviewed P47992; accessed 2026-06-11",
            "UniProt identifies human XCL1/lymphotactin as a C motif chemokine ligand, establishing the sequence target and biological context.",
            "Target identity for the metamorphic state-switch dossier.",
        ),
        _row(
            target,
            "XCL1_LTN10_XCR1_FUNCTION",
            "state_A_function_context",
            "literature-derived state/function annotations",
            "XCL1/Ltn10 chemokine-like monomer XCR1 receptor-binding state",
            "Tuinstra RL et al. Interconversion between two unrelated protein folds in the lymphotactin native state. Proc Natl Acad Sci U S A. 2008;105:5057-5062.",
            "literature_annotation_accessed_2026-06-11",
            "The chemokine-like monomeric XCL1/Ltn10 state is the XCR1 receptor-binding chemokine state.",
            "State A grammar: chemokine-like monomer receptor-binding state.",
        ),
        _row(
            target,
            "XCL1_LTN40_GAG_HIV_FUNCTION",
            "state_B_function_context",
            "literature-derived state/function annotations",
            "XCL1/Ltn40 beta-sandwich dimer GAG/HIV-inhibitory state",
            "Tuinstra RL et al. Proc Natl Acad Sci U S A. 2008;105:5057-5062; Dishman AF and Volkman BF. Curr Opin Struct Biol. 2018;50:90-98.",
            "literature_annotation_accessed_2026-06-11",
            "The beta-sandwich dimeric XCL1/Ltn40 state is described as a distinct biologically relevant state with glycosaminoglycan-binding/HIV-inhibitory context.",
            "State B grammar: beta-sandwich dimer state with distinct function.",
        ),
        _row(
            target,
            "XCL1_NO_MIXED_STATE_POOLING",
            "no_mixed_state_pooling_rule",
            "literature-derived state/function annotations",
            "XCL1 metamorphic two-state separation rule",
            "Tuinstra RL et al. Proc Natl Acad Sci U S A. 2008;105:5057-5062; Dishman AF and Volkman BF. Curr Opin Struct Biol. 2018;50:90-98.",
            "literature_annotation_accessed_2026-06-11",
            "Because XCL1 interconverts between two unrelated native-state folds with different biological roles, V36 forbids pooling state-A and state-B evidence into one single-fold grammar.",
            "State-switch grammar firewall: do not mix incompatible state evidence.",
        ),
    ]


def _snca_rows() -> list[dict[str, Any]]:
    target = "alpha_synuclein_SNCA"
    return [
        _row(
            target,
            "SNCA_UNIPROT_P37840_IDENTITY",
            "ensemble_context",
            "UniProt sequence/features/function annotations",
            "UniProtKB P37840 alpha-synuclein / SNCA",
            "https://www.uniprot.org/uniprotkb/P37840/entry",
            "UniProtKB reviewed P37840; last annotation update 2026-06-10; accessed 2026-06-11",
            "UniProt identifies human SNCA as alpha-synuclein, a neuronal protein with membrane-associated multimeric functional context.",
            "Target identity and biological ensemble context.",
        ),
        _row(
            target,
            "SNCA_DISPROT_DISORDER_FREE_STATE",
            "intrinsic_disorder_context",
            "DisProt disorder annotations",
            "DisProt alpha-synuclein disorder annotation",
            "https://disprot.org/ ; alpha-synuclein / UniProtKB P37840 disorder annotation",
            "DisProt accessed 2026-06-11",
            "DisProt records alpha-synuclein as intrinsically disordered/unstructured in the free state.",
            "Free-state IDP grammar: disorder is the default, not a missing coordinate model.",
        ),
        _row(
            target,
            "SNCA_MEMBRANE_BOUND_HELIX_CONTEXT",
            "disorder_to_order_context",
            "literature-derived state/function annotations",
            "Alpha-synuclein membrane-bound alpha-helical transition",
            "Davidson WS et al. Stabilization of alpha-synuclein secondary structure upon binding to synthetic membranes. J Biol Chem. 1998;273:9443-9449.",
            "literature_annotation_accessed_2026-06-11",
            "Alpha-synuclein can undergo disorder-to-order transition toward alpha-helical structure when bound to lipid membranes.",
            "Disorder-to-order grammar without selecting a single native fold.",
        ),
        _row(
            target,
            "SNCA_NO_SINGLE_NATIVE_FOLD_RULE",
            "no_single_native_fold_rule",
            "experimentally described motifs/state labels",
            "Alpha-synuclein ensemble-not-single-fold rule",
            "DisProt alpha-synuclein / UniProtKB P37840; Davidson WS et al. J Biol Chem. 1998;273:9443-9449.",
            "accessed_2026-06-11",
            "SNCA is treated as an intrinsically disordered ensemble/disorder-to-order problem; V36 forbids forcing it into a single compact native-fold grammar.",
            "IDP firewall: ensemble grammar, no single native fold claim.",
        ),
    ]


TARGET_SPECS: dict[str, dict[str, Any]] = {
    "KcsA": {
        "display_name": "KcsA pH-gated potassium channel",
        "problem_class": "membrane_pore_filter_oligomeric_interface_problem",
        "grammar_focus": "membrane channel pore/filter/interface grammar",
        "required_buckets": [
            "sequence_or_family_identity",
            "ion_selectivity_context",
            "membrane_topology_context",
            "filter_or_signature_context",
        ],
        "grammar_rules": {
            "coordinate_contact_tables_allowed": False,
            "native_metrics_before_selection_allowed": False,
            "annotation_promoted_to_constraints": False,
        },
        "rows": _kcsa_rows,
        "rejected_sources": [
            {
                "source_name": "KcsA 1BL8 pore/filter contact CSV",
                "file_path": "data/external_constraints/KcsA/pore_filter/kcsa_1bl8_pore_filter_external_contacts.csv",
                "rejection_reason": "coordinate-derived V33/V34 contact table; not allowed in V36 non-coordinate dossier",
            },
            {
                "source_name": "KcsA 1BL8 assembly/interface contact CSV",
                "file_path": "data/external_constraints/KcsA/assembly_interface/kcsa_1bl8_assembly_interface_external_contacts.csv",
                "rejection_reason": "coordinate-derived V33/V34 contact table; not allowed in V36 non-coordinate dossier",
            },
        ],
    },
    "XCL1_lymphotactin": {
        "display_name": "XCL1 / lymphotactin metamorphic chemokine",
        "problem_class": "metamorphic_two_fold_state_switch_problem",
        "grammar_focus": "state-switch grammar with state-specific evidence separation",
        "required_buckets": [
            "state_A_function_context",
            "state_B_function_context",
            "metamorphic_two_state_context",
            "no_mixed_state_pooling_rule",
        ],
        "grammar_rules": {
            "mixed_state_pooling_allowed": False,
            "state_A_label": "chemokine_like_monomer_XCR1_binding_state",
            "state_B_label": "beta_sandwich_dimer_GAG_HIV_inhibitory_state",
            "native_metrics_before_selection_allowed": False,
        },
        "rows": _xcl1_rows,
        "rejected_sources": [
            {
                "source_name": "Mixed XCL1 state pool",
                "file_path": None,
                "rejection_reason": "state-A and state-B evidence must remain separated; mixed-state pooling is a blocked grammar error",
            },
        ],
    },
    "alpha_synuclein_SNCA": {
        "display_name": "Alpha-synuclein / SNCA",
        "problem_class": "intrinsically_disordered_ensemble_problem",
        "grammar_focus": "IDP ensemble/disorder-to-order grammar",
        "required_buckets": [
            "intrinsic_disorder_context",
            "ensemble_context",
            "disorder_to_order_context",
            "no_single_native_fold_rule",
        ],
        "grammar_rules": {
            "single_native_fold_model_allowed": False,
            "ensemble_not_single_fold": True,
            "aggregation_context_is_not_folding_solved_claim": True,
            "native_metrics_before_selection_allowed": False,
        },
        "rows": _snca_rows,
        "rejected_sources": [
            {
                "source_name": "Single compact native-fold SNCA model",
                "file_path": None,
                "rejection_reason": "SNCA is an IDP/ensemble problem; forcing one native-fold grammar is blocked",
            },
        ],
    },
}


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "target",
        "source_id",
        "evidence_bucket",
        "source_type",
        "source_name",
        "source_url_or_citation",
        "source_date_or_version",
        "evidence_statement",
        "grammar_role",
        "coordinate_derived",
        "internal_runtime_source",
        "native_metrics_used_for_selection",
        "coordinate_truth_used_before_selection",
        "claim_allowed",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})


def _bucket_status(rows: list[dict[str, Any]], required: list[str]) -> dict[str, Any]:
    present = sorted({str(row.get("evidence_bucket")) for row in rows})
    missing = [bucket for bucket in required if bucket not in present]
    return {
        "required_buckets": required,
        "present_buckets": present,
        "missing_buckets": missing,
        "complete": not missing,
    }


def build_target_dossier(target: str) -> dict[str, Any]:
    spec = TARGET_SPECS[target]
    rows = spec["rows"]()
    bucket_status = _bucket_status(rows, spec["required_buckets"])
    source_manifest = {
        "kind": "V36_TARGET_SOURCE_MANIFEST_v0",
        "target": target,
        "display_name": spec["display_name"],
        "claim_allowed": False,
        "new_MD_allowed": False,
        "source_row_count": len(rows),
        "rows": rows,
    }
    evidence_dossier = {
        "kind": "V36_TARGET_EVIDENCE_DOSSIER_v0",
        "target": target,
        "display_name": spec["display_name"],
        "problem_class": spec["problem_class"],
        "grammar_focus": spec["grammar_focus"],
        "required_buckets": spec["required_buckets"],
        "bucket_status": bucket_status,
        "grammar_rules": spec["grammar_rules"],
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "locked_interpretation": (
            "This dossier is target-specific operator grammar evidence only. It does not convert annotations into constraints, "
            "does not use native coordinates for selection, does not run MD, and does not make a folding-solved claim."
        ),
    }
    acquisition_log = {
        "kind": "V36_ACQUISITION_LOG_v0",
        "target": target,
        "generated_by": "scripts/build_v36_real_evidence_dossiers_v0.py",
        "acquisition_mode": "curated_external_source_dossier_small_auditable_no_large_downloads",
        "accessed_date": "2026-06-11",
        "source_count": len(rows),
        "source_ids": [row["source_id"] for row in rows],
        "notes": [
            "Sources are external annotations/literature/state labels recorded as provenance rows.",
            "No PDB coordinate contact tables, AlphaFold/ESMFold/RoseTTAFold coordinates, native metrics, internal runtime reports, or MD outputs are used as evidence.",
        ],
    }
    rejected_sources = {
        "kind": "V36_REJECTED_SOURCES_v0",
        "target": target,
        "rows": spec["rejected_sources"],
    }
    return {
        "source_manifest": source_manifest,
        "evidence_dossier": evidence_dossier,
        "evidence_rows": rows,
        "acquisition_log": acquisition_log,
        "rejected_sources": rejected_sources,
    }


def write_target_dossier(root: Path, target: str) -> dict[str, Path]:
    package = build_target_dossier(target)
    out_dir = root / target
    paths = {
        "source_manifest": out_dir / "source_manifest.json",
        "evidence_dossier": out_dir / "evidence_dossier.json",
        "evidence_table": out_dir / "evidence_table.csv",
        "acquisition_log": out_dir / "acquisition_log.json",
        "rejected_sources": out_dir / "rejected_sources.json",
    }
    _write_json(paths["source_manifest"], package["source_manifest"])
    _write_json(paths["evidence_dossier"], package["evidence_dossier"])
    _write_csv(paths["evidence_table"], package["evidence_rows"])
    _write_json(paths["acquisition_log"], package["acquisition_log"])
    _write_json(paths["rejected_sources"], package["rejected_sources"])
    return paths


def build_all_dossiers(root: Path = DATA_ROOT) -> dict[str, Any]:
    artifacts: dict[str, dict[str, str]] = {}
    source_counts: dict[str, int] = {}
    for target in TARGET_ORDER:
        paths = write_target_dossier(root, target)
        artifacts[target] = {key: str(path) for key, path in paths.items()}
        source_counts[target] = len(TARGET_SPECS[target]["rows"]())
    return {
        "kind": "V36_REAL_EVIDENCE_DOSSIER_ACQUISITION_v0",
        "target_count": len(TARGET_ORDER),
        "targets": TARGET_ORDER,
        "source_counts_by_target": source_counts,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
        "artifacts": artifacts,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build V36 real evidence dossiers.")
    parser.add_argument("--out-root", type=Path, default=DATA_ROOT)
    args = parser.parse_args()
    summary = build_all_dossiers(args.out_root)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
