from __future__ import annotations

"""Coarse Protein Esperanto operator grammar and simulator.

This module is intentionally not an atomistic simulator.  It turns allowed
sequence/evidence into a coarse state field, selects a mechanism grammar,
builds operator activations, and evolves mechanism-level observables that can
be sealed before holdout validation.
"""

from copy import deepcopy
from hashlib import sha256
import json
from math import sin
from statistics import mean
from typing import Any, Iterable


PURE_NON_COORDINATE = "pure_non_coordinate"
SPATIAL_PROXY_NON_COORDINATE = "spatial_proxy_non_coordinate"
COORDINATE_DERIVED = "coordinate_derived"
INTERNAL_RUNTIME = "internal_runtime"

EVIDENCE_CLASSES = [
    PURE_NON_COORDINATE,
    SPATIAL_PROXY_NON_COORDINATE,
    COORDINATE_DERIVED,
    INTERNAL_RUNTIME,
]

MECHANISM_CLASSES = [
    "globular_closure",
    "intrinsic_disorder_phase_separation",
    "disorder_boundary_and_fold_upon_binding",
    "beta_closure_topology",
    "multidomain_allosteric_architecture",
    "membrane_multidomain_folding_proteostasis",
    "metamorphic_fold_switching",
    "short_region_host_interface_hijacking",
    "fold_upon_binding_disorder",
    "cofactor_ligand_assisted_stabilization",
    "metal_cluster_and_ligand_locked_basin",
    "assembly_required_folding",
    "oligomerization_controlled_folding",
    "insufficient_evidence_clean_abstain",
]

PROCESS_CLASSES = [
    "two_state",
    "intermediate_bearing",
    "multi_basin",
    "disorder_biased",
    "disorder_boundary",
    "beta_closure",
    "multidomain_allostery",
    "fold_upon_binding",
]

UNIVERSAL_MARKS = [
    "charge",
    "hydrophobicity",
    "aromatic_density",
    "low_complexity",
    "disorder_tendency",
    "secondary_structure_tendency",
    "motif_presence",
    "proline_glycine_cysteine_effects",
    "membrane_tendency",
    "domain_boundary_tendency",
    "interface_motif_tendency",
    "evolutionary_conservation_if_allowed",
]

UNIVERSAL_OPERATORS = [
    "closure_operator",
    "repulsion_operator",
    "frustration_operator",
    "disorder_operator",
    "phase_operator",
    "interface_operator",
    "membrane_pressure_operator",
    "dual_basin_switch_operator",
    "proteostasis_operator",
    "host_hijack_operator",
]

STATE_VARIABLES = [
    "residue_exposure",
    "segment_compaction",
    "contact_probability",
    "operator_activation",
    "frustration",
    "state_basin_occupancy",
    "interface_readiness",
    "disorder_order_balance",
    "proteostasis_routing",
    "assembly_required_core",
    "partner_completed_core",
    "interface_buried_hydrophobicity",
    "monomer_incomplete_topology",
    "coiled_coil_register_dependency",
    "domain_swap_candidate",
    "biological_oligomer_context",
    "generic_complex_only",
    "assembly_ambiguous",
    "metal_cluster_geometry",
    "coordination_shell_integrity",
    "ligand_locked_basin",
    "apo_holo_basin_shift",
    "generic_ligand_only",
    "metal_ligand_ambiguous",
    "IDR_boundary",
    "structured_domain_plus_IDR_tail",
    "fold_upon_binding_region",
    "phase_prone_low_complexity",
    "flexible_loop_not_disorder",
    "disorder_with_local_motif",
    "closed_beta_topology",
    "strand_register",
    "beta_sheet_closure",
    "soluble_beta_barrel",
    "membrane_beta_barrel",
    "beta_propeller_repeat_closure",
    "beta_sandwich_core",
    "jelly_roll_wrap",
    "greek_key_beta_lock",
    "beta_helix_solenoid_stack",
    "alpha_beta_barrel_distinction",
    "open_beta_sheet_ambiguous",
    "closed_beta_confident",
    "closed_beta_ambiguous",
    "strand_register_insufficient",
    "beta_topology_conflict",
    "multidomain_allostery",
    "domain_boundary",
    "hinge_region",
    "interdomain_lock",
    "allosteric_basin_shift",
    "domain_reorientation",
    "modular_architecture",
    "domain_swapping",
]

ENVIRONMENTAL_PRESSURES = [
    "concentration",
    "salt",
    "rna_partner",
    "membrane",
    "proteostasis_quality_control",
    "binding_partner",
    "ligand_or_cofactor",
    "release_or_activation_context",
]

COFACTOR_CONTEXT_TOKENS = [
    "cofactor_context",
    "ligand_context",
    "metal_context",
    "heme_context",
    "nucleotide_context",
]

METAL_CLUSTER_CONTEXT_TOKENS = [
    "metal_cluster_geometry",
    "metal cluster",
    "metal-cluster",
    "iron-sulfur",
    "iron sulfur",
    "fe-s",
    "metalloprotein",
    "metal-binding",
    "metal binding",
    "coordination_shell_integrity",
]

LIGAND_LOCKED_CONTEXT_TOKENS = [
    "ligand_locked_basin",
    "ligand-locked basin",
    "cofactor_locked_basin",
    "cofactor-locked basin",
    "nucleotide_context",
    "flavin",
    "substrate-bound",
    "substrate bound",
    "apo_holo_basin_shift",
]

DISORDER_BOUNDARY_CONTEXT_TOKENS = [
    "disorder_context",
    "IDR_boundary",
    "idr boundary",
    "structured_domain_plus_IDR_tail",
    "structured domain plus idr tail",
    "fold_upon_binding_region",
    "fold-upon-binding region",
    "phase_prone_low_complexity",
    "phase-prone low complexity",
    "flexible_loop_not_disorder",
    "flexible loop not disorder",
    "disorder_with_local_motif",
    "disorder with local motif",
    "disordered region tendency",
    "low complexity tendency",
]

BETA_CLOSURE_CONTEXT_TOKENS = [
    "closed_beta_topology",
    "closed beta topology",
    "strand_register",
    "strand register",
    "beta_sheet_closure",
    "beta sheet closure",
    "soluble_beta_barrel",
    "soluble beta barrel",
    "membrane_beta_barrel",
    "membrane beta barrel",
    "outer membrane beta barrel",
    "beta_propeller_repeat_closure",
    "beta propeller",
    "beta-propeller",
    "beta_sandwich_core",
    "beta sandwich",
    "jelly_roll_wrap",
    "jelly roll",
    "jelly-roll",
    "greek_key_beta_lock",
    "greek key",
    "greek-key",
    "beta_helix_solenoid_stack",
    "beta helix",
    "beta-helix",
    "beta solenoid",
    "beta-solenoid",
    "alpha_beta_barrel_distinction",
    "alpha beta barrel",
    "alpha-beta barrel",
]

BETA_AMBIGUOUS_CONTEXT_TOKENS = [
    "open_beta_sheet_ambiguous",
    "open beta sheet ambiguous",
    "strand_register_insufficient",
    "strand register insufficient",
    "beta_topology_conflict",
    "beta topology conflict",
    "beta propensity only",
    "weak beta closure evidence",
]

BETA_TOPOLOGY_STATE_VARIABLES = [
    "closed_beta_topology",
    "soluble_beta_barrel",
    "membrane_beta_barrel",
    "beta_propeller_repeat_closure",
    "beta_sandwich_core",
    "jelly_roll_wrap",
    "greek_key_beta_lock",
    "beta_helix_solenoid_stack",
    "alpha_beta_barrel_distinction",
]

MULTIDOMAIN_ALLOSTERIC_CONTEXT_TOKENS = [
    "multidomain_allostery",
    "multidomain allostery",
    "multidomain architecture",
    "multi-domain architecture",
    "domain_boundary",
    "domain boundary",
    "hinge_region",
    "hinge region",
    "interdomain_lock",
    "interdomain lock",
    "allosteric_basin_shift",
    "allosteric basin shift",
    "domain_reorientation",
    "domain reorientation",
    "modular_architecture",
    "modular architecture",
    "domain_swapping",
    "domain swapping",
    "domain-swapped",
    "domain swapped",
]

MULTIDOMAIN_ALLOSTERIC_STATE_VARIABLES = [
    "multidomain_allostery",
    "domain_boundary",
    "hinge_region",
    "interdomain_lock",
    "allosteric_basin_shift",
    "domain_reorientation",
    "modular_architecture",
    "domain_swapping",
]

OLIGOMER_CONTEXT_TOKENS = [
    "oligomer_context",
    "assembly_context",
    "partner_copy_context",
    "heteromeric_context",
    "homomeric_context",
]

ASSEMBLY_REQUIRED_CONTEXT_TOKENS = [
    "assembly_required_core",
    "assembly_required_folding",
    "partner_completed_core",
    "partner-completed fold",
    "partner-stabilized fold",
    "partner stabilized fold",
    "interface_buried_hydrophobicity",
    "monomer_incomplete_topology",
    "coiled_coil_register_dependency",
    "coiled-coil register",
    "coiled coil register",
    "leucine zipper",
    "domain_swap_candidate",
    "domain-swapped fold",
    "domain swapped fold",
    "obligate assembly",
    "obligate oligomer",
    "biological_oligomer_context",
    "assembly-dependent folding",
    "assembly dependent folding",
    "multimer-stabilized basin",
    "multimer stabilized basin",
]

NEGATIVE_ASSEMBLY_CONTEXT_TOKENS = [
    "not assembly_required_core",
    "not assembly_required_folding",
    "not assembly required",
    "not assembly-required",
    "not obligate assembly",
    "not partner_completed_core",
    "not partner-completed",
    "not biological_oligomer_context",
]

GENERIC_COMPLEX_ONLY_TOKENS = [
    "generic_complex_only",
    "generic complex only",
    "generic complex_not_assembly",
    "generic complex not assembly",
]

ASSEMBLY_AMBIGUITY_TOKENS = [
    "assembly_ambiguous",
    "assembly_topology_ambiguous",
    "weak assembly evidence",
    "ambiguous assembly evidence",
]

STRONG_MEMBRANE_CONTEXT_TOKENS = [
    "membrane_context_strong",
    "transmembrane_context",
    "channel_context",
    "transporter_context",
    "receptor_membrane_context",
]

NEGATIVE_MEMBRANE_TOPOLOGY_TOKENS = [
    "no membrane topology",
    "no transmembrane",
    "not transmembrane",
    "without transmembrane",
    "no bilayer-spanning",
    "no bilayer spanning",
    "hydrophobicity-alone",
    "soluble_hydrophobic_core_context",
    "cofactor-buried hydrophobic pocket",
    "oligomeric interface hydrophobicity",
    "no opm/pdbtm/memprotmd transmembrane assignment",
]

PERIPHERAL_MEMBRANE_CONTEXT_TOKENS = [
    "peripheral membrane",
    "membrane-associated",
    "membrane associated",
    "monotopic/peripheral",
    "monotopic",
    "amphipathic peripheral helix",
    "lipid anchor",
    "lipid-facing surface",
]

HYDROPHOBIC = frozenset("AILMFWVYC")
AROMATIC = frozenset("FWY")
POSITIVE = frozenset("KRH")
NEGATIVE = frozenset("DE")
POLAR = frozenset("STNQ")
BREAKERS = frozenset("PG")
HELIX = frozenset("AEHKLMQR")
BETA = frozenset("VIFYWT")
LOW_COMPLEXITY_BIASED = frozenset("GSQNYP")


GRAMMAR_RULES: dict[str, dict[str, Any]] = {
    "globular_closure": {
        "marks": ["hydrophobic", "aromatic", "secondary_structure_tendency"],
        "pressures": ["solvent", "folding_nucleus"],
        "operators": ["closure_operator", "frustration_operator"],
        "state_change": "expanded chain to compact nucleus/contact-enriched ensemble",
        "testable_effect": "core-region perturbation weakens compaction and long-range topology",
        "null_control": "generic hydrophobicity alone cannot validate target-specific closure",
        "falsification_rule": "independent holdout shows no compaction or wrong operator regions",
    },
    "intrinsic_disorder_phase_separation": {
        "marks": ["low_complexity", "aromatic", "charged", "glycine/proline-rich"],
        "pressures": ["concentration", "salt", "rna_partner"],
        "operators": ["disorder_operator", "phase_operator", "repulsion_operator"],
        "state_change": "expanded disorder to dynamic phase-prone ensemble without single stable fold",
        "testable_effect": "aromatic/charge or condition perturbations shift condensation tendency",
        "null_control": "forcing compact single-fold grammar fails",
        "falsification_rule": "holdout shows condition-independent stable globular dominance",
    },
    "disorder_boundary_and_fold_upon_binding": {
        "marks": [
            "IDR_boundary",
            "structured_domain_plus_IDR_tail",
            "fold_upon_binding_region",
            "phase_prone_low_complexity",
            "disorder_with_local_motif",
        ],
        "pressures": ["binding_partner", "concentration", "salt"],
        "operators": ["disorder_operator", "interface_operator", "phase_operator"],
        "state_change": "generic oligomer or compact-fold readout to explicit disorder-boundary ensemble with partner-conditioned local order",
        "testable_effect": "boundary, local motif, or partner perturbation changes ordering without promoting a whole-chain solved fold",
        "null_control": "generic oligomer copy metadata cannot override explicit IDR/low-complexity boundary evidence",
        "falsification_rule": "holdout shows a complete stable globular or assembly-required fold without persistent disorder boundary",
    },
    "beta_closure_topology": {
        "marks": [
            "closed_beta_topology",
            "strand_register",
            "beta_sheet_closure",
            "soluble_beta_barrel",
            "membrane_beta_barrel",
            "beta_propeller_repeat_closure",
            "beta_sandwich_core",
            "jelly_roll_wrap",
            "greek_key_beta_lock",
            "beta_helix_solenoid_stack",
            "alpha_beta_barrel_distinction",
        ],
        "pressures": ["solvent", "membrane", "strand_register", "repeat_closure"],
        "operators": ["closure_operator", "interface_operator", "frustration_operator"],
        "state_change": "beta propensity to subtype-specific strand-register closure without confusing membrane, propeller, sandwich, solenoid, or alpha/beta barrels",
        "testable_effect": "strand-register or blade/repeat perturbations weaken the selected beta closure while preserving non-beta priority classes",
        "null_control": "beta-rich sequence alone is not topology and must abstain when closure evidence is weak",
        "falsification_rule": "holdout shows open beta sheet, wrong beta subtype, or protected membrane/assembly/metal/disorder context",
    },
    "multidomain_allosteric_architecture": {
        "marks": [
            "multidomain_allostery",
            "domain_boundary",
            "hinge_region",
            "interdomain_lock",
            "allosteric_basin_shift",
            "domain_reorientation",
            "modular_architecture",
            "domain_swapping",
        ],
        "pressures": ["binding_partner", "interdomain_coupling", "release_or_activation_context"],
        "operators": ["interface_operator", "dual_basin_switch_operator", "closure_operator", "frustration_operator"],
        "state_change": "single-domain closure shortcut to domain-boundary, hinge, lock, and allosteric-basin architecture",
        "testable_effect": "hinge, interface, or interdomain-lock perturbation shifts domain orientation and allosteric basin occupancy",
        "null_control": "generic domain text or copy-count metadata cannot validate an allosteric architecture without explicit domain-boundary/hinge/swap context",
        "falsification_rule": "holdout shows independent single-domain closure, solved prior membrane/assembly/metal/disorder/beta context, or no interdomain operator dependence",
    },
    "membrane_multidomain_folding_proteostasis": {
        "marks": ["membrane_segment", "domain_boundary", "mutation_site", "interface"],
        "pressures": ["membrane", "proteostasis_quality_control", "interdomain_coupling"],
        "operators": ["membrane_pressure_operator", "closure_operator", "proteostasis_operator", "interface_operator"],
        "state_change": "domain-local destabilization to interface/proteostasis routing defect",
        "testable_effect": "NBD1-only correction is partial; interface/proteostasis correction strengthens rescue",
        "null_control": "single local mutation-only shortcut fails",
        "falsification_rule": "holdout shows purely local defect with no interface/proteostasis dependence",
    },
    "metamorphic_fold_switching": {
        "marks": ["state_basin", "partner_context", "secondary_structure_conflict"],
        "pressures": ["release_or_activation_context", "binding_partner"],
        "operators": ["dual_basin_switch_operator", "interface_operator", "frustration_operator"],
        "state_change": "one context-favored basin to competing alpha/beta state distribution",
        "testable_effect": "context or partner perturbations shift basin occupancy",
        "null_control": "single averaged consensus fold is rejected",
        "falsification_rule": "holdout supports only one state under all contexts",
    },
    "short_region_host_interface_hijacking": {
        "marks": ["short_linear_motif", "c_terminal_interface", "host_binding_surface"],
        "pressures": ["host_binding_partner", "localization_context"],
        "operators": ["host_hijack_operator", "interface_operator", "disorder_operator"],
        "state_change": "short exposed motif to host-interface capture without global fold requirement",
        "testable_effect": "C-terminal/interface perturbations weaken host hijacking consequences",
        "null_control": "generic viral accessory annotation cannot validate the packet",
        "falsification_rule": "wrong-region perturbations dominate while predicted interface is inert",
    },
    "fold_upon_binding_disorder": {
        "marks": ["disorder_tendency", "partner_motif", "amphipathic_segment"],
        "pressures": ["binding_partner"],
        "operators": ["disorder_operator", "interface_operator", "closure_operator"],
        "state_change": "free disorder to partner-conditioned local ordering",
        "testable_effect": "partner removal restores disorder and weakens interface readiness",
        "null_control": "bound-state order cannot be promoted to free-state solved fold",
        "falsification_rule": "free-state holdout shows stable standalone fold",
    },
    "cofactor_ligand_assisted_stabilization": {
        "marks": ["cofactor_motif", "polar_cluster", "burial_defect"],
        "pressures": ["ligand_or_cofactor"],
        "operators": ["closure_operator", "interface_operator"],
        "state_change": "weak apo basin to ligand-stabilized basin",
        "testable_effect": "cofactor removal destabilizes the predicted basin",
        "null_control": "generic ligand annotation alone cannot define full grammar",
        "falsification_rule": "ligand-independent stability dominates all conditions",
    },
    "metal_cluster_and_ligand_locked_basin": {
        "marks": ["metal_cluster_geometry", "coordination_shell_integrity", "ligand_locked_basin"],
        "pressures": ["ligand_or_cofactor", "redox_or_coordination_context"],
        "operators": ["interface_operator", "closure_operator", "frustration_operator"],
        "state_change": "generic cofactor basin to metal/ligand-locked basin with apo-holo separation",
        "testable_effect": "metal/cofactor removal or coordination-shell disruption unlocks the basin and increases frustration",
        "null_control": "generic cofactor annotation without metal, heme, nucleotide, or ligand-lock evidence remains generic cofactor",
        "falsification_rule": "holdout shows ligand-independent folding or no geometry-locked basin effect",
    },
    "assembly_required_folding": {
        "marks": [
            "assembly_required_core",
            "partner_completed_core",
            "interface_buried_hydrophobicity",
            "monomer_incomplete_topology",
        ],
        "pressures": ["binding_partner", "partner_copy", "concentration"],
        "operators": ["interface_operator", "closure_operator", "frustration_operator"],
        "state_change": "monomer-incomplete topology to partner-completed assembly basin",
        "testable_effect": "partner/interface disruption exposes incomplete monomer topology and weakens the assembled basin",
        "null_control": "generic complex annotation alone cannot validate assembly-required folding",
        "falsification_rule": "holdout shows a complete monomeric fold, true membrane topology, or ligand-stabilized core without assembly dependence",
    },
    "oligomerization_controlled_folding": {
        "marks": ["repeat_interface", "oligomer_motif", "burial_surface"],
        "pressures": ["concentration", "partner_copy"],
        "operators": ["interface_operator", "closure_operator"],
        "state_change": "monomeric partial order to oligomer-stabilized interface basin",
        "testable_effect": "interface or concentration perturbation shifts folding/assembly",
        "null_control": "monomer-only grammar fails",
        "falsification_rule": "holdout shows no assembly dependence",
    },
}


def stable_hash(data: Any) -> str:
    payload = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return sha256(payload).hexdigest()


def bounded(value: float) -> float:
    return round(max(0.0, min(1.0, float(value))), 6)


def _avg(values: Iterable[float]) -> float:
    values = list(values)
    return round(mean(values), 6) if values else 0.0


def normalize_sequence(sequence: str) -> str:
    normalized = "".join(aa for aa in sequence.upper() if "A" <= aa <= "Z")
    if not normalized:
        raise ValueError("sequence must contain at least one residue")
    return normalized


def _window(sequence: str, index0: int, radius: int = 3) -> str:
    return sequence[max(0, index0 - radius) : min(len(sequence), index0 + radius + 1)]


def _charge(aa: str) -> float:
    if aa in POSITIVE:
        return 1.0
    if aa in NEGATIVE:
        return -1.0
    return 0.0


def _secondary_propensity(aa: str, window: str) -> str:
    helix_score = sum(1.0 if residue in HELIX else 0.25 for residue in window) / len(window)
    beta_score = sum(1.0 if residue in BETA else 0.25 for residue in window) / len(window)
    breaker_score = sum(1.0 if residue in BREAKERS else 0.0 for residue in window) / len(window)
    if breaker_score >= 0.25:
        return "coil_or_disorder"
    if helix_score >= beta_score + 0.10:
        return "helix_prone"
    if beta_score >= helix_score + 0.05:
        return "beta_prone"
    if aa in POLAR or aa in LOW_COMPLEXITY_BIASED:
        return "coil_or_disorder"
    return "mixed"


def _residue_mark(sequence: str, index0: int) -> dict[str, Any]:
    aa = sequence[index0]
    window = _window(sequence, index0, radius=4)
    unique_fraction = len(set(window)) / len(window)
    biased_fraction = sum(1 for residue in window if residue in LOW_COMPLEXITY_BIASED) / len(window)
    hydrophobic_fraction = sum(1 for residue in window if residue in HYDROPHOBIC) / len(window)
    charged_fraction = sum(1 for residue in window if residue in POSITIVE or residue in NEGATIVE) / len(window)
    aromatic_fraction = sum(1 for residue in window if residue in AROMATIC) / len(window)
    pro_gly_fraction = sum(1 for residue in window if residue in BREAKERS) / len(window)
    low_complexity = unique_fraction <= 0.48 or biased_fraction >= 0.52
    disorder = bounded(0.45 * biased_fraction + 0.30 * pro_gly_fraction + 0.25 * charged_fraction - 0.15 * hydrophobic_fraction)
    membrane = bounded((hydrophobic_fraction - 0.45) / 0.35)
    interface = bounded(0.45 * hydrophobic_fraction + 0.35 * aromatic_fraction + 0.20 * charged_fraction)
    return {
        "position_index": index0 + 1,
        "residue_identity": aa,
        "charge_mark": _charge(aa),
        "hydrophobic_mark": aa in HYDROPHOBIC,
        "aromatic_mark": aa in AROMATIC,
        "polar_mark": aa in POLAR,
        "glycine_mark": aa == "G",
        "proline_mark": aa == "P",
        "cysteine_mark": aa == "C",
        "low_complexity_mark": low_complexity,
        "disorder_mark": disorder,
        "secondary_propensity_mark": _secondary_propensity(aa, window),
        "membrane_mark": membrane,
        "interface_mark": interface,
        "domain_boundary_tendency": bounded(pro_gly_fraction + charged_fraction * 0.3),
        "motif_presence": "aromatic_charged" if aa in AROMATIC or aa in POSITIVE or aa in NEGATIVE else "none",
        "evolutionary_conservation_if_allowed": "not_used_in_mvp",
    }


def build_sequence_field(sequence: str, *, segment_size: int = 12) -> dict[str, Any]:
    sequence = normalize_sequence(sequence)
    residues = [_residue_mark(sequence, index0) for index0 in range(len(sequence))]
    segments: list[dict[str, Any]] = []
    for start0 in range(0, len(sequence), segment_size):
        segment_residues = residues[start0 : start0 + segment_size]
        segment_seq = sequence[start0 : start0 + segment_size]
        segment_index = len(segments) + 1
        segments.append({
            "segment_id": f"S{segment_index:03d}",
            "start": start0 + 1,
            "end": start0 + len(segment_seq),
            "sequence": segment_seq,
            "charge_density": round(sum(row["charge_mark"] for row in segment_residues) / len(segment_residues), 6),
            "hydrophobic_density": _avg(1.0 if row["hydrophobic_mark"] else 0.0 for row in segment_residues),
            "aromatic_density": _avg(1.0 if row["aromatic_mark"] else 0.0 for row in segment_residues),
            "beta_propensity_density": _avg(1.0 if row["secondary_propensity_mark"] == "beta_prone" else 0.0 for row in segment_residues),
            "low_complexity_density": _avg(1.0 if row["low_complexity_mark"] else 0.0 for row in segment_residues),
            "disorder_density": _avg(row["disorder_mark"] for row in segment_residues),
            "membrane_density": _avg(row["membrane_mark"] for row in segment_residues),
            "interface_density": _avg(row["interface_mark"] for row in segment_residues),
            "pro_gly_density": _avg(1.0 if row["glycine_mark"] or row["proline_mark"] else 0.0 for row in segment_residues),
        })
    global_metrics = {
        "length": len(sequence),
        "net_charge_per_residue": round(sum(row["charge_mark"] for row in residues) / len(sequence), 6),
        "hydrophobic_density": _avg(1.0 if row["hydrophobic_mark"] else 0.0 for row in residues),
        "aromatic_density": _avg(1.0 if row["aromatic_mark"] else 0.0 for row in residues),
        "beta_propensity_density": _avg(1.0 if row["secondary_propensity_mark"] == "beta_prone" else 0.0 for row in residues),
        "low_complexity_density": _avg(1.0 if row["low_complexity_mark"] else 0.0 for row in residues),
        "mean_disorder": _avg(row["disorder_mark"] for row in residues),
        "mean_membrane": _avg(row["membrane_mark"] for row in residues),
        "mean_interface": _avg(row["interface_mark"] for row in residues),
        "coordinate_truth_used": False,
    }
    for residue in residues:
        residue["segment_id"] = f"S{((residue['position_index'] - 1) // segment_size) + 1:03d}"
        residue["current_state"] = "uninitialized"
        residue["state_confidence"] = 0.0
        residue["operator_activations"] = {}
    return {
        "kind": "PROTEIN_ESPERANTO_SEQUENCE_FIELD_v0",
        "sequence_length": len(sequence),
        "segment_size": segment_size,
        "required_marks": UNIVERSAL_MARKS,
        "residues": residues,
        "segments": segments,
        "global_metrics": global_metrics,
        "coordinate_truth_used": False,
    }


def _bounded_span(sequence_length: int, start: int, end: int) -> list[int]:
    return [max(1, min(sequence_length, int(start))), max(1, min(sequence_length, int(end)))]


def _route_span(sequence_length: int, route: str) -> list[int]:
    if route == "central":
        return _bounded_span(sequence_length, max(1, sequence_length // 4), max(1, (sequence_length * 3) // 4))
    if route == "beta_core":
        return _bounded_span(sequence_length, max(1, sequence_length // 3), max(1, sequence_length - max(8, sequence_length // 6)))
    if route == "mini_alpha_core":
        return _bounded_span(sequence_length, max(1, sequence_length // 10), max(1, (sequence_length * 7) // 10))
    if route == "alpha_core":
        return _bounded_span(sequence_length, max(1, sequence_length // 4), max(1, sequence_length - max(6, sequence_length // 8)))
    if route == "late_c":
        return _bounded_span(sequence_length, max(1, sequence_length - max(10, sequence_length // 4)), sequence_length)
    return _bounded_span(sequence_length, 1, min(sequence_length, 12))


def _process_region(label: str, family: str, span: list[int]) -> dict[str, Any]:
    return {"label": label, "family": family, "span": span}


def predict_process_route_profile(
    *,
    sequence: str,
    target_name: str = "",
    structural_class: str = "",
    architecture_hint: str = "",
    selected_mechanism_class: str = "",
) -> dict[str, Any]:
    """Predict a coarse process route without opening process holdouts.

    This is still a Protein Esperanto readout, not atomistic MD.  It chooses a
    route-level process profile from sequence plus allowed non-coordinate
    target metadata and intentionally returns abstention when the mechanism
    grammar itself abstains.
    """

    normalized = normalize_sequence(sequence)
    length = len(normalized)
    text = f"{target_name} {structural_class} {architecture_hint}".lower()
    if selected_mechanism_class == "insufficient_evidence_clean_abstain":
        process_class = "abstain_no_process_claim"
        decision = "abstain_recommended"
        early_family = "none"
        late_family = "none"
        early_span = _bounded_span(length, 1, min(length, 1))
        late_span = early_span
        nucleus = "no_nucleus_decision"
        sensitive_families: list[str] = []
    elif "ww" in text:
        process_class = "multi_basin"
        decision = "accepted_with_caution"
        early_family = "alternative_beta_hairpin"
        late_family = "sheet_locking"
        early_span = _route_span(length, "central")
        late_span = _route_span(length, "late_c")
        nucleus = "distributed_or_alternative_nucleus"
        sensitive_families = ["aromatic_beta_core", "alternative_beta_hairpin"]
    elif "larger" in text and ("alpha beta" in text or "alpha_beta" in text or "helical bundle" in text or "enzyme-like" in text):
        process_class = "intermediate_bearing"
        decision = "accepted_with_caution"
        early_family = "subdomain_intermediate_core"
        late_family = "delayed_helix_or_terminal_docking"
        early_span = _bounded_span(length, 1, max(1, (length * 9) // 10))
        late_span = _bounded_span(length, max(1, (length * 3) // 5), length)
        nucleus = "intermediate_with_nucleus"
        sensitive_families = ["intermediate_stabilizing_core", "hydrophobic_core", "helix_bundle_core"]
    elif structural_class == "beta" or "beta-sheet" in text or "beta sheet" in text or "sh3" in text:
        process_class = "two_state"
        decision = "accepted"
        early_family = "beta_hairpin_core"
        late_family = "terminal_strand_lock"
        early_span = _route_span(length, "beta_core")
        late_span = _route_span(length, "early_n")
        nucleus = "nucleus_present"
        sensitive_families = ["beta_hairpin_core", "turn_or_sheet_core"]
    elif structural_class == "alpha_beta":
        process_class = "two_state"
        decision = "accepted"
        early_family = "mixed_alpha_beta_core"
        late_family = "terminal_consolidation"
        early_span = _route_span(length, "central")
        late_span = _route_span(length, "late_c" if "ubiquitin" in text else "early_n")
        nucleus = "nucleus_present"
        sensitive_families = ["mixed_alpha_beta_core", "hydrophobic_core"]
    elif structural_class.startswith("alpha") or "helix" in text or "caged aromatic core" in text:
        process_class = "two_state"
        decision = "accepted"
        early_family = "helix_bundle_core"
        late_family = "terminal_consolidation"
        early_span = _route_span(length, "mini_alpha_core" if length <= 45 else "alpha_core")
        late_span = _route_span(length, "late_c" if length <= 28 else "early_n")
        nucleus = "nucleus_present"
        sensitive_families = ["helix_bundle_core", "hydrophobic_core", "aromatic_core"]
    else:
        process_class = "two_state"
        decision = "accepted"
        early_family = "mixed_alpha_beta_core"
        late_family = "terminal_consolidation"
        early_span = _route_span(length, "central")
        late_span = _route_span(length, "late_c" if "ubiquitin" in text else "early_n")
        nucleus = "nucleus_present"
        sensitive_families = ["mixed_alpha_beta_core", "hydrophobic_core"]
    return {
        "kind": "PROTEIN_ESPERANTO_PROCESS_ROUTE_PROFILE_v0",
        "process_decision": decision,
        "predicted_process_class": process_class,
        "predicted_early_forming_region": _process_region("predicted early-forming route region", early_family, early_span),
        "predicted_late_forming_region": _process_region("predicted late-forming consolidation region", late_family, late_span),
        "predicted_folding_nucleus_decision": nucleus,
        "predicted_mutation_sensitive_region_families": sensitive_families,
        "available_process_classes": PROCESS_CLASSES,
        "coordinate_truth_used_before_seal": False,
        "atomistic_md_executed": False,
    }


def _flatten_text(value: Any) -> str:
    if isinstance(value, dict):
        return " ".join(_flatten_text(item) for item in value.values())
    if isinstance(value, list):
        return " ".join(_flatten_text(item) for item in value)
    return str(value)


def classify_evidence_source(source: dict[str, Any]) -> str:
    if source.get("source_class") in EVIDENCE_CLASSES:
        return str(source["source_class"])
    text = _flatten_text(source).lower()
    if source.get("coordinate_derived") or source.get("coordinate_truth_used_before_prediction") or source.get("native_metrics_used_for_selection"):
        return COORDINATE_DERIVED
    explicit_noncoordinate = source.get("coordinate_derived") is False
    if not explicit_noncoordinate and any(token in text for token in ["pdb", "mmcif", "coordinate", "contact map", "native contact", "alphafold", "esmfold", "rosettafold"]):
        return COORDINATE_DERIVED
    if source.get("internal_runtime") or source.get("internal_runtime_source"):
        return INTERNAL_RUNTIME
    explicit_not_runtime = source.get("internal_runtime_source") is False and source.get("internal_runtime") is False
    if not explicit_not_runtime and ("first_contact_clean_pharmacotopology_layer_run" in text or "runtime artifact" in text):
        return INTERNAL_RUNTIME
    if source.get("spatial_proxy"):
        return SPATIAL_PROXY_NON_COORDINATE
    if any(token in text for token in ["dca", "fret", "crosslink", "chemical shift", "saxs", "epr", "deer", "distance restraint"]):
        return SPATIAL_PROXY_NON_COORDINATE
    return PURE_NON_COORDINATE


def _evidence_role(source: dict[str, Any]) -> str:
    if source.get("source_role"):
        return str(source["source_role"])
    allowed_use = str(source.get("allowed_use", "")).lower()
    if source.get("holdout_source") or "holdout" in allowed_use or "post" in allowed_use:
        return "holdout_validation"
    if "prediction_input" in allowed_use or source.get("allowed_for_prediction"):
        return "prediction_input"
    if source.get("blocked"):
        return "blocked"
    return "prediction_input"


def evidence_boundary_gate(
    sources: list[dict[str, Any]],
    *,
    allow_spatial_proxy_prediction: bool = True,
    seal_already_created: bool = False,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for source in sources:
        source_class = classify_evidence_source(source)
        role = _evidence_role(source)
        explicitly_tagged_spatial = source.get("source_class") == SPATIAL_PROXY_NON_COORDINATE or source.get("spatial_proxy") is True
        holdout_before_seal = role == "holdout_validation" and not seal_already_created
        allowed_to_initialize = False
        blocked_reason = ""
        if holdout_before_seal:
            blocked_reason = "holdout_opened_before_seal"
        elif source_class == PURE_NON_COORDINATE and role == "prediction_input":
            allowed_to_initialize = True
        elif source_class == SPATIAL_PROXY_NON_COORDINATE and role == "prediction_input" and allow_spatial_proxy_prediction and explicitly_tagged_spatial:
            allowed_to_initialize = True
        elif source_class == SPATIAL_PROXY_NON_COORDINATE:
            blocked_reason = "spatial_proxy_requires_explicit_tag_and_prediction_role"
        elif source_class == COORDINATE_DERIVED:
            blocked_reason = "coordinate_derived_blocked_before_seal"
        elif source_class == INTERNAL_RUNTIME:
            blocked_reason = "internal_runtime_never_biological_evidence"
        else:
            blocked_reason = "not_a_prediction_initialization_source"
        rows.append({
            "source_id": str(source.get("source_id", source.get("accession", "unknown_source"))),
            "source_class": source_class,
            "source_role": role,
            "spatial_proxy": source_class == SPATIAL_PROXY_NON_COORDINATE,
            "coordinate_derived": source_class == COORDINATE_DERIVED,
            "internal_runtime": source_class == INTERNAL_RUNTIME,
            "allowed_to_initialize_language_field": allowed_to_initialize,
            "blocked_reason": blocked_reason,
            "allowed_for_holdout_after_seal": role == "holdout_validation" and source_class != INTERNAL_RUNTIME,
            "allowed_for_claim": allowed_to_initialize or (seal_already_created and role == "holdout_validation" and source_class != INTERNAL_RUNTIME),
        })
    return {
        "kind": "PROTEIN_ESPERANTO_EVIDENCE_BOUNDARY_GATE_v0",
        "rows": rows,
        "allowed_initialization_source_ids": [row["source_id"] for row in rows if row["allowed_to_initialize_language_field"]],
        "blocked_source_ids": [row["source_id"] for row in rows if row["blocked_reason"]],
        "coordinate_derived_source_count_before_prediction": sum(1 for row in rows if row["coordinate_derived"] and row["source_role"] == "prediction_input"),
        "internal_runtime_source_count_for_prediction": sum(1 for row in rows if row["internal_runtime"] and row["source_role"] == "prediction_input"),
        "spatial_proxy_untagged_or_misused_count": sum(1 for row in rows if row["blocked_reason"] == "spatial_proxy_requires_explicit_tag_and_prediction_role"),
        "holdout_opened_before_seal": any(row["blocked_reason"] == "holdout_opened_before_seal" for row in rows),
        "coordinate_truth_used_before_prediction": any(row["coordinate_derived"] and row["source_role"] == "prediction_input" for row in rows),
        "internal_runtime_used_as_biological_evidence": any(row["internal_runtime"] and row["source_role"] == "prediction_input" for row in rows),
    }


def _allowed_source_text(sources: list[dict[str, Any]], gate: dict[str, Any]) -> str:
    allowed = set(gate["allowed_initialization_source_ids"])
    selected = []
    for source in sources:
        source_id = str(source.get("source_id", source.get("accession", "unknown_source")))
        if source_id in allowed:
            visible_source = {
                key: value
                for key, value in source.items()
                if not str(key).startswith("withheld") and str(key) not in {"blocked_prediction_inputs"}
            }
            selected.append(_flatten_text(visible_source))
    return " ".join(selected).lower()


def _contains_any(text: str, tokens: list[str]) -> bool:
    return any(token in text for token in tokens)


def _contains_standalone_word(text: str, word: str) -> bool:
    normalized = "".join(ch if ch.isalnum() else " " for ch in text)
    return word in normalized.split()


def _beta_topology_word_from_text(text: str) -> str | None:
    if any(token in text for token in ["membrane_beta_barrel", "membrane beta barrel", "outer membrane beta barrel"]):
        return "membrane_beta_barrel"
    if any(token in text for token in ["beta_propeller_repeat_closure", "beta propeller", "beta-propeller", "blade repeat", "repeat blade"]):
        return "beta_propeller_repeat_closure"
    if any(token in text for token in ["soluble_beta_barrel", "soluble beta barrel"]):
        return "soluble_beta_barrel"
    if any(token in text for token in ["beta_sandwich_core", "beta sandwich", "immunoglobulin sandwich", "ig-like sandwich"]):
        return "beta_sandwich_core"
    if any(token in text for token in ["jelly_roll_wrap", "jelly roll", "jelly-roll"]):
        return "jelly_roll_wrap"
    if any(token in text for token in ["greek_key_beta_lock", "greek key", "greek-key"]):
        return "greek_key_beta_lock"
    if any(token in text for token in ["beta_helix_solenoid_stack", "beta helix", "beta-helix", "beta solenoid", "beta-solenoid"]):
        return "beta_helix_solenoid_stack"
    if any(token in text for token in ["alpha_beta_barrel_distinction", "alpha beta barrel", "alpha-beta barrel", "tim barrel"]):
        return "alpha_beta_barrel_distinction"
    if any(token in text for token in ["closed_beta_topology", "closed beta topology", "strand_register", "strand register", "beta_sheet_closure", "beta sheet closure"]):
        return "closed_beta_topology"
    return None


def _multidomain_allosteric_word_from_text(text: str) -> str | None:
    if any(token in text for token in ["domain_swapping", "domain swapping", "domain-swapped", "domain swapped", "swapped dimer"]):
        return "domain_swapping"
    if any(token in text for token in ["allosteric_basin_shift", "allosteric basin shift", "allosteric", "allostery"]):
        return "allosteric_basin_shift"
    if any(token in text for token in ["hinge_region", "hinge region", "hinge"]):
        return "hinge_region"
    if any(token in text for token in ["interdomain_lock", "interdomain lock", "inter-domain lock"]):
        return "interdomain_lock"
    if any(token in text for token in ["domain_reorientation", "domain reorientation", "domain movement"]):
        return "domain_reorientation"
    if any(token in text for token in ["modular_architecture", "modular architecture"]):
        return "modular_architecture"
    if any(token in text for token in ["multidomain_allostery", "multidomain allostery", "multidomain", "multi-domain"]):
        return "multidomain_allostery"
    if any(token in text for token in ["domain_boundary", "domain boundary"]):
        return "domain_boundary"
    return None


def select_mechanism_grammar(
    *,
    sequence_field: dict[str, Any],
    evidence_manifest: dict[str, Any],
    sources: list[dict[str, Any]],
    forced_grammar: str | None = None,
) -> dict[str, Any]:
    if not evidence_manifest["allowed_initialization_source_ids"]:
        natural = "insufficient_evidence_clean_abstain"
        reason = "no_allowed_prediction_evidence"
        beta_topology_word = None
        multidomain_word = None
    else:
        text = _allowed_source_text(sources, evidence_manifest)
        metrics = sequence_field["global_metrics"]
        cofactor_context = _contains_any(text, COFACTOR_CONTEXT_TOKENS)
        metal_cluster_context = _contains_any(text, METAL_CLUSTER_CONTEXT_TOKENS)
        ligand_locked_context = _contains_any(text, LIGAND_LOCKED_CONTEXT_TOKENS)
        disorder_boundary_context = _contains_any(text, DISORDER_BOUNDARY_CONTEXT_TOKENS)
        beta_topology_word = _beta_topology_word_from_text(text)
        beta_closure_context = beta_topology_word is not None or _contains_any(text, BETA_CLOSURE_CONTEXT_TOKENS)
        beta_ambiguous_context = _contains_any(text, BETA_AMBIGUOUS_CONTEXT_TOKENS)
        multidomain_word = _multidomain_allosteric_word_from_text(text)
        multidomain_context = multidomain_word is not None or _contains_any(text, MULTIDOMAIN_ALLOSTERIC_CONTEXT_TOKENS)
        negative_assembly_context = _contains_any(text, NEGATIVE_ASSEMBLY_CONTEXT_TOKENS)
        explicit_assembly_required_context = _contains_any(text, ASSEMBLY_REQUIRED_CONTEXT_TOKENS) and not negative_assembly_context
        biological_oligomer_context = _contains_any(text, OLIGOMER_CONTEXT_TOKENS)
        generic_complex_only_context = _contains_any(text, GENERIC_COMPLEX_ONLY_TOKENS) or (
            _contains_standalone_word(text, "complex")
            and not explicit_assembly_required_context
            and not biological_oligomer_context
            and not cofactor_context
        )
        assembly_ambiguity_context = _contains_any(text, ASSEMBLY_AMBIGUITY_TOKENS)
        soluble_monomeric_core_context = any(
            token in text
            for token in [
                "soluble_monomeric_core_context",
                "monomeric soluble core",
                "complete soluble monomer",
                "standalone soluble fold",
            ]
        )
        negative_membrane_topology_context = _contains_any(text, NEGATIVE_MEMBRANE_TOPOLOGY_TOKENS)
        peripheral_membrane_context = _contains_any(text, PERIPHERAL_MEMBRANE_CONTEXT_TOKENS)
        topology_conflict_context = negative_membrane_topology_context or peripheral_membrane_context
        explicit_membrane_context = _contains_any(text, STRONG_MEMBRANE_CONTEXT_TOKENS) and not topology_conflict_context
        membrane_text_context = any(
            token in text
            for token in ["cftr", "f508", "proteostasis", "trafficking", "nbd1", "transmembrane", "membrane"]
        )
        if (
            any(token in text for token in ["orf6", "rae1", "nup98", "host hijack", "nucleocytoplasmic", "interferon"])
            and not explicit_membrane_context
            and not explicit_assembly_required_context
            and not metal_cluster_context
            and not ligand_locked_context
        ):
            natural = "short_region_host_interface_hijacking"
            reason = "short_region_host_interface_evidence"
        elif any(token in text for token in [
            "rfah",
            "fold switch",
            "fold-switch",
            "metamorphic",
            "autoinhibited",
            "released ctd",
            "alpha-state",
            "beta-state",
            "alpha state",
            "beta state",
            "dual basin",
            "dual-basin",
            "competing basin",
        ]):
            natural = "metamorphic_fold_switching"
            reason = "state_separated_partner_context_evidence"
        elif beta_topology_word == "membrane_beta_barrel":
            natural = "beta_closure_topology"
            reason = "explicit_membrane_beta_barrel_closure_context_preserves_true_membrane_beta_topology"
        elif explicit_membrane_context:
            if (
                "f508" in text
                or "proteostasis" in text
                or "trafficking" in text
                or metrics["mean_membrane"] >= 0.10
                or explicit_membrane_context
            ):
                natural = "membrane_multidomain_folding_proteostasis"
                reason = "strong_membrane_topology_context_prioritized_over_incidental_ligand_or_assembly_context"
            else:
                natural = "insufficient_evidence_clean_abstain"
                reason = "generic_membrane_annotation_without_specific_operator"
        elif explicit_assembly_required_context:
            natural = "assembly_required_folding"
            reason = "explicit_assembly_required_core_or_partner_completed_topology_context"
        elif metal_cluster_context or ligand_locked_context:
            natural = "metal_cluster_and_ligand_locked_basin"
            reason = "explicit_metal_cluster_geometry_or_ligand_locked_basin_context"
        elif disorder_boundary_context:
            natural = "disorder_boundary_and_fold_upon_binding"
            reason = "explicit_disorder_boundary_or_fold_upon_binding_context_prioritized_over_generic_oligomer"
        elif beta_ambiguous_context:
            natural = "insufficient_evidence_clean_abstain"
            reason = "beta_propensity_or_open_sheet_without_closure_logic_requires_abstention"
        elif beta_closure_context:
            natural = "beta_closure_topology"
            reason = "explicit_beta_closure_topology_context"
        elif multidomain_context:
            natural = "multidomain_allosteric_architecture"
            reason = "explicit_multidomain_allosteric_architecture_context"
        elif cofactor_context:
            natural = "cofactor_ligand_assisted_stabilization"
            reason = "explicit_ligand_cofactor_or_metal_context"
        elif soluble_monomeric_core_context and metrics["mean_disorder"] < 0.32:
            natural = "globular_closure"
            reason = "explicit_soluble_monomeric_core_context"
        elif biological_oligomer_context and (
            assembly_ambiguity_context
            or peripheral_membrane_context
            or generic_complex_only_context
            or (membrane_text_context and not topology_conflict_context)
            or cofactor_context
        ):
            natural = "insufficient_evidence_clean_abstain"
            reason = "assembly_topology_ambiguous_competing_weak_explanations"
        elif biological_oligomer_context:
            natural = "oligomerization_controlled_folding"
            reason = "explicit_biological_oligomer_context_without_assembly_required_claim"
        elif topology_conflict_context:
            natural = "insufficient_evidence_clean_abstain"
            reason = "membrane_topology_conflict_or_peripheral_context_requires_abstention"
        elif generic_complex_only_context:
            natural = "insufficient_evidence_clean_abstain"
            reason = "generic_complex_only_is_not_obligate_assembly_evidence"
        elif membrane_text_context:
            if (
                "f508" in text
                or "proteostasis" in text
                or "trafficking" in text
                or metrics["mean_membrane"] >= 0.10
            ):
                natural = "membrane_multidomain_folding_proteostasis"
                reason = "membrane_multidomain_mutation_interface_evidence"
            else:
                natural = "insufficient_evidence_clean_abstain"
                reason = "generic_membrane_annotation_without_specific_operator"
        elif any(token in text for token in ["low complexity", "disprot", "disordered", "phase", "prion", "lcd", "fus", "tdp-43", "tdp43"]):
            if metrics["low_complexity_density"] >= 0.25 or metrics["mean_disorder"] >= 0.18:
                natural = "intrinsic_disorder_phase_separation"
                reason = "low_complexity_disorder_phase_evidence"
            else:
                natural = "insufficient_evidence_clean_abstain"
                reason = "disorder_text_without_sequence_support"
        elif any(token in text for token in ["binding", "partner", "motif"]) and metrics["mean_disorder"] >= 0.18:
            natural = "fold_upon_binding_disorder"
            reason = "partner_conditioned_disorder_evidence"
        elif any(token in text for token in ["cofactor", "ligand"]):
            natural = "cofactor_ligand_assisted_stabilization"
            reason = "ligand_or_cofactor_context"
        elif any(token in text for token in ["oligomer", "assembly", "multimer"]):
            if _contains_standalone_word(text, "complex") or membrane_text_context:
                natural = "insufficient_evidence_clean_abstain"
                reason = "assembly_ambiguous_generic_text_requires_abstention"
            else:
                natural = "oligomerization_controlled_folding"
                reason = "oligomerization_context"
        elif any(token in text for token in [
            "folding kinetics",
            "folding-kinetics",
            "single-domain",
            "single domain",
            "protein folding target",
            "three-helix",
            "helix bundle",
            "beta-sheet",
            "beta sheet",
            "ww domain",
            "sh3",
            "mini-protein",
            "miniprotein",
            "fast folder",
            "fast mini",
            "caged aromatic core",
        ]) and metrics["mean_disorder"] < 0.32:
            natural = "globular_closure"
            reason = "noncoordinate_process_metadata_and_sequence_support_globular_closure"
        elif any(token in text for token in ["mini-protein", "miniprotein", "trp-cage", "caged aromatic core"]):
            natural = "globular_closure"
            reason = "small_fast_folder_metadata_supports_globular_closure"
        elif metrics["hydrophobic_density"] >= 0.32 and metrics["mean_disorder"] < 0.20:
            natural = "globular_closure"
            reason = "sequence_only_hydrophobic_closure_support"
        else:
            natural = "insufficient_evidence_clean_abstain"
            reason = "insufficient_specific_operator_evidence"
    forced_rejected = forced_grammar is not None and forced_grammar != natural
    mechanism_class = "insufficient_evidence_clean_abstain" if forced_rejected else natural
    return {
        "kind": "PROTEIN_ESPERANTO_MECHANISM_CLASSIFIER_v0",
        "mechanism_class": mechanism_class,
        "natural_mechanism_class": natural,
        "selected_beta_topology_word": beta_topology_word if mechanism_class == "beta_closure_topology" else None,
        "selected_multidomain_word": multidomain_word if mechanism_class == "multidomain_allosteric_architecture" else None,
        "forced_grammar": forced_grammar,
        "forced_grammar_rejected": forced_rejected,
        "selection_reason": "forced_wrong_grammar_rejected" if forced_rejected else reason,
        "grammar_rule": GRAMMAR_RULES.get(mechanism_class),
        "available_mechanism_classes": MECHANISM_CLASSES,
        "folding_problem_solved": False,
    }


def _strong_segments(sequence_field: dict[str, Any], key: str, *, limit: int = 3) -> list[dict[str, Any]]:
    return sorted(sequence_field["segments"], key=lambda row: row.get(key, 0.0), reverse=True)[:limit]


def _span_from_segment(segment: dict[str, Any]) -> str:
    return f"{segment['start']}-{segment['end']}"


def _operator(
    name: str,
    acts_on: str,
    strength: float,
    evidence: list[str],
    pushes: str,
    perturbation: str,
    falsifier: str,
    state_variable: str,
) -> dict[str, Any]:
    return {
        "operator": name,
        "acts_on": acts_on,
        "activation_strength": bounded(strength),
        "activated_by_evidence_ids": evidence,
        "pushes_toward": pushes,
        "perturbation_should": perturbation,
        "falsified_by": falsifier,
        "state_variable": state_variable,
    }


def build_operator_field(
    *,
    sequence_field: dict[str, Any],
    mechanism: dict[str, Any],
    evidence_manifest: dict[str, Any],
    focus_regions: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    mechanism_class = mechanism["mechanism_class"]
    evidence = list(evidence_manifest["allowed_initialization_source_ids"])
    metrics = sequence_field["global_metrics"]
    focus_regions = focus_regions or []
    operators: list[dict[str, Any]] = []
    hydrophobic_segments = _strong_segments(sequence_field, "hydrophobic_density")
    disorder_segments = _strong_segments(sequence_field, "disorder_density")
    phase_segments = sorted(
        sequence_field["segments"],
        key=lambda row: row["low_complexity_density"] + row["aromatic_density"] + abs(row["charge_density"]),
        reverse=True,
    )[:3]
    interface_segments = _strong_segments(sequence_field, "interface_density")
    membrane_segments = _strong_segments(sequence_field, "membrane_density")
    if mechanism_class == "globular_closure":
        operators.append(_operator(
            "closure_operator",
            ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
            metrics["hydrophobic_density"] + 0.25,
            evidence,
            "local compaction and contact probability increase",
            "core hydrophobic/aromatic substitutions weaken compaction",
            "no compaction under independent SAXS/FRET/crosslink/structure holdout",
            "segment_compaction",
        ))
    elif mechanism_class == "intrinsic_disorder_phase_separation":
        operators.extend([
            _operator(
                "disorder_operator",
                ", ".join(_span_from_segment(row) for row in disorder_segments),
                metrics["mean_disorder"] + metrics["low_complexity_density"],
                evidence,
                "expanded dynamic ensemble with no stable single fold",
                "charge/proline/glycine changes tune expansion",
                "stable globular fold dominates without condition dependence",
                "disorder_order_balance",
            ),
            _operator(
                "phase_operator",
                ", ".join(_span_from_segment(row) for row in phase_segments),
                metrics["low_complexity_density"] + metrics["aromatic_density"] + 0.20,
                evidence,
                "weak multivalent attraction and phase-prone condensate threshold",
                "aromatic/charge/salt/RNA perturbations shift phase threshold",
                "phase behavior unaffected by predicted stickers or conditions",
                "state_basin_occupancy",
            ),
            _operator(
                "repulsion_operator",
                "charged low-complexity windows",
                abs(metrics["net_charge_per_residue"]) + 0.25,
                evidence,
                "expanded ensemble pressure balancing weak attraction",
                "charge screening changes ensemble size",
                "charge/salt manipulations have no directional effect",
                "residue_exposure",
            ),
        ])
    elif mechanism_class == "disorder_boundary_and_fold_upon_binding":
        operators.extend([
            _operator(
                "disorder_operator",
                ", ".join(_span_from_segment(row) for row in disorder_segments),
                metrics["mean_disorder"] + metrics["low_complexity_density"] + 0.18,
                evidence,
                "persistent IDR boundary and structured-domain/tail separation",
                "boundary truncation or charge/proline/glycine edits shift disorder persistence",
                "generic oligomer interface explains the holdout without an IDR boundary",
                "IDR_boundary",
            ),
            _operator(
                "interface_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments),
                metrics["mean_interface"] + metrics["mean_disorder"] + 0.20,
                evidence,
                "local fold-upon-binding motif readiness without whole-chain compaction",
                "partner or motif removal weakens local ordering",
                "partner context has no directional effect on the predicted local motif",
                "fold_upon_binding_region",
            ),
            _operator(
                "phase_operator",
                ", ".join(_span_from_segment(row) for row in phase_segments),
                metrics["low_complexity_density"] + metrics["aromatic_density"] + 0.22,
                evidence,
                "phase-prone low-complexity pressure bounded by local motif or domain context",
                "sticker/charge/salt perturbations shift phase-prone basin",
                "low-complexity context behaves as an ordinary compact loop",
                "phase_prone_low_complexity",
            ),
        ])
    elif mechanism_class == "beta_closure_topology":
        beta_word = mechanism.get("selected_beta_topology_word") or "closed_beta_topology"
        operators.extend([
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in _strong_segments(sequence_field, "beta_propensity_density")),
                metrics["hydrophobic_density"] + metrics["beta_propensity_density"] + 0.18,
                evidence,
                f"{beta_word} strand-register closure rather than generic beta propensity",
                "register-shift or blade/repeat perturbation weakens closed beta topology",
                "open beta sheet or wrong beta subtype explains the holdout better",
                beta_word,
            ),
            _operator(
                "interface_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments),
                metrics["mean_interface"] + 0.30,
                evidence,
                "inter-strand edge pairing and closure interface readiness",
                "edge/interface perturbation opens the beta topology",
                "beta-rich region stays open without closure contacts",
                "beta_sheet_closure",
            ),
            _operator(
                "frustration_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                metrics["aromatic_density"] + metrics["mean_interface"] + 0.14,
                evidence,
                "closure conflict separates membrane, propeller, sandwich, solenoid, and alpha/beta barrel classes",
                "wrong closure class increases topology conflict",
                "all beta-rich classes collapse to the same topology",
                "beta_topology_conflict",
            ),
        ])
    elif mechanism_class == "multidomain_allosteric_architecture":
        multidomain_word = mechanism.get("selected_multidomain_word") or "multidomain_allostery"
        operators.extend([
            _operator(
                "interface_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments),
                metrics["mean_interface"] + 0.42,
                evidence,
                "interdomain lock and domain-boundary interface readiness",
                "interface or lock perturbation weakens interdomain coupling",
                "single-domain closure explains the holdout without interdomain dependence",
                "interdomain_lock",
            ),
            _operator(
                "dual_basin_switch_operator",
                "hinge/allosteric domain-coupling axis",
                metrics["mean_interface"] + metrics["aromatic_density"] + 0.30,
                evidence,
                "domain reorientation shifts an allosteric basin rather than one averaged fold",
                "hinge or allosteric push changes basin occupancy and orientation",
                "domain orientation is condition-independent",
                "allosteric_basin_shift",
            ),
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                metrics["hydrophobic_density"] + 0.18,
                evidence,
                "modular architecture closes through coupled domain packing",
                "domain-boundary or packing damage lowers modular compaction",
                "a standalone compact domain explains the full architecture",
                "modular_architecture",
            ),
            _operator(
                "frustration_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments[:2]),
                metrics["mean_interface"] + metrics["aromatic_density"] + 0.12,
                evidence,
                f"{multidomain_word} separates hinge, lock, reorientation, and swapped-domain topology from generic domains",
                "wrong domain-boundary operator raises allosteric or swapping conflict",
                "generic domain text alone predicts all states equally",
                multidomain_word,
            ),
        ])
    elif mechanism_class == "membrane_multidomain_folding_proteostasis":
        f508 = next((row for row in focus_regions if "F508" in str(row.get("name", "")) or row.get("position") == 508), None)
        f508_span = "508" if f508 is None else str(f508.get("span", f508.get("position", "508")))
        operators.extend([
            _operator(
                "membrane_pressure_operator",
                ", ".join(_span_from_segment(row) for row in membrane_segments),
                metrics["mean_membrane"] + 0.35,
                evidence,
                "membrane-buried routing pressure and topology context",
                "membrane/context disruption weakens maturation route",
                "membrane context irrelevant to holdout rescue logic",
                "proteostasis_routing",
            ),
            _operator(
                "closure_operator",
                f"NBD1/local deletion focus {f508_span}",
                0.70,
                evidence,
                "NBD1 local stability and partial domain closure",
                "F508del weakens; NBD1-only correction partially rescues",
                "NBD1 stability evidence fails to track F508del",
                "segment_compaction",
            ),
            _operator(
                "interface_operator",
                "NBD1-MSD/interdomain correction axis",
                0.74,
                evidence,
                "interdomain interface readiness",
                "interface correction strengthens rescue beyond NBD1-only correction",
                "interface/proteostasis correction has no added effect",
                "interface_readiness",
            ),
            _operator(
                "proteostasis_operator",
                "folding quality-control and trafficking route",
                0.78,
                evidence,
                "maturation and trafficking escape from quality control",
                "corrector/proteostasis condition improves routing",
                "maturation defect absent despite local destabilization",
                "proteostasis_routing",
            ),
        ])
    elif mechanism_class == "metamorphic_fold_switching":
        operators.extend([
            _operator(
                "dual_basin_switch_operator",
                "state-separated CTD/core region",
                0.86,
                evidence,
                "two incompatible alpha/beta state basins",
                "partner/release perturbations shift alpha/beta occupancy",
                "only one stable state under all contexts",
                "state_basin_occupancy",
            ),
            _operator(
                "interface_operator",
                "partner-context interface",
                0.72,
                evidence,
                "context-conditioned state selection",
                "partner removal or stabilization shifts state occupancy",
                "partner context has no directional effect",
                "interface_readiness",
            ),
            _operator(
                "frustration_operator",
                "secondary-structure conflict region",
                0.68,
                evidence,
                "state conflict rather than averaged consensus fold",
                "forcing one fold fails",
                "one averaged state explains all holdouts",
                "frustration",
            ),
        ])
    elif mechanism_class == "short_region_host_interface_hijacking":
        length = sequence_field["sequence_length"]
        cterm = f"{max(1, length - 23)}-{length}"
        operators.extend([
            _operator(
                "host_hijack_operator",
                f"C-terminal host-interface region {cterm}",
                0.90,
                evidence,
                "host-interface capture without global fold requirement",
                "C-terminal disruption weakens host binding and transport/IFN consequences",
                "wrong-region perturbations dominate while C-terminus is inert",
                "interface_readiness",
            ),
            _operator(
                "interface_operator",
                f"short linear motif/interface window {cterm}",
                0.78,
                evidence,
                "RAE1/NUP98 or host-surface readiness",
                "host partner removal weakens the operator",
                "host partner evidence does not localize to predicted region",
                "contact_probability",
            ),
            _operator(
                "disorder_operator",
                "short exposed accessory protein context",
                max(metrics["mean_disorder"], 0.35),
                evidence,
                "exposed motif rather than stable globular fold",
                "forcing compact fold is rejected",
                "global stable fold is required for function",
                "residue_exposure",
            ),
        ])
    elif mechanism_class == "cofactor_ligand_assisted_stabilization":
        operators.extend([
            _operator(
                "interface_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments),
                metrics["mean_interface"] + 0.32,
                evidence,
                "ligand/cofactor pocket readiness and local stabilization",
                "cofactor removal or pocket disruption weakens the stabilized basin",
                "apo and ligand-bound contexts are indistinguishable",
                "interface_readiness",
            ),
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                metrics["hydrophobic_density"] + 0.18,
                evidence,
                "cofactor-assisted compaction of a weak apo basin",
                "cofactor pocket perturbation reduces compaction",
                "ligand context has no directional effect on stability",
                "segment_compaction",
            ),
        ])
    elif mechanism_class == "metal_cluster_and_ligand_locked_basin":
        operators.extend([
            _operator(
                "interface_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments),
                metrics["mean_interface"] + 0.42,
                evidence,
                "coordination shell and ligand pocket complete the locked basin",
                "metal/cofactor removal or coordinating-side-chain disruption unlocks the basin",
                "generic cofactor annotation explains the holdout without a geometry-locked basin",
                "coordination_shell_integrity",
            ),
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                metrics["hydrophobic_density"] + 0.20,
                evidence,
                "ligand-locked compaction rather than free apo closure",
                "apo conversion lowers basin occupancy and compaction",
                "ligand-free closure remains equally stable",
                "ligand_locked_basin",
            ),
            _operator(
                "frustration_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments[:2]),
                metrics["mean_interface"] + metrics["aromatic_density"] + 0.16,
                evidence,
                "apo/holo separation and coordination geometry resolve local frustration",
                "coordination loss increases frustration and shifts basin occupancy",
                "metal or ligand context has no directional effect",
                "apo_holo_basin_shift",
            ),
        ])
    elif mechanism_class == "assembly_required_folding":
        operators.extend([
            _operator(
                "interface_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments),
                metrics["mean_interface"] + 0.42,
                evidence,
                "partner-completed core and biological assembly readiness",
                "partner/interface disruption exposes incomplete monomer topology",
                "complete monomer, membrane, or ligand grammar explains the holdout without assembly",
                "partner_completed_core",
            ),
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                metrics["hydrophobic_density"] + 0.16,
                evidence,
                "hydrophobic core closure only after assembly context resolves the surface",
                "assembly-interface mutation lowers partner-completed closure",
                "monomer-only closure remains stable and complete",
                "assembly_required_core",
            ),
            _operator(
                "frustration_operator",
                ", ".join(_span_from_segment(row) for row in membrane_segments),
                metrics["mean_membrane"] + metrics["mean_interface"] + 0.10,
                evidence,
                "unresolved monomer topology rather than clean membrane insertion",
                "partner completion lowers topological frustration",
                "true membrane topology or soluble monomer explains the signal directly",
                "monomer_incomplete_topology",
            ),
        ])
    elif mechanism_class == "oligomerization_controlled_folding":
        operators.extend([
            _operator(
                "interface_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments),
                metrics["mean_interface"] + 0.38,
                evidence,
                "partner-copy interface readiness and assembly-stabilized folding",
                "interface or concentration perturbation weakens the assembled basin",
                "monomer-only grammar explains the holdout equally well",
                "interface_readiness",
            ),
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                metrics["hydrophobic_density"] + 0.12,
                evidence,
                "local closure coupled to oligomeric surface burial",
                "assembly-interface mutation lowers closure support",
                "folding is independent of oligomeric context",
                "segment_compaction",
            ),
        ])
    elif mechanism_class != "insufficient_evidence_clean_abstain":
        operators.append(_operator(
            "interface_operator",
            ", ".join(_span_from_segment(row) for row in interface_segments),
            metrics["mean_interface"] + 0.20,
            evidence,
            "context-conditioned interaction readiness",
            "context perturbation shifts mechanism",
            "context has no directional effect",
            "interface_readiness",
        ))
    return {
        "kind": "PROTEIN_ESPERANTO_OPERATOR_FIELD_v0",
        "mechanism_class": mechanism_class,
        "operators": operators,
        "operator_names": [row["operator"] for row in operators],
        "active_operator_count": len(operators),
        "coordinate_truth_used": False,
    }


def _operator_strength(field: dict[str, Any], name: str) -> float:
    values = [row["activation_strength"] for row in field["operators"] if row["operator"] == name]
    return max(values) if values else 0.0


def _with_scaled_operators(operator_field: dict[str, Any], scales: dict[str, float]) -> dict[str, Any]:
    field = deepcopy(operator_field)
    for row in field["operators"]:
        row["activation_strength"] = bounded(row["activation_strength"] * scales.get(row["operator"], 1.0))
    return field


def _jitter(seed: str, step: int) -> float:
    digest = int(stable_hash({"seed": seed, "step": step})[:8], 16)
    return round((sin(digest) + 1.0) * 0.01, 6)


def simulate_operator_trajectory(
    *,
    sequence_field: dict[str, Any],
    operator_field: dict[str, Any],
    mechanism_class: str,
    perturbation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if perturbation:
        operator_field = _with_scaled_operators(operator_field, perturbation.get("operator_scales", {}))
    seed = stable_hash({
        "length": sequence_field["sequence_length"],
        "mechanism": mechanism_class,
        "operators": operator_field["operators"],
        "perturbation": perturbation or {},
    })
    strengths = {name: _operator_strength(operator_field, name) for name in UNIVERSAL_OPERATORS}
    metrics = sequence_field["global_metrics"]
    beta_subtype = next(
        (
            row["state_variable"]
            for row in operator_field["operators"]
            if row["state_variable"] in BETA_TOPOLOGY_STATE_VARIABLES
        ),
        "closed_beta_topology",
    )
    multidomain_state_variables = [
        row["state_variable"]
        for row in operator_field["operators"]
        if row["state_variable"] in MULTIDOMAIN_ALLOSTERIC_STATE_VARIABLES
    ]
    multidomain_subtype = (
        "domain_swapping"
        if "domain_swapping" in multidomain_state_variables
        else (multidomain_state_variables[-1] if multidomain_state_variables else "multidomain_allostery")
    )
    timepoints: list[dict[str, Any]] = []
    for timepoint in [0, 100, 500, 1000]:
        progress = timepoint / 1000.0
        noise = _jitter(seed, timepoint)
        closure = strengths["closure_operator"]
        disorder = strengths["disorder_operator"]
        phase = strengths["phase_operator"]
        membrane = strengths["membrane_pressure_operator"]
        switch = strengths["dual_basin_switch_operator"]
        proteostasis = strengths["proteostasis_operator"]
        host = strengths["host_hijack_operator"]
        interface = strengths["interface_operator"]
        frustration_strength = strengths["frustration_operator"]
        if mechanism_class == "intrinsic_disorder_phase_separation":
            basin = {
                "expanded_disordered": bounded(0.62 - 0.10 * progress + disorder * 0.10),
                "phase_prone_dynamic": bounded(0.22 + phase * 0.36 * progress),
                "compact_single_fold": bounded(0.04 + closure * 0.05),
            }
            segment_compaction = bounded(0.10 + phase * 0.18 * progress + closure * 0.05 + noise)
            contact_probability = bounded(0.08 + phase * 0.16 * progress + noise)
            interface_readiness = bounded(0.10 + interface * 0.15)
            proteostasis_routing = 0.0
        elif mechanism_class == "disorder_boundary_and_fold_upon_binding":
            partner_loss = float((perturbation or {}).get("partner_loss", 0.0))
            motif_damage = float((perturbation or {}).get("motif_damage", 0.0))
            idr_boundary = bounded(0.34 + disorder * 0.42 + phase * 0.12 - 0.26 * motif_damage)
            local_order = bounded(0.18 + interface * 0.44 * progress + disorder * 0.14 - 0.44 * partner_loss - 0.32 * motif_damage)
            low_complexity_basin = bounded(0.22 + phase * 0.38 * progress + disorder * 0.10)
            basin = {
                "disorder_boundary_ensemble": idr_boundary,
                "fold_upon_binding_basin": local_order,
                "phase_prone_low_complexity": low_complexity_basin,
                "compact_single_fold": bounded(0.05 + closure * 0.04),
            }
            segment_compaction = bounded(0.10 + local_order * 0.18 + closure * 0.04 + noise)
            contact_probability = bounded(0.10 + local_order * 0.32 + phase * 0.08)
            interface_readiness = local_order
            proteostasis_routing = 0.0
        elif mechanism_class == "beta_closure_topology":
            register_damage = float((perturbation or {}).get("register_damage", 0.0))
            closure_conflict = float((perturbation or {}).get("closure_conflict", 0.0))
            register = bounded(0.28 + closure * 0.34 * progress + interface * 0.20 - 0.45 * register_damage)
            closure_state = bounded(
                0.24
                + closure * 0.36 * progress
                + interface * 0.22
                + metrics["beta_propensity_density"] * 0.12
                - 0.34 * closure_conflict
            )
            conflict = bounded(0.14 + frustration_strength * 0.20 + closure_conflict + register_damage * 0.25)
            basin = {
                "closed_beta_topology": closure_state,
                "strand_register": register,
                "open_beta_sheet_ambiguous": bounded(0.46 - closure_state * 0.28 + conflict * 0.18),
                "wrong_beta_topology": conflict,
            }
            segment_compaction = bounded(0.16 + closure_state * 0.40 + closure * 0.12 + noise)
            contact_probability = bounded(0.14 + closure_state * 0.42 + register * 0.18)
            interface_readiness = bounded(0.16 + interface * 0.34 + register * 0.20)
            proteostasis_routing = 0.0
        elif mechanism_class == "multidomain_allosteric_architecture":
            hinge_damage = float((perturbation or {}).get("hinge_damage", 0.0))
            lock_damage = float((perturbation or {}).get("lock_damage", 0.0))
            allosteric_push = float((perturbation or {}).get("allosteric_push", 0.0))
            domain_boundary_state = bounded(
                0.24
                + interface * 0.20
                + frustration_strength * 0.12
                + metrics["mean_interface"] * 0.10
                + progress * 0.08
                - hinge_damage * 0.12
            )
            interdomain_lock = bounded(
                0.18
                + interface * 0.44 * progress
                + closure * 0.14
                + domain_boundary_state * 0.18
                - 0.42 * lock_damage
            )
            allosteric_shift = bounded(
                0.14
                + switch * 0.44 * progress
                + frustration_strength * 0.14
                + allosteric_push
                - 0.18 * lock_damage
            )
            reorientation = bounded(
                0.12
                + switch * 0.34 * progress
                + interface * 0.20
                + allosteric_shift * 0.16
                - 0.30 * hinge_damage
            )
            hinge = bounded(
                0.20
                + frustration_strength * 0.24
                + domain_boundary_state * 0.18
                + switch * 0.08
                - 0.25 * hinge_damage
            )
            modular = bounded(0.18 + closure * 0.30 * progress + interface * 0.18 + domain_boundary_state * 0.16)
            swapped = bounded(interdomain_lock * 0.86 if multidomain_subtype == "domain_swapping" else max(0.0, interdomain_lock - 0.16) * 0.20)
            basin = {
                "multidomain_allosteric_basin": bounded(max(allosteric_shift, reorientation, interdomain_lock)),
                "domain_boundary": domain_boundary_state,
                "hinge_region": hinge,
                "interdomain_locked_basin": interdomain_lock,
                "allosteric_basin_shift": allosteric_shift,
                "domain_reorientation_basin": reorientation,
                "modular_architecture": modular,
                "domain_swapping": swapped,
                "single_domain_shortcut": bounded(0.46 - modular * 0.20 - interdomain_lock * 0.20 + lock_damage * 0.12),
            }
            segment_compaction = bounded(0.18 + modular * 0.30 + closure * 0.20 * progress + noise)
            contact_probability = bounded(0.16 + interdomain_lock * 0.36 + allosteric_shift * 0.16 + closure * 0.12)
            interface_readiness = interdomain_lock
            proteostasis_routing = 0.0
        elif mechanism_class == "membrane_multidomain_folding_proteostasis":
            damage = float((perturbation or {}).get("damage", 0.0))
            rescue = float((perturbation or {}).get("rescue", 0.0))
            stability = bounded(0.55 + closure * 0.22 * progress - 0.30 * damage + 0.20 * rescue)
            interface_ready = bounded(0.45 + interface * 0.25 * progress - 0.22 * damage + 0.24 * rescue)
            routing = bounded(0.40 + proteostasis * 0.28 * progress + membrane * 0.10 - 0.30 * damage + 0.30 * rescue)
            basin = {
                "mature_membrane_routed": routing,
                "qc_retained_misfolded": bounded(0.62 - routing + 0.30 * damage - 0.15 * rescue),
                "partial_nbd1_rescue": stability,
            }
            segment_compaction = stability
            contact_probability = bounded(0.20 + interface_ready * 0.45)
            interface_readiness = interface_ready
            proteostasis_routing = routing
        elif mechanism_class == "metamorphic_fold_switching":
            release = float((perturbation or {}).get("release", 0.0))
            alpha_bias = float((perturbation or {}).get("alpha_bias", 0.0))
            beta_bias = float((perturbation or {}).get("beta_bias", 0.0))
            alpha = bounded(0.64 - 0.24 * progress * switch - 0.20 * release + 0.22 * alpha_bias)
            beta = bounded(0.22 + 0.24 * progress * switch + 0.22 * release + 0.22 * beta_bias)
            basin = {
                "alpha_context_basin": alpha,
                "beta_released_basin": beta,
                "averaged_single_fold": bounded(0.06 + (1.0 - switch) * 0.06),
            }
            segment_compaction = bounded(0.36 + 0.10 * switch)
            contact_probability = bounded(0.28 + 0.12 * switch)
            interface_readiness = bounded(0.35 + interface * 0.30)
            proteostasis_routing = 0.0
        elif mechanism_class == "short_region_host_interface_hijacking":
            disruption = float((perturbation or {}).get("interface_disruption", 0.0))
            host_ready = bounded(0.20 + host * 0.65 * progress + interface * 0.20 - 0.55 * disruption)
            basin = {
                "host_interface_engaged": host_ready,
                "exposed_short_region": bounded(0.55 + disorder * 0.20 - 0.20 * disruption),
                "compact_single_fold": bounded(0.04 + closure * 0.04),
            }
            segment_compaction = bounded(0.08 + closure * 0.07)
            contact_probability = bounded(0.12 + host_ready * 0.52)
            interface_readiness = host_ready
            proteostasis_routing = 0.0
        elif mechanism_class == "globular_closure":
            basin = {
                "compact_folded": bounded(0.18 + closure * 0.62 * progress),
                "expanded_unfolded": bounded(0.78 - closure * 0.50 * progress),
            }
            segment_compaction = bounded(0.18 + closure * 0.56 * progress)
            contact_probability = bounded(0.14 + closure * 0.62 * progress)
            interface_readiness = bounded(0.12 + interface * 0.18)
            proteostasis_routing = 0.0
        elif mechanism_class == "cofactor_ligand_assisted_stabilization":
            cofactor_loss = float((perturbation or {}).get("cofactor_loss", 0.0))
            rescue = float((perturbation or {}).get("rescue", 0.0))
            pocket_ready = bounded(0.24 + interface * 0.52 * progress + closure * 0.14 - 0.45 * cofactor_loss + 0.22 * rescue)
            compact = bounded(0.22 + closure * 0.42 * progress + pocket_ready * 0.16 - 0.24 * cofactor_loss)
            basin = {
                "ligand_stabilized_basin": pocket_ready,
                "apo_weak_basin": bounded(0.62 - pocket_ready + 0.28 * cofactor_loss - 0.10 * rescue),
                "generic_compact_basin": bounded(0.12 + closure * 0.10),
            }
            segment_compaction = compact
            contact_probability = bounded(0.18 + compact * 0.35 + pocket_ready * 0.18)
            interface_readiness = pocket_ready
            proteostasis_routing = 0.0
        elif mechanism_class == "metal_cluster_and_ligand_locked_basin":
            cofactor_loss = float((perturbation or {}).get("cofactor_loss", 0.0))
            coordination_damage = float((perturbation or {}).get("coordination_damage", 0.0))
            rescue = float((perturbation or {}).get("rescue", 0.0))
            coordination_shell = bounded(
                0.18
                + interface * 0.50 * progress
                + frustration_strength * 0.16
                - 0.50 * coordination_damage
                - 0.34 * cofactor_loss
                + 0.20 * rescue
            )
            locked_basin = bounded(
                0.20
                + closure * 0.34 * progress
                + coordination_shell * 0.34
                - 0.42 * cofactor_loss
                + 0.16 * rescue
            )
            basin = {
                "metal_cluster_locked_basin": coordination_shell,
                "ligand_locked_basin": locked_basin,
                "apo_unlocked_basin": bounded(0.64 - locked_basin + 0.26 * cofactor_loss + 0.22 * coordination_damage),
            }
            segment_compaction = bounded(0.18 + closure * 0.38 * progress + locked_basin * 0.22)
            contact_probability = bounded(0.16 + segment_compaction * 0.30 + coordination_shell * 0.26)
            interface_readiness = coordination_shell
            proteostasis_routing = 0.0
        elif mechanism_class == "assembly_required_folding":
            interface_disruption = float((perturbation or {}).get("interface_disruption", 0.0))
            concentration_rescue = float((perturbation or {}).get("concentration_rescue", 0.0))
            partner_completion = bounded(
                0.16
                + interface * 0.58 * progress
                + closure * 0.18
                + frustration_strength * 0.08
                - 0.54 * interface_disruption
                + 0.24 * concentration_rescue
            )
            monomer_gap = bounded(0.52 + frustration_strength * 0.18 - partner_completion * 0.34 + 0.28 * interface_disruption)
            basin = {
                "assembly_required_basin": partner_completion,
                "monomer_incomplete_topology": monomer_gap,
                "assembly_ambiguous_basin": bounded(0.24 + 0.24 * interface_disruption - 0.16 * partner_completion),
            }
            segment_compaction = bounded(0.18 + closure * 0.34 * progress + partner_completion * 0.22)
            contact_probability = bounded(0.14 + partner_completion * 0.46 + closure * 0.14)
            interface_readiness = partner_completion
            proteostasis_routing = 0.0
        elif mechanism_class == "oligomerization_controlled_folding":
            interface_disruption = float((perturbation or {}).get("interface_disruption", 0.0))
            concentration_rescue = float((perturbation or {}).get("concentration_rescue", 0.0))
            assembly_ready = bounded(0.18 + interface * 0.62 * progress + closure * 0.10 - 0.48 * interface_disruption + 0.22 * concentration_rescue)
            basin = {
                "assembly_stabilized_basin": assembly_ready,
                "monomer_partial_order": bounded(0.46 + closure * 0.18 - assembly_ready * 0.22),
                "interface_rejected_basin": bounded(0.10 + 0.35 * interface_disruption),
            }
            segment_compaction = bounded(0.20 + closure * 0.36 * progress + assembly_ready * 0.18)
            contact_probability = bounded(0.16 + assembly_ready * 0.42 + closure * 0.16)
            interface_readiness = assembly_ready
            proteostasis_routing = 0.0
        else:
            basin = {"clean_abstain": 1.0}
            segment_compaction = 0.0
            contact_probability = 0.0
            interface_readiness = 0.0
            proteostasis_routing = 0.0
        exposure = bounded(0.70 + metrics["mean_disorder"] * 0.18 - segment_compaction * 0.22)
        disorder_order = bounded(0.45 + disorder * 0.35 - closure * 0.16 - segment_compaction * 0.10)
        timepoints.append({
            "timepoint": timepoint,
            "residue_exposure": exposure,
            "segment_compaction": segment_compaction,
            "contact_probability": contact_probability,
            "operator_activation": bounded(sum(strengths.values()) / max(1, len(operator_field["operators"]))),
            "frustration": bounded(strengths["frustration_operator"] + max(0.0, 0.35 - contact_probability)),
            "state_basin_occupancy": basin,
            "interface_readiness": interface_readiness,
            "disorder_order_balance": disorder_order,
            "proteostasis_routing": proteostasis_routing,
            "assembly_required_core": bounded(segment_compaction if mechanism_class == "assembly_required_folding" else 0.0),
            "partner_completed_core": bounded(interface_readiness if mechanism_class == "assembly_required_folding" else 0.0),
            "interface_buried_hydrophobicity": bounded(
                metrics["hydrophobic_density"] * interface_readiness
                if mechanism_class == "assembly_required_folding"
                else 0.0
            ),
            "monomer_incomplete_topology": bounded(
                basin.get("monomer_incomplete_topology", 0.0)
                if mechanism_class == "assembly_required_folding"
                else 0.0
            ),
            "assembly_ambiguous": bounded(
                basin.get("assembly_ambiguous_basin", 0.0)
                if mechanism_class == "assembly_required_folding"
                else 0.0
            ),
            "metal_cluster_geometry": bounded(
                basin.get("metal_cluster_locked_basin", 0.0)
                if mechanism_class == "metal_cluster_and_ligand_locked_basin"
                else 0.0
            ),
            "coordination_shell_integrity": bounded(
                interface_readiness
                if mechanism_class == "metal_cluster_and_ligand_locked_basin"
                else 0.0
            ),
            "ligand_locked_basin": bounded(
                basin.get("ligand_locked_basin", 0.0)
                if mechanism_class == "metal_cluster_and_ligand_locked_basin"
                else 0.0
            ),
            "apo_holo_basin_shift": bounded(
                basin.get("apo_unlocked_basin", 0.0)
                if mechanism_class == "metal_cluster_and_ligand_locked_basin"
                else 0.0
            ),
            "metal_ligand_ambiguous": 0.0,
            "IDR_boundary": bounded(
                basin.get("disorder_boundary_ensemble", 0.0)
                if mechanism_class == "disorder_boundary_and_fold_upon_binding"
                else 0.0
            ),
            "structured_domain_plus_IDR_tail": bounded(
                basin.get("disorder_boundary_ensemble", 0.0) * (1.0 - basin.get("compact_single_fold", 0.0))
                if mechanism_class == "disorder_boundary_and_fold_upon_binding"
                else 0.0
            ),
            "fold_upon_binding_region": bounded(
                basin.get("fold_upon_binding_basin", 0.0)
                if mechanism_class == "disorder_boundary_and_fold_upon_binding"
                else 0.0
            ),
            "phase_prone_low_complexity": bounded(
                basin.get("phase_prone_low_complexity", 0.0)
                if mechanism_class == "disorder_boundary_and_fold_upon_binding"
                else 0.0
            ),
            "flexible_loop_not_disorder": bounded(
                0.0 if mechanism_class == "disorder_boundary_and_fold_upon_binding" else max(0.0, 0.22 - metrics["mean_disorder"])
            ),
            "disorder_with_local_motif": bounded(
                min(basin.get("disorder_boundary_ensemble", 0.0), basin.get("fold_upon_binding_basin", 0.0))
                if mechanism_class == "disorder_boundary_and_fold_upon_binding"
                else 0.0
            ),
            "closed_beta_topology": bounded(
                basin.get("closed_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology"
                else 0.0
            ),
            "strand_register": bounded(
                basin.get("strand_register", 0.0)
                if mechanism_class == "beta_closure_topology"
                else 0.0
            ),
            "beta_sheet_closure": bounded(
                basin.get("closed_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology"
                else 0.0
            ),
            "soluble_beta_barrel": bounded(
                basin.get("closed_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology" and beta_subtype == "soluble_beta_barrel"
                else 0.0
            ),
            "membrane_beta_barrel": bounded(
                basin.get("closed_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology" and beta_subtype == "membrane_beta_barrel"
                else 0.0
            ),
            "beta_propeller_repeat_closure": bounded(
                basin.get("closed_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology" and beta_subtype == "beta_propeller_repeat_closure"
                else 0.0
            ),
            "beta_sandwich_core": bounded(
                basin.get("closed_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology" and beta_subtype == "beta_sandwich_core"
                else 0.0
            ),
            "jelly_roll_wrap": bounded(
                basin.get("closed_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology" and beta_subtype == "jelly_roll_wrap"
                else 0.0
            ),
            "greek_key_beta_lock": bounded(
                basin.get("closed_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology" and beta_subtype == "greek_key_beta_lock"
                else 0.0
            ),
            "beta_helix_solenoid_stack": bounded(
                basin.get("closed_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology" and beta_subtype == "beta_helix_solenoid_stack"
                else 0.0
            ),
            "alpha_beta_barrel_distinction": bounded(
                basin.get("closed_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology" and beta_subtype == "alpha_beta_barrel_distinction"
                else 0.0
            ),
            "open_beta_sheet_ambiguous": bounded(
                basin.get("open_beta_sheet_ambiguous", 0.0)
                if mechanism_class == "beta_closure_topology"
                else 0.0
            ),
            "closed_beta_confident": bounded(
                basin.get("closed_beta_topology", 0.0) - basin.get("wrong_beta_topology", 0.0) * 0.30
                if mechanism_class == "beta_closure_topology"
                else 0.0
            ),
            "closed_beta_ambiguous": bounded(
                max(basin.get("open_beta_sheet_ambiguous", 0.0), basin.get("wrong_beta_topology", 0.0))
                if mechanism_class == "beta_closure_topology"
                else 0.0
            ),
            "strand_register_insufficient": bounded(
                0.36 - basin.get("strand_register", 0.0)
                if mechanism_class == "beta_closure_topology"
                else 0.0
            ),
            "beta_topology_conflict": bounded(
                basin.get("wrong_beta_topology", 0.0)
                if mechanism_class == "beta_closure_topology"
                else 0.0
            ),
            "multidomain_allostery": bounded(
                basin.get("multidomain_allosteric_basin", 0.0)
                if mechanism_class == "multidomain_allosteric_architecture"
                else 0.0
            ),
            "domain_boundary": bounded(
                basin.get("domain_boundary", 0.0)
                if mechanism_class == "multidomain_allosteric_architecture"
                else 0.0
            ),
            "hinge_region": bounded(
                basin.get("hinge_region", 0.0)
                if mechanism_class == "multidomain_allosteric_architecture"
                else 0.0
            ),
            "interdomain_lock": bounded(
                basin.get("interdomain_locked_basin", 0.0)
                if mechanism_class == "multidomain_allosteric_architecture"
                else 0.0
            ),
            "allosteric_basin_shift": bounded(
                basin.get("allosteric_basin_shift", 0.0)
                if mechanism_class == "multidomain_allosteric_architecture"
                else 0.0
            ),
            "domain_reorientation": bounded(
                basin.get("domain_reorientation_basin", 0.0)
                if mechanism_class == "multidomain_allosteric_architecture"
                else 0.0
            ),
            "modular_architecture": bounded(
                basin.get("modular_architecture", 0.0)
                if mechanism_class == "multidomain_allosteric_architecture"
                else 0.0
            ),
            "domain_swapping": bounded(
                basin.get("domain_swapping", 0.0)
                if mechanism_class == "multidomain_allosteric_architecture" and multidomain_subtype == "domain_swapping"
                else 0.0
            ),
        })
    final = timepoints[-1]
    return {
        "kind": "PROTEIN_ESPERANTO_COARSE_OPERATOR_TRAJECTORY_v0",
        "mechanism_class": mechanism_class,
        "perturbation_id": (perturbation or {}).get("perturbation_id", "wild_type_or_reference"),
        "timepoints": timepoints,
        "final_state_summary": final,
        "predicted_contact_interaction_probability_map": contact_probability_map(sequence_field, operator_field, mechanism_class),
        "predicted_state_basin_occupancy": final["state_basin_occupancy"],
        "coordinate_truth_used": False,
        "atomistic_md_executed": False,
    }


def contact_probability_map(sequence_field: dict[str, Any], operator_field: dict[str, Any], mechanism_class: str) -> list[dict[str, Any]]:
    segments = sequence_field["segments"]
    if not segments:
        return []
    pairs: list[dict[str, Any]] = []
    if mechanism_class == "intrinsic_disorder_phase_separation":
        candidates = _strong_segments(sequence_field, "aromatic_density", limit=4)
        for left in candidates:
            for right in candidates:
                if left["segment_id"] >= right["segment_id"]:
                    continue
                pairs.append({
                    "segment_a": left["segment_id"],
                    "segment_b": right["segment_id"],
                    "probability": bounded(0.12 + 0.25 * (left["aromatic_density"] + right["aromatic_density"])),
                    "interaction_type": "weak_multivalent_phase_contact",
                })
    elif mechanism_class == "disorder_boundary_and_fold_upon_binding":
        disorder_top = _strong_segments(sequence_field, "disorder_density", limit=2)
        interface_top = _strong_segments(sequence_field, "interface_density", limit=2)
        for left in disorder_top:
            for right in interface_top:
                if left["segment_id"] == right["segment_id"]:
                    continue
                pairs.append({
                    "segment_a": left["segment_id"],
                    "segment_b": right["segment_id"],
                    "probability": bounded(0.18 + _operator_strength(operator_field, "interface_operator") * 0.28),
                    "interaction_type": "idr_boundary_fold_upon_binding_contact",
                })
    elif mechanism_class == "beta_closure_topology":
        beta_top = _strong_segments(sequence_field, "beta_propensity_density", limit=4)
        if len(beta_top) < 2:
            beta_top = _strong_segments(sequence_field, "hydrophobic_density", limit=4)
        for left in beta_top[:2]:
            for right in beta_top[2:]:
                if left["segment_id"] == right["segment_id"]:
                    continue
                pairs.append({
                    "segment_a": left["segment_id"],
                    "segment_b": right["segment_id"],
                    "probability": bounded(0.20 + _operator_strength(operator_field, "closure_operator") * 0.32),
                    "interaction_type": "beta_strand_register_closure",
                })
    elif mechanism_class == "multidomain_allosteric_architecture":
        interface_top = _strong_segments(sequence_field, "interface_density", limit=4)
        if len(interface_top) < 2:
            interface_top = _strong_segments(sequence_field, "hydrophobic_density", limit=4)
        for left in interface_top[:2]:
            for right in interface_top[2:]:
                if left["segment_id"] == right["segment_id"]:
                    continue
                pairs.append({
                    "segment_a": left["segment_id"],
                    "segment_b": right["segment_id"],
                    "probability": bounded(
                        0.22
                        + _operator_strength(operator_field, "interface_operator") * 0.26
                        + _operator_strength(operator_field, "dual_basin_switch_operator") * 0.12
                    ),
                    "interaction_type": "interdomain_allosteric_lock",
                })
    elif mechanism_class == "short_region_host_interface_hijacking":
        cterm = segments[-1]
        pairs.append({
            "segment_a": cterm["segment_id"],
            "segment_b": "host_RAE1_NUP98_surface",
            "probability": bounded(0.72 + _operator_strength(operator_field, "host_hijack_operator") * 0.18),
            "interaction_type": "host_interface_capture",
        })
    elif mechanism_class == "metamorphic_fold_switching":
        pairs.append({
            "segment_a": segments[0]["segment_id"],
            "segment_b": segments[-1]["segment_id"],
            "probability": bounded(0.45 + _operator_strength(operator_field, "dual_basin_switch_operator") * 0.22),
            "interaction_type": "state_basin_contact_rewrite",
        })
    elif mechanism_class == "membrane_multidomain_folding_proteostasis":
        top_membrane = _strong_segments(sequence_field, "membrane_density", limit=2)
        focus = segments[min(len(segments) - 1, 42)] if len(segments) > 42 else segments[len(segments) // 2]
        for segment in top_membrane:
            pairs.append({
                "segment_a": focus["segment_id"],
                "segment_b": segment["segment_id"],
                "probability": bounded(0.34 + _operator_strength(operator_field, "interface_operator") * 0.28),
                "interaction_type": "interdomain_membrane_interface",
            })
    elif mechanism_class == "cofactor_ligand_assisted_stabilization":
        top = _strong_segments(sequence_field, "interface_density", limit=3)
        for segment in top:
            pairs.append({
                "segment_a": segment["segment_id"],
                "segment_b": "ligand_or_cofactor_pocket",
                "probability": bounded(0.36 + _operator_strength(operator_field, "interface_operator") * 0.34),
                "interaction_type": "cofactor_stabilized_interface",
            })
    elif mechanism_class == "metal_cluster_and_ligand_locked_basin":
        top = _strong_segments(sequence_field, "interface_density", limit=3)
        for segment in top:
            pairs.append({
                "segment_a": segment["segment_id"],
                "segment_b": "metal_cluster_or_ligand_locked_pocket",
                "probability": bounded(0.38 + _operator_strength(operator_field, "interface_operator") * 0.36),
                "interaction_type": "coordination_shell_locked_interface",
            })
    elif mechanism_class == "assembly_required_folding":
        top = _strong_segments(sequence_field, "interface_density", limit=3)
        for segment in top:
            pairs.append({
                "segment_a": segment["segment_id"],
                "segment_b": "partner_completed_core_interface",
                "probability": bounded(0.32 + _operator_strength(operator_field, "interface_operator") * 0.38),
                "interaction_type": "assembly_required_partner_completion",
            })
    elif mechanism_class == "oligomerization_controlled_folding":
        top = _strong_segments(sequence_field, "interface_density", limit=3)
        for segment in top:
            pairs.append({
                "segment_a": segment["segment_id"],
                "segment_b": "partner_copy_interface",
                "probability": bounded(0.34 + _operator_strength(operator_field, "interface_operator") * 0.36),
                "interaction_type": "assembly_stabilized_interface",
            })
    else:
        top = _strong_segments(sequence_field, "hydrophobic_density", limit=4)
        for left in top[:2]:
            for right in top[2:]:
                pairs.append({
                    "segment_a": left["segment_id"],
                    "segment_b": right["segment_id"],
                    "probability": bounded(0.24 + _operator_strength(operator_field, "closure_operator") * 0.42),
                    "interaction_type": "coarse_closure_contact",
                })
    return pairs[:8]


def perturbation_table(
    *,
    sequence_field: dict[str, Any],
    operator_field: dict[str, Any],
    mechanism_class: str,
    perturbations: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    baseline = simulate_operator_trajectory(
        sequence_field=sequence_field,
        operator_field=operator_field,
        mechanism_class=mechanism_class,
    )
    baseline_final = baseline["final_state_summary"]
    rows: list[dict[str, Any]] = []
    for perturbation in perturbations:
        trajectory = simulate_operator_trajectory(
            sequence_field=sequence_field,
            operator_field=operator_field,
            mechanism_class=mechanism_class,
            perturbation=perturbation,
        )
        final = trajectory["final_state_summary"]
        metric = str(perturbation.get("metric", "operator_activation"))
        baseline_value = _extract_metric(baseline_final, metric)
        value = _extract_metric(final, metric)
        expected_direction = str(perturbation.get("expected_direction", "change"))
        observed_direction = _direction(value - baseline_value)
        expected_met = expected_direction == observed_direction or (
            expected_direction == "change" and observed_direction != "unchanged"
        )
        rows.append({
            "perturbation_id": perturbation["perturbation_id"],
            "perturbation": perturbation["description"],
            "metric": metric,
            "baseline_value": baseline_value,
            "perturbed_value": value,
            "expected_direction": expected_direction,
            "observed_direction": observed_direction,
            "direction_passed": expected_met,
            "trajectory_final_state": final,
        })
    return rows


def _extract_metric(final_state: dict[str, Any], metric: str) -> float:
    if metric.startswith("basin:"):
        basin = metric.split(":", 1)[1]
        return float(final_state["state_basin_occupancy"].get(basin, 0.0))
    return float(final_state.get(metric, 0.0))


def _direction(delta: float) -> str:
    if delta >= 0.04:
        return "increase"
    if delta <= -0.04:
        return "decrease"
    return "unchanged"


def validate_against_holdout(
    *,
    sealed_packet: dict[str, Any],
    holdout: dict[str, Any],
) -> dict[str, Any]:
    final = sealed_packet["trajectory_summary"]["final_state_summary"]
    mechanism = sealed_packet["selected_mechanism_grammar"]["mechanism_class"]
    expected = holdout["expected_mechanism_class"]
    checks: list[dict[str, Any]] = []
    checks.append({
        "check_id": "mechanism_class_matches_holdout",
        "passed": mechanism == expected,
        "observed": mechanism,
        "expected": expected,
    })
    for check in holdout.get("expected_observables", []):
        value = _extract_metric(final, check["metric"])
        threshold = float(check["threshold"])
        comparator = check.get("comparator", ">=")
        passed = value >= threshold if comparator == ">=" else value <= threshold
        checks.append({
            "check_id": check["check_id"],
            "passed": passed,
            "metric": check["metric"],
            "observed": value,
            "comparator": comparator,
            "threshold": threshold,
        })
    score_label = "supported" if all(row["passed"] for row in checks) else "contradicted"
    if sealed_packet["evidence_manifest"]["coordinate_derived_source_count_before_prediction"] > 0:
        score_label = "blocked_for_leakage"
    return {
        "kind": "PROTEIN_ESPERANTO_POST_SEAL_VALIDATION_v0",
        "target_id": sealed_packet["target_id"],
        "sealed_prediction_hash": sealed_packet["prediction_hash"],
        "holdout_opened_after_prediction_hash": holdout.get("holdout_opened_after_prediction_hash"),
        "score_label": score_label,
        "checks": checks,
        "coordinate_derived_sources_used_before_prediction": sealed_packet["evidence_manifest"]["coordinate_derived_source_count_before_prediction"],
        "internal_runtime_sources_used_for_prediction": sealed_packet["evidence_manifest"]["internal_runtime_source_count_for_prediction"],
        "postseal_sources": holdout.get("postseal_sources", []),
    }


def make_openmm_bridge_spec() -> dict[str, Any]:
    mappings = [
        {
            "operator": "closure_operator",
            "custom_force_term": "weak attractive bias between predicted coarse regions",
            "guard": "bias may be compared with unbiased baseline; it must not encode holdout coordinates",
        },
        {
            "operator": "disorder_operator",
            "custom_force_term": "penalty against over-collapse for high-disorder regions",
            "guard": "does not force a target fold",
        },
        {
            "operator": "interface_operator",
            "custom_force_term": "partner-proximity or surface-readiness bias",
            "guard": "partner context must come from sealed operator grammar",
        },
        {
            "operator": "dual_basin_switch_operator",
            "custom_force_term": "two-state basin potential with context-dependent weights",
            "guard": "must preserve competing basins rather than average them",
        },
        {
            "operator": "membrane_pressure_operator",
            "custom_force_term": "environment-dependent burial/exposure term",
            "guard": "membrane field cannot replace post-seal validation",
        },
        {
            "operator": "proteostasis_operator",
            "custom_force_term": "routing/state-weight observable outside atomistic force field core",
            "guard": "quality-control readout remains an observable, not a coordinate restraint",
        },
        {
            "operator": "host_hijack_operator",
            "custom_force_term": "short motif host-interface bias",
            "guard": "host-interface coordinates remain holdout-only before sealing",
        },
    ]
    return {
        "kind": "V56_OPERATOR_TO_CUSTOM_FORCE_BRIDGE_SPEC_v0",
        "openmm_execution_required_for_v56": False,
        "full_atom_md_first_proof": False,
        "custom_force_mapping_count": len(mappings),
        "mappings": mappings,
        "pass_condition": "biased physical simulation improves sealed observables over unbiased baseline without encoding the answer",
        "fail_condition": "custom forces force the holdout answer or use coordinate truth before sealing",
    }


def build_sealed_simulation_packet(
    *,
    target_id: str,
    target_name: str,
    sequence: str,
    sources: list[dict[str, Any]],
    focus_regions: list[dict[str, Any]] | None = None,
    perturbations: list[dict[str, Any]] | None = None,
    forced_grammar: str | None = None,
) -> dict[str, Any]:
    sequence_field = build_sequence_field(sequence)
    evidence_manifest = evidence_boundary_gate(sources)
    mechanism = select_mechanism_grammar(
        sequence_field=sequence_field,
        evidence_manifest=evidence_manifest,
        sources=sources,
        forced_grammar=forced_grammar,
    )
    operator_field = build_operator_field(
        sequence_field=sequence_field,
        mechanism=mechanism,
        evidence_manifest=evidence_manifest,
        focus_regions=focus_regions,
    )
    trajectory = simulate_operator_trajectory(
        sequence_field=sequence_field,
        operator_field=operator_field,
        mechanism_class=mechanism["mechanism_class"],
    )
    perturbation_rows = perturbation_table(
        sequence_field=sequence_field,
        operator_field=operator_field,
        mechanism_class=mechanism["mechanism_class"],
        perturbations=perturbations or [],
    )
    packet = {
        "kind": "V52_COARSE_OPERATOR_FOLDING_SIMULATOR_MVP_PACKET_v0",
        "target_id": target_id,
        "target_name": target_name,
        "input_evidence_manifest": {
            "source_count": len(sources),
            "source_ids": [str(source.get("source_id", source.get("accession", "unknown_source"))) for source in sources],
        },
        "evidence_manifest": evidence_manifest,
        "selected_mechanism_grammar": mechanism,
        "operator_field": operator_field,
        "initial_sequence_field_map": sequence_field,
        "trajectory_summary": trajectory,
        "predicted_contact_interaction_probability_map": trajectory["predicted_contact_interaction_probability_map"],
        "predicted_state_basin_occupancy": trajectory["predicted_state_basin_occupancy"],
        "predicted_perturbation_table": perturbation_rows,
        "predicted_falsifiers": [
            row["falsified_by"] for row in operator_field["operators"]
        ] or ["insufficient evidence should remain abstained until allowed sources exist"],
        "sealed_before_holdout": True,
        "coordinate_truth_used_before_prediction": False,
        "atomistic_md_executed": False,
        "folding_problem_solved": False,
    }
    packet["prediction_hash"] = stable_hash({key: value for key, value in packet.items() if key != "prediction_hash"})
    return packet


def shuffled_sequence(sequence: str) -> str:
    sequence = normalize_sequence(sequence)
    buckets = [sequence[index::3] for index in range(3)]
    return "".join(bucket[::-1] for bucket in buckets)


def deterministic_random_sequence(length: int) -> str:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(alphabet[(index * 7 + 3) % len(alphabet)] for index in range(length))


def sequence_operator_coherence(packet: dict[str, Any]) -> float:
    mechanism = packet["selected_mechanism_grammar"]["mechanism_class"]
    field = packet["initial_sequence_field_map"]["global_metrics"]
    segments = packet["initial_sequence_field_map"]["segments"]
    operators = packet["operator_field"]["operators"]
    activation = _avg(row["activation_strength"] for row in operators)
    if mechanism == "intrinsic_disorder_phase_separation":
        local = max(
            (
                row["low_complexity_density"]
                + row["aromatic_density"]
                + row["disorder_density"]
                for row in segments
            ),
            default=0.0,
        )
        support = 0.35 * (field["low_complexity_density"] + field["aromatic_density"] + field["mean_disorder"]) + 0.65 * local
    elif mechanism == "disorder_boundary_and_fold_upon_binding":
        local_disorder = max((row["low_complexity_density"] + row["disorder_density"] for row in segments), default=0.0)
        local_interface = max((row["interface_density"] for row in segments), default=0.0)
        support = (
            0.30 * (field["low_complexity_density"] + field["mean_disorder"])
            + 0.35 * local_disorder
            + 0.20 * local_interface
            + 0.15 * activation
        )
    elif mechanism == "beta_closure_topology":
        local_beta = max((row["beta_propensity_density"] for row in segments), default=0.0)
        local_interface = max((row["interface_density"] for row in segments), default=0.0)
        local_aromatic = max((row["aromatic_density"] for row in segments), default=0.0)
        support = (
            0.24 * field["hydrophobic_density"]
            + 0.26 * field["beta_propensity_density"]
            + 0.18 * field["mean_interface"]
            + 0.20 * local_beta
            + 0.12 * max(local_interface, local_aromatic)
            + activation
        )
    elif mechanism == "multidomain_allosteric_architecture":
        local_interface = max((row["interface_density"] for row in segments), default=0.0)
        local_hydrophobic = max((row["hydrophobic_density"] for row in segments), default=0.0)
        local_aromatic = max((row["aromatic_density"] for row in segments), default=0.0)
        support = (
            0.20 * field["hydrophobic_density"]
            + 0.26 * field["mean_interface"]
            + 0.26 * local_interface
            + 0.12 * max(local_hydrophobic, local_aromatic)
            + activation
        )
    elif mechanism == "membrane_multidomain_folding_proteostasis":
        local = max((row["membrane_density"] for row in segments), default=0.0)
        support = 0.50 * field["mean_membrane"] + 0.50 * local + activation
    elif mechanism == "metamorphic_fold_switching":
        support = activation
    elif mechanism == "short_region_host_interface_hijacking":
        cterminal = segments[-2:] if len(segments) >= 2 else segments
        local = _avg(row["interface_density"] for row in cterminal)
        support = 0.35 * field["mean_interface"] + 0.65 * local + activation
    elif mechanism == "cofactor_ligand_assisted_stabilization":
        local = max((row["interface_density"] for row in segments), default=0.0)
        support = 0.30 * field["hydrophobic_density"] + 0.35 * field["mean_interface"] + 0.35 * local + activation
    elif mechanism == "metal_cluster_and_ligand_locked_basin":
        local = max((row["interface_density"] + row["aromatic_density"] for row in segments), default=0.0)
        support = 0.24 * field["hydrophobic_density"] + 0.36 * field["mean_interface"] + 0.40 * local + activation
    elif mechanism == "assembly_required_folding":
        local = max((row["interface_density"] for row in segments), default=0.0)
        membrane_like = max((row["membrane_density"] for row in segments), default=0.0)
        support = 0.25 * field["hydrophobic_density"] + 0.25 * field["mean_interface"] + 0.30 * local + 0.20 * membrane_like + activation
    elif mechanism == "oligomerization_controlled_folding":
        local = max((row["interface_density"] for row in segments), default=0.0)
        support = 0.30 * field["hydrophobic_density"] + 0.25 * field["mean_interface"] + 0.45 * local + activation
    elif mechanism == "globular_closure":
        local = max((row["hydrophobic_density"] for row in segments), default=0.0)
        support = 0.50 * field["hydrophobic_density"] + 0.50 * local + activation
    else:
        support = 0.0
    return bounded(0.5 * activation + 0.5 * support)
