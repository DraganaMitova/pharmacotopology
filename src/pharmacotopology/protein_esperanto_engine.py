from __future__ import annotations

"""Coarse Protein Esperanto operator grammar and simulator.

This module is intentionally not an atomistic simulator.  It turns allowed
sequence/evidence into a coarse state field, selects a mechanism grammar,
builds operator activations, and evolves mechanism-level observables that can
be sealed before holdout validation.
"""

from collections import Counter
from copy import deepcopy
from functools import lru_cache
from hashlib import sha256
import importlib.util
from itertools import permutations
import json
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
    "secretory_disulfide_redox_topology",
    "signal_peptide_vs_true_tm_routing",
    "coiled_coil_register_topology",
    "repeat_solenoid_topology",
    "knotted_topology",
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
    "secretory_redox_topology",
    "signal_peptide_tm_boundary",
    "coiled_coil_register",
    "repeat_solenoid",
    "threaded_knot_topology",
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
    "disulfide_pairing_operator",
    "secretory_redox_operator",
    "signal_peptide_routing_operator",
    "tm_insertion_operator",
    "cleavage_context_operator",
    "secretory_routing_operator",
    "heptad_register_operator",
    "coiled_coil_interface_operator",
    "oligomeric_register_operator",
    "register_shift_frustration_operator",
    "repeat_phase_operator",
    "solenoid_axis_operator",
    "local_repeat_closure_operator",
    "global_repeat_stack_operator",
    "repeat_boundary_frustration_operator",
    "threading_operator",
    "topological_closure_operator",
    "long_range_threading_operator",
    "slipknot_intermediate_operator",
    "knotting_frustration_operator",
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
    "secretory_redox_context",
    "disulfide_pairing_topology",
    "cysteine_pairing_constraint",
    "extracellular_stabilized_fold",
    "glycosylation_context",
    "redox_mispaired_frustration",
    "signal_peptide_removed_context",
    "secretory_quality_control",
    "signal_peptide_routing_context",
    "cleavage_site_context",
    "n_terminal_secretory_hydrophobic_patch",
    "true_transmembrane_span_context",
    "single_pass_tm_conflict",
    "multi_pass_tm_conflict",
    "secretory_lumenal_routing",
    "membrane_insertion_routing",
    "signal_anchor_ambiguity",
    "heptad_register_context",
    "hydrophobic_repeat_phase",
    "parallel_antiparallel_register",
    "oligomeric_coiled_coil_core",
    "register_shift_frustration",
    "coiled_coil_assembly_dependency",
    "leucine_zipper_context",
    "repeat_unit_context",
    "solenoid_axis_context",
    "curved_repeat_stack",
    "local_repeat_closure",
    "global_repeat_topology",
    "repeat_phase_alignment",
    "repeat_boundary_frustration",
    "ankyrin_armadillo_tpr_lrr_context",
    "knot_core_context",
    "threading_loop_context",
    "slipknot_intermediate_context",
    "topological_closure_constraint",
    "long_range_threading_dependency",
    "knotting_frustration",
    "unknotted_decoy_dominance",
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

SECRETORY_DISULFIDE_CONTEXT_TOKENS = [
    "disulfide_secretory_redox_context",
    "disulfide secretory redox context",
    "disulfide_bond_topology",
    "disulfide bond topology",
    "disulphide bond topology",
    "secretory_redox_context",
    "secretory redox context",
    "cysteine_pairing_constraint",
    "cysteine pairing constraint",
    "extracellular_stabilized_fold",
    "extracellular stabilized fold",
    "glycosylation_context",
    "glycosylation context",
    "signal_peptide_removed_context",
    "signal peptide removed context",
    "secretory_quality_control",
    "secretory quality control",
    "disulfide",
    "disulphide",
    "secreted",
    "secretory",
    "extracellular",
    "cysteine-rich",
    "cysteine rich",
]

SIGNAL_PEPTIDE_ONLY_CONTEXT_TOKENS = [
    "signal peptide only",
    "signal-peptide only",
    "cleaved signal peptide",
    "signal sequence only",
    "signal peptide without disulfide",
    "signal sequence without cysteine pairing",
]

SIGNAL_PEPTIDE_ROUTING_CONTEXT_TOKENS = [
    "signal_peptide_vs_true_TM",
    "signal_peptide_vs_true_tm",
    "signal peptide vs true tm",
    "signal peptide",
    "signal-peptide",
    "signal sequence",
    "cleaved signal peptide",
    "cleavable signal peptide",
    "secretory signal peptide",
    "n-terminal signal peptide",
    "n terminal signal peptide",
    "signal_peptide_routing_context",
    "signal peptide routing context",
    "cleavage_site_context",
    "cleavage site",
    "signal peptidase",
    "secretory_lumenal_routing",
    "secretory lumenal routing",
    "lumenal routing",
    "er lumen",
    "export signal",
    "sec signal",
]

TRUE_TM_ROUTING_CONTEXT_TOKENS = [
    "true_transmembrane_span_context",
    "true transmembrane span context",
    "true transmembrane",
    "transmembrane span",
    "tm span",
    "transmembrane helix",
    "tm helix",
    "multi-pass transmembrane",
    "multipass transmembrane",
    "multi pass transmembrane",
    "integral membrane",
    "membrane_insertion_routing",
    "membrane insertion routing",
]

SIGNAL_ANCHOR_AMBIGUITY_TOKENS = [
    "signal_anchor_ambiguity",
    "signal anchor ambiguity",
    "signal anchor",
    "signal-anchor",
    "uncleaved signal anchor",
    "uncleaved signal-anchor",
]

CYS_HIS_METAL_COORDINATION_TOKENS = [
    "cys-his coordination",
    "cys his coordination",
    "cys2his2",
    "c2h2",
    "zinc finger",
    "zinc-finger",
    "zinc-binding",
    "zinc binding",
    "iron-sulfur cysteine",
    "cysteine ligated metal",
]

SECRETORY_DISULFIDE_STATE_VARIABLES = [
    "secretory_redox_context",
    "disulfide_pairing_topology",
    "cysteine_pairing_constraint",
    "extracellular_stabilized_fold",
    "glycosylation_context",
    "redox_mispaired_frustration",
    "signal_peptide_removed_context",
    "secretory_quality_control",
]

SIGNAL_PEPTIDE_ROUTING_STATE_VARIABLES = [
    "signal_peptide_routing_context",
    "cleavage_site_context",
    "n_terminal_secretory_hydrophobic_patch",
    "true_transmembrane_span_context",
    "single_pass_tm_conflict",
    "multi_pass_tm_conflict",
    "secretory_lumenal_routing",
    "membrane_insertion_routing",
    "signal_anchor_ambiguity",
]

COILED_COIL_REGISTER_CONTEXT_TOKENS = [
    "coiled_coil_register",
    "coiled-coil register",
    "coiled coil register",
    "coiled-coil",
    "coiled coil",
    "heptad_repeat",
    "heptad repeat",
    "register_alignment",
    "register alignment",
    "leucine zipper",
    "leucine_zipper_context",
    "parallel_vs_antiparallel_coil",
    "parallel antiparallel coil",
    "oligomeric_coiled_coil_core",
    "oligomeric coiled coil core",
    "heptad_register_context",
    "hydrophobic_repeat_phase",
]

REPEAT_SOLENOID_CONTEXT_TOKENS = [
    "repeat_solenoid_topology",
    "repeat solenoid topology",
    "repeat_unit",
    "repeat unit",
    "solenoid_axis",
    "solenoid axis",
    "curved_repeat_stack",
    "curved repeat stack",
    "local_repeat_closure",
    "local repeat closure",
    "global_repeat_topology",
    "global repeat topology",
    "repeat_phase_alignment",
    "ankyrin repeat",
    "ankyrin_armadillo_tpr_lrr_context",
    "armadillo repeat",
    "tpr repeat",
    "leucine-rich repeat",
    "leucine rich repeat",
    "pentapeptide repeat",
]

KNOTTED_TOPOLOGY_CONTEXT_TOKENS = [
    "knotted_topology",
    "knotted topology",
    "knot_core_context",
    "knot core",
    "threading_loop_context",
    "threading loop",
    "slipknot_intermediate_context",
    "slipknot",
    "slip-knot",
    "topological_closure_constraint",
    "topological closure constraint",
    "long_range_threading_dependency",
    "long-range threading",
    "knotting_frustration",
]

COILED_COIL_REGISTER_STATE_VARIABLES = [
    "heptad_register_context",
    "hydrophobic_repeat_phase",
    "parallel_antiparallel_register",
    "oligomeric_coiled_coil_core",
    "register_shift_frustration",
    "coiled_coil_assembly_dependency",
    "leucine_zipper_context",
]

REPEAT_SOLENOID_STATE_VARIABLES = [
    "repeat_unit_context",
    "solenoid_axis_context",
    "curved_repeat_stack",
    "local_repeat_closure",
    "global_repeat_topology",
    "repeat_phase_alignment",
    "repeat_boundary_frustration",
    "ankyrin_armadillo_tpr_lrr_context",
]

KNOTTED_TOPOLOGY_STATE_VARIABLES = [
    "knot_core_context",
    "threading_loop_context",
    "slipknot_intermediate_context",
    "topological_closure_constraint",
    "long_range_threading_dependency",
    "knotting_frustration",
    "unknotted_decoy_dominance",
]

SELF_DECISION_CANDIDATE_GRAMMARS: dict[str, dict[str, Any]] = {}

E73_WORD_LIFECYCLE = [
    "unseen_pattern",
    "pressure_cluster",
    "candidate_word",
    "proto_grammar",
    "learned_grammar",
    "rejected_or_merged_word",
]

E73_PRESSURE_CHANNELS = [
    "repeated_clean_abstention_pressure",
    "wrong_grammar_pressure",
    "contradiction_pressure",
    "sentinel_pressure",
    "metadata_masking_pressure",
    "perturbation_pressure",
    "physical_execution_mismatch_pressure",
    "compression_pressure",
]

NEGATIVE_EVIDENCE_PRESSURE_CHANNELS = {
    "not_globular_pressure": ["not globular", "not soluble monomer", "no compact globular fold"],
    "not_membrane_pressure": ["not membrane", "not transmembrane", "no transmembrane", "no membrane topology"],
    "not_secretory_pressure": ["not secretory", "not extracellular", "no signal peptide", "not disulfide"],
    "not_metal_pressure": ["not metal", "no metal", "not ligand", "no cofactor", "not cys-his"],
    "not_repeat_pressure": ["not repeat", "not solenoid", "no repeat topology"],
    "not_assembly_pressure": ["not assembly", "not oligomer", "not obligate assembly"],
    "not_beta_pressure": ["not beta", "no beta barrel", "not beta propeller", "not beta topology"],
    "not_knotted_pressure": ["not knot", "not knotted", "unknotted", "no threading topology"],
}

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

SELF_DECISION_LEARNED_GRAMMAR_FAMILIES = {
    "globular_closure": [
        "soluble_monomeric_core_context",
        "complete soluble monomer",
        "standalone soluble fold",
        "single-domain",
        "single domain",
        "folding kinetics",
        "protein folding target",
    ],
    "intrinsic_disorder_phase_separation": [
        "low complexity",
        "phase separation",
        "prion",
        "lcd",
        "fus",
        "tdp-43",
        "tdp43",
    ],
    "disorder_boundary_and_fold_upon_binding": DISORDER_BOUNDARY_CONTEXT_TOKENS,
    "beta_closure_topology": BETA_CLOSURE_CONTEXT_TOKENS,
    "multidomain_allosteric_architecture": MULTIDOMAIN_ALLOSTERIC_CONTEXT_TOKENS,
    "secretory_disulfide_redox_topology": SECRETORY_DISULFIDE_CONTEXT_TOKENS,
    "signal_peptide_vs_true_tm_routing": SIGNAL_PEPTIDE_ROUTING_CONTEXT_TOKENS
    + TRUE_TM_ROUTING_CONTEXT_TOKENS
    + SIGNAL_ANCHOR_AMBIGUITY_TOKENS,
    "coiled_coil_register_topology": COILED_COIL_REGISTER_CONTEXT_TOKENS,
    "repeat_solenoid_topology": REPEAT_SOLENOID_CONTEXT_TOKENS,
    "knotted_topology": KNOTTED_TOPOLOGY_CONTEXT_TOKENS,
    "membrane_multidomain_folding_proteostasis": STRONG_MEMBRANE_CONTEXT_TOKENS
    + ["cftr", "f508", "proteostasis", "trafficking", "transmembrane", "membrane"],
    "metamorphic_fold_switching": [
        "fold switch",
        "fold-switch",
        "metamorphic",
        "dual basin",
        "dual-basin",
        "competing basin",
    ],
    "short_region_host_interface_hijacking": ["host hijack", "rae1", "nup98", "nucleocytoplasmic", "interferon"],
    "fold_upon_binding_disorder": ["fold_upon_binding", "fold upon binding", "partner motif", "binding motif"],
    "cofactor_ligand_assisted_stabilization": COFACTOR_CONTEXT_TOKENS + ["cofactor", "ligand"],
    "metal_cluster_and_ligand_locked_basin": METAL_CLUSTER_CONTEXT_TOKENS + LIGAND_LOCKED_CONTEXT_TOKENS,
    "assembly_required_folding": ASSEMBLY_REQUIRED_CONTEXT_TOKENS,
    "oligomerization_controlled_folding": OLIGOMER_CONTEXT_TOKENS + ["oligomer", "multimer"],
}

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
    "secretory_disulfide_redox_topology": {
        "marks": [
            "secretory_redox_context",
            "disulfide_pairing_topology",
            "cysteine_pairing_constraint",
            "extracellular_stabilized_fold",
            "glycosylation_context",
            "signal_peptide_removed_context",
        ],
        "pressures": ["secretory_redox_environment", "extracellular_quality_control", "glycosylation_context"],
        "operators": [
            "disulfide_pairing_operator",
            "secretory_redox_operator",
            "frustration_operator",
            "proteostasis_operator",
            "closure_operator",
        ],
        "state_change": "generic cysteine-rich or secretory text to paired disulfide/redox topology with extracellular stabilization",
        "testable_effect": "cysteine-pairing or redox perturbation weakens disulfide topology and raises mispaired frustration",
        "null_control": "signal peptide text, odd cysteine noise, true TM topology, or Cys-His metal coordination cannot validate disulfide grammar",
        "falsification_rule": "matched controls beat the real target, or metal/TM/beta/repeat grammar explains the evidence better",
    },
    "signal_peptide_vs_true_tm_routing": {
        "marks": [
            "signal_peptide_routing_context",
            "cleavage_site_context",
            "n_terminal_secretory_hydrophobic_patch",
            "true_transmembrane_span_context",
            "secretory_lumenal_routing",
            "membrane_insertion_routing",
            "signal_anchor_ambiguity",
        ],
        "pressures": ["secretory_translocon", "membrane_insertion", "cleavage_context"],
        "operators": [
            "signal_peptide_routing_operator",
            "tm_insertion_operator",
            "cleavage_context_operator",
            "secretory_routing_operator",
            "membrane_pressure_operator",
            "frustration_operator",
        ],
        "state_change": "generic N-terminal hydrophobic signal to cleavable secretory route while separating true TM and signal-anchor ambiguity",
        "testable_effect": "N-terminal masking or cleavage-context masking weakens secretory routing more than true-TM decoys",
        "null_control": "multi-pass TM, signal-anchor ambiguity, secretory disulfide, coiled/repeat/knotted candidates, and membrane sentinels cannot validate signal-peptide grammar",
        "falsification_rule": "matched controls beat the real signal, true-TM insertion dominates, or wrong grammar explains the boundary better",
    },
    "coiled_coil_register_topology": {
        "marks": [
            "heptad_register_context",
            "hydrophobic_repeat_phase",
            "parallel_antiparallel_register",
            "oligomeric_coiled_coil_core",
            "register_shift_frustration",
            "coiled_coil_assembly_dependency",
            "leucine_zipper_context",
        ],
        "pressures": ["register_alignment", "oligomeric_interface", "hydrophobic_repeat_phase"],
        "operators": [
            "heptad_register_operator",
            "coiled_coil_interface_operator",
            "oligomeric_register_operator",
            "register_shift_frustration_operator",
            "closure_operator",
        ],
        "state_change": "generic helix, assembly, or globular closure to explicit heptad/register-controlled coiled-coil topology",
        "testable_effect": "heptad shuffle or hydrophobic phase shift damages register state more than matched enemy grammars",
        "null_control": "generic oligomer, ordinary helix bundle, or globular hydrophobic core cannot validate coiled-coil register",
        "falsification_rule": "assembly/globular grammar preserves the evidence better, matched heptad controls beat the real target, or register perturbation is inert",
    },
    "repeat_solenoid_topology": {
        "marks": [
            "repeat_unit_context",
            "solenoid_axis_context",
            "curved_repeat_stack",
            "local_repeat_closure",
            "global_repeat_topology",
            "repeat_phase_alignment",
            "repeat_boundary_frustration",
            "ankyrin_armadillo_tpr_lrr_context",
        ],
        "pressures": ["repeat_phase", "solenoid_axis", "local_to_global_repeat_stack"],
        "operators": [
            "repeat_phase_operator",
            "solenoid_axis_operator",
            "local_repeat_closure_operator",
            "global_repeat_stack_operator",
            "repeat_boundary_frustration_operator",
        ],
        "state_change": "generic beta/multidomain/globular readout to repeat-unit phase and solenoid-axis topology",
        "testable_effect": "repeat-order shuffle or boundary masking damages phase-aligned solenoid state",
        "null_control": "long multidomain text, beta closure, or generic repeat words without phase/axis support cannot validate solenoid grammar",
        "falsification_rule": "beta/multidomain/globular grammar explains the target better, repeat controls beat the real target, or repeat perturbation is inert",
    },
    "knotted_topology": {
        "marks": [
            "knot_core_context",
            "threading_loop_context",
            "slipknot_intermediate_context",
            "topological_closure_constraint",
            "long_range_threading_dependency",
            "knotting_frustration",
            "unknotted_decoy_dominance",
        ],
        "pressures": ["threading", "topological_closure", "long_range_order"],
        "operators": [
            "threading_operator",
            "topological_closure_operator",
            "long_range_threading_operator",
            "slipknot_intermediate_operator",
            "knotting_frustration_operator",
        ],
        "state_change": "generic long-range closure to explicit knot/slipknot threading topology",
        "testable_effect": "threading or long-range-order perturbation damages knot state while unknotted decoys fail",
        "null_control": "sequence alone cannot validate knotted topology; explicit non-coordinate topology context is required",
        "falsification_rule": "globular/beta/repeat grammar explains the target better, topology-masked controls beat the real target, or threading perturbation is inert",
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


def protein_esperanto_epistemological_status() -> dict[str, Any]:
    return {
        "language_layer": "protein_esperanto_mechanism_language",
        "is_physical_simulation": False,
        "is_atomistic_md": False,
        "trajectory_interpretation": "operator_state_propagation_not_time_physical_dynamics",
        "contact_map_interpretation": "hypothesized_interaction_language_map_not_native_contact_probability",
        "physical_basis_claim_allowed": False,
        "folding_problem_solved": False,
        "requires_independent_physical_validation": True,
    }


def _active_pressure_names(channels: dict[str, Any]) -> list[str]:
    return sorted(name for name, value in channels.items() if bool(value))


def _negative_evidence_pressure(text: str, final_state: dict[str, Any]) -> dict[str, dict[str, Any]]:
    channels: dict[str, dict[str, Any]] = {}
    state_keys = {
        "not_globular_pressure": ["contact_probability", "segment_compaction"],
        "not_membrane_pressure": ["proteostasis_routing", "membrane_insertion_routing", "true_transmembrane_span_context"],
        "not_secretory_pressure": ["secretory_redox_context", "signal_peptide_routing_context", "secretory_lumenal_routing"],
        "not_metal_pressure": ["metal_cluster_geometry", "ligand_locked_basin", "coordination_shell_integrity"],
        "not_repeat_pressure": ["repeat_unit_context", "global_repeat_topology", "repeat_phase_alignment"],
        "not_assembly_pressure": ["assembly_required_core", "partner_completed_core"],
        "not_beta_pressure": ["closed_beta_topology", "strand_register", "beta_sheet_closure"],
        "not_knotted_pressure": ["knot_core_context", "threading_loop_context", "topological_closure_constraint"],
    }
    for channel, tokens in NEGATIVE_EVIDENCE_PRESSURE_CHANNELS.items():
        explicit_hits = _token_hits(text, tokens)
        opposed_readouts = [
            key
            for key in state_keys.get(channel, [])
            if float(final_state.get(key, 0.0)) > 0.0
        ]
        channels[channel] = {
            "explicit_negative_hits": explicit_hits,
            "opposed_positive_readouts": opposed_readouts,
            "pressure_active": bool(explicit_hits or opposed_readouts),
            "channel_interpretation": "inhibitory_language_pressure_not_a_threshold",
        }
    return channels


def language_acquisition_observation_from_packet(
    packet: dict[str, Any],
    *,
    visible_context_text: str = "",
    matched_control_dominance_passed: bool | None = None,
    sentinel_regression: bool = False,
    physical_execution_mismatch: bool = False,
) -> dict[str, Any]:
    judge = packet["self_decision_judge"]
    final_state = packet["trajectory_summary"]["final_state_summary"]
    selected = packet["selected_mechanism_grammar"]["mechanism_class"]
    final_decision = judge["final_self_decision"]
    pressure_channels = {
        "repeated_clean_abstention_pressure": final_decision.startswith("clean_abstain"),
        "wrong_grammar_pressure": judge["wrong_grammar_separation"] == "wrong_grammar_competes",
        "contradiction_pressure": bool(judge["contradictions"]),
        "sentinel_pressure": bool(sentinel_regression),
        "metadata_masking_pressure": judge["masking_stability"] == "unstable_under_nondefining_mask",
        "perturbation_pressure": judge["operator_basis_stability"] == "nondefining_operator_basis_sensitive",
        "physical_execution_mismatch_pressure": bool(physical_execution_mismatch),
        "compression_pressure": (
            final_decision.startswith("clean_abstain")
            or judge["wrong_grammar_separation"] == "wrong_grammar_competes"
            or bool(judge["contradictions"])
            or bool(physical_execution_mismatch)
        ),
    }
    negative_pressure = _negative_evidence_pressure(visible_context_text.lower(), final_state)
    active_negative = [
        name
        for name, row in negative_pressure.items()
        if row["pressure_active"]
    ]
    active_pressure = _active_pressure_names(pressure_channels)
    runner_ups = list(judge.get("runner_up_mechanisms") or [])
    fingerprint = stable_hash({
        "selected": selected,
        "reason": judge["self_decision_reason"],
        "active_pressure": active_pressure,
        "active_negative": active_negative,
        "runner_ups": runner_ups,
        "cross_view_binding": judge["cross_view_binding"],
        "operator_basis_stability": judge["operator_basis_stability"],
        "temporal_binding": judge["temporal_binding"],
    })[:12]
    return {
        "kind": "E73_LANGUAGE_ACQUISITION_OBSERVATION_v0",
        "target_id": packet["target_id"],
        "selected_mechanism": selected,
        "natural_mechanism": packet["selected_mechanism_grammar"]["natural_mechanism_class"],
        "final_self_decision": final_decision,
        "acceptance_decision": judge["acceptance_decision"],
        "self_decision_reason": judge["self_decision_reason"],
        "pressure_channels": pressure_channels,
        "active_pressure_channels": active_pressure,
        "negative_evidence_pressure": negative_pressure,
        "active_negative_evidence_channels": active_negative,
        "runner_up_mechanisms": runner_ups,
        "matched_control_dominance_passed": matched_control_dominance_passed,
        "wrong_grammar_separation": judge["wrong_grammar_separation"],
        "operator_basis_stability": judge["operator_basis_stability"],
        "temporal_binding": judge["temporal_binding"],
        "coordinate_truth_used_before_prediction": packet["coordinate_truth_used_before_prediction"],
        "physical_basis_claim_allowed": judge["physical_basis_claim_allowed"],
        "folding_problem_solved": packet["folding_problem_solved"],
        "pressure_fingerprint": fingerprint,
    }


def _candidate_word_name(fingerprint: str, observations: list[dict[str, Any]]) -> str:
    pressure_names = Counter(
        pressure
        for observation in observations
        for pressure in observation["active_pressure_channels"]
    )
    negative_names = Counter(
        pressure
        for observation in observations
        for pressure in observation["active_negative_evidence_channels"]
    )
    if pressure_names:
        stem = pressure_names.most_common(1)[0][0].replace("_pressure", "")
    elif negative_names:
        stem = negative_names.most_common(1)[0][0].replace("_pressure", "")
    else:
        stem = "unseen_pattern"
    return f"{stem}_{fingerprint}"


def build_proto_grammar_from_observations(
    *,
    candidate_word: str,
    observations: list[dict[str, Any]],
) -> dict[str, Any]:
    pressure_counter = Counter(
        pressure
        for observation in observations
        for pressure in observation["active_pressure_channels"]
    )
    negative_counter = Counter(
        pressure
        for observation in observations
        for pressure in observation["active_negative_evidence_channels"]
    )
    enemy_counter = Counter(
        enemy
        for observation in observations
        for enemy in observation["runner_up_mechanisms"]
        if enemy != observation["selected_mechanism"]
    )
    selected_counter = Counter(observation["selected_mechanism"] for observation in observations)
    matched_controls_pass = all(
        observation["matched_control_dominance_passed"] is not False
        for observation in observations
    )
    wrong_grammar_challenge_fails = all(
        observation["wrong_grammar_separation"] != "wrong_grammar_competes"
        for observation in observations
    )
    perturbation_paired = all(
        observation["operator_basis_stability"] != "nondefining_operator_basis_sensitive"
        for observation in observations
    )
    sentinels_do_not_regress = not any(
        observation["pressure_channels"]["sentinel_pressure"]
        for observation in observations
    )
    truth_boundary_sealed = not any(
        observation["coordinate_truth_used_before_prediction"]
        for observation in observations
    )
    physical_claim_blocked = not any(
        observation["physical_basis_claim_allowed"] or observation["folding_problem_solved"]
        for observation in observations
    )
    promotion_tests = {
        "selected_grammar_top_internal_readout": all(
            observation["acceptance_decision"] == "accepted"
            for observation in observations
        ),
        "required_views_derived_from_grammar": True,
        "matched_control_dominance_passes": matched_controls_pass,
        "wrong_grammar_challenge_fails": wrong_grammar_challenge_fails,
        "perturbation_response_is_paired_baseline": perturbation_paired,
        "sentinels_do_not_regress": sentinels_do_not_regress,
        "coordinate_native_truth_stays_sealed": truth_boundary_sealed,
        "physical_basis_claim_remains_blocked": physical_claim_blocked,
    }
    promoted = all(promotion_tests.values())
    if promoted:
        lifecycle_state = "learned_grammar"
        promotion_outcome = "promoted_through_proto_grammar_tests"
    elif len(observations) == 1:
        lifecycle_state = "pressure_cluster"
        promotion_outcome = "cleanly_abstained_single_pressure_observation"
    else:
        lifecycle_state = "proto_grammar"
        promotion_outcome = "cleanly_abstained_pressure_support_not_yet_learned"
    return {
        "kind": "E73_PROTO_GRAMMAR_v0",
        "candidate_word": candidate_word,
        "lifecycle_state": lifecycle_state,
        "observation_count": len(observations),
        "pressure_channels": dict(pressure_counter),
        "negative_evidence_channels": dict(negative_counter),
        "proposed_mechanism_class": f"{candidate_word}_mechanism",
        "state_variables": sorted(set(pressure_counter) | set(negative_counter)),
        "operators": [f"{pressure}_operator" for pressure in sorted(pressure_counter)],
        "enemy_grammars": [name for name, _count in enemy_counter.most_common()],
        "definition_by_known_words": {
            "existing_selected_grammars": dict(selected_counter),
            "enemy_grammars": dict(enemy_counter),
            "pressure_components": dict(pressure_counter),
            "negative_pressure_components": dict(negative_counter),
        },
        "usage_by_context": {
            "pressure_fingerprints": sorted({observation["pressure_fingerprint"] for observation in observations}),
            "observation_target_ids": [observation["target_id"] for observation in observations],
            "context_route": "repeated_pressure_support_observation_without_forced_label",
        },
        "matched_controls": "required_before_promotion",
        "perturbation_tests": "paired_baseline_required_before_promotion",
        "abstention_conditions": [
            "single_observation_without_replay",
            "wrong_grammar_competes",
            "sentinel_regression",
            "native_truth_leakage",
            "physical_claim_requested_without_independent_holdout",
        ],
        "physical_execution_expectations": [
            "selected_proto_grammar_must_beat_wrong_grammar",
            "selected_proto_grammar_must_beat_masked_grammar",
            "physical_basis_claim_remains_separate_until_independent_holdout",
        ],
        "promotion_tests": promotion_tests,
        "promotion_outcome": promotion_outcome,
        "merge_rule": "merge_into_existing_word_when_existing_grammar_compresses_pressure_without_new_state_variables",
        "retire_rule": "retire_when_replay_causes_sentinel_regression_or_enemy_grammar_stealing",
    }


def protein_language_acquisition_cortex(observations: list[dict[str, Any]]) -> dict[str, Any]:
    pressure_observations = [
        observation
        for observation in observations
        if observation["active_pressure_channels"] or observation["active_negative_evidence_channels"]
    ]
    grouped: dict[str, list[dict[str, Any]]] = {}
    for observation in pressure_observations:
        grouped.setdefault(observation["pressure_fingerprint"], []).append(observation)
    proto_grammars = [
        build_proto_grammar_from_observations(
            candidate_word=_candidate_word_name(fingerprint, rows),
            observations=rows,
        )
        for fingerprint, rows in sorted(grouped.items())
    ]
    ranked = sorted(
        proto_grammars,
        key=lambda row: (
            row["observation_count"],
            len(row["pressure_channels"]),
            len(row["negative_evidence_channels"]),
            row["candidate_word"],
        ),
        reverse=True,
    )
    learned = [row for row in ranked if row["lifecycle_state"] == "learned_grammar"]
    abstained = [
        row
        for row in ranked
        if row["promotion_outcome"].startswith("cleanly_abstained")
    ]
    return {
        "kind": "E73_PROTEIN_LANGUAGE_ACQUISITION_CORTEX_v0",
        "engine_revision": "E73",
        "word_lifecycle": E73_WORD_LIFECYCLE,
        "pressure_channel_names": E73_PRESSURE_CHANNELS,
        "negative_evidence_pressure_channel_names": sorted(NEGATIVE_EVIDENCE_PRESSURE_CHANNELS),
        "observation_count": len(observations),
        "pressure_observation_count": len(pressure_observations),
        "candidate_words": ranked,
        "candidate_word_count": len(ranked),
        "learned_grammar_promotions": len(learned),
        "cleanly_abstained_candidate_words": len(abstained),
        "candidate_words_ranked_by_endogenous_pressure_support": True,
        "candidate_word_proposal_hash": stable_hash([
            {
                "candidate_word": row["candidate_word"],
                "lifecycle_state": row["lifecycle_state"],
                "pressure_channels": row["pressure_channels"],
                "negative_evidence_channels": row["negative_evidence_channels"],
                "observation_count": row["observation_count"],
            }
            for row in ranked
        ]),
        "physical_basis_claim_allowed": False,
        "folding_problem_solved": False,
    }


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


def _present(value: float) -> bool:
    return float(value) > 0.0


def _dominates(value: float, *others: float) -> bool:
    return float(value) > max((float(other) for other in others), default=0.0)


def _not_dominated_by(value: float, *others: float) -> bool:
    return float(value) >= max((float(other) for other in others), default=0.0)


def _secondary_propensity(aa: str, window: str) -> str:
    helix_score = sum(1.0 if residue in HELIX else 0.0 for residue in window) / len(window)
    beta_score = sum(1.0 if residue in BETA else 0.0 for residue in window) / len(window)
    breaker_score = sum(1.0 if residue in BREAKERS else 0.0 for residue in window) / len(window)
    if _dominates(breaker_score, helix_score, beta_score):
        return "coil_or_disorder"
    if _dominates(helix_score, beta_score):
        return "helix_prone"
    if _dominates(beta_score, helix_score):
        return "beta_prone"
    if aa in POLAR or aa in LOW_COMPLEXITY_BIASED:
        return "coil_or_disorder"
    return "mixed"


def _periodic_hydrophobic_signal(sequence: str, period: int = 7) -> float:
    offsets = []
    for offset in range(period):
        positions = [sequence[index] for index in range(offset, len(sequence), period)]
        if len(positions) < 2:
            continue
        offsets.append(sum(1 for residue in positions if residue in HYDROPHOBIC) / len(positions))
    return round(max(offsets, default=0.0), 6)


def _repeat_signature(sequence: str, kmer_size: int = 4) -> float:
    kmers = [sequence[index : index + kmer_size] for index in range(0, len(sequence) - kmer_size + 1)]
    repeated_positions = sum(1 for kmer in kmers if kmers.count(kmer) > 1)
    return round(repeated_positions / len(kmers), 6) if kmers else 0.0


def _longest_run_fraction(sequence: str, alphabet: frozenset[str]) -> float:
    longest = 0
    current = 0
    for residue in sequence:
        if residue in alphabet:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return round(longest / max(1, len(sequence)), 6)


def _spread_n_terminal_hydrophobes_for_control(sequence: str, window_size: int = 35) -> str:
    window = sequence[:window_size]
    rest = sequence[window_size:]
    hydrophobic = [aa for aa in window if aa in HYDROPHOBIC]
    other = [aa for aa in window if aa not in HYDROPHOBIC]
    slots: list[str | None] = [None] * len(window)
    positions = list(range(0, len(window), 2)) + list(range(1, len(window), 2))
    for position, residue in zip(positions, hydrophobic):
        slots[position] = residue
    other_iter = iter(other)
    for index, residue in enumerate(slots):
        if residue is None:
            slots[index] = next(other_iter)
    return "".join(str(residue) for residue in slots) + rest


def _mask_n_terminal_hydrophobes_for_control(sequence: str, window_size: int = 35) -> str:
    window = sequence[:window_size]
    rest = sequence[window_size:]
    return "".join("S" if aa in HYDROPHOBIC else aa for aa in window) + rest


def _n_terminal_hydrophobic_signal(sequence: str) -> float:
    window = sequence[: min(len(sequence), 35)]
    if not window:
        return 0.0
    hydrophobic = sum(1 for residue in window if residue in HYDROPHOBIC) / len(window)
    positive = sum(1 for residue in window[:8] if residue in POSITIVE) / max(1, len(window[:8]))
    hydrophobic_core_run = _longest_run_fraction(window, HYDROPHOBIC)
    cleavage_tail = window[min(len(window), 18) :]
    cleavage_tail_small_or_polar = (
        sum(1 for residue in cleavage_tail if residue in POLAR or residue in {"A", "G", "S"}) / len(cleavage_tail)
        if cleavage_tail
        else 0.0
    )
    return bounded(_avg([hydrophobic, positive, hydrophobic_core_run, cleavage_tail_small_or_polar]))


def _sequence_from_field(sequence_field: dict[str, Any]) -> str:
    return "".join(row["residue_identity"] for row in sequence_field["residues"])


def _n_terminal_signal_beats_counterfactuals(sequence: str) -> bool:
    real = _n_terminal_hydrophobic_signal(sequence)
    shuffled = _n_terminal_hydrophobic_signal(_spread_n_terminal_hydrophobes_for_control(sequence))
    masked = _n_terminal_hydrophobic_signal(_mask_n_terminal_hydrophobes_for_control(sequence))
    return real > shuffled and real > masked


def _cysteine_topology_label(sequence_field: dict[str, Any]) -> str:
    metrics = sequence_field["global_metrics"]
    cysteine_count = int(metrics.get("cysteine_count", 0))
    if cysteine_count < 2:
        return "cysteine_pairing_signal_absent"
    if cysteine_count % 2 != 0:
        return "odd_cysteine_pairing_conflict"
    cysteine_segments = [row for row in sequence_field["segments"] if row.get("cysteine_density", 0.0) > 0.0]
    if len(cysteine_segments) >= 2:
        return "paired_cysteine_topology_visible"
    return "local_cysteine_pairing_signal_visible"


def _has_paired_cysteine_topology(sequence_field: dict[str, Any]) -> bool:
    return "visible" in _cysteine_topology_label(sequence_field)


def _signal_peptide_routing_label(sequence_field: dict[str, Any]) -> str:
    metrics = sequence_field["global_metrics"]
    segments = sequence_field["segments"]
    sequence = _sequence_from_field(sequence_field)
    n_terminal = float(metrics.get("n_terminal_hydrophobic_signal", 0.0))
    internal_membrane = max((row["membrane_density"] for row in segments[2:]), default=0.0)
    first_membrane = max((row["membrane_density"] for row in segments[:2]), default=0.0)
    if n_terminal <= 0.0 and first_membrane <= 0.0:
        return "n_terminal_signal_absent"
    if internal_membrane > max(n_terminal, first_membrane):
        return "internal_true_tm_signal_dominates"
    if not _n_terminal_signal_beats_counterfactuals(sequence):
        return "n_terminal_signal_counterfactual_not_dominant"
    if n_terminal >= internal_membrane or first_membrane >= internal_membrane:
        return "n_terminal_signal_over_internal_tm_visible"
    return "signal_tm_boundary_unresolved"


def _has_signal_peptide_routing_topology(sequence_field: dict[str, Any]) -> bool:
    return "visible" in _signal_peptide_routing_label(sequence_field)


def _residue_mark(sequence: str, index0: int) -> dict[str, Any]:
    aa = sequence[index0]
    window = _window(sequence, index0, radius=4)
    unique_fraction = len(set(window)) / len(window)
    biased_fraction = sum(1 for residue in window if residue in LOW_COMPLEXITY_BIASED) / len(window)
    hydrophobic_fraction = sum(1 for residue in window if residue in HYDROPHOBIC) / len(window)
    charged_fraction = sum(1 for residue in window if residue in POSITIVE or residue in NEGATIVE) / len(window)
    aromatic_fraction = sum(1 for residue in window if residue in AROMATIC) / len(window)
    pro_gly_fraction = sum(1 for residue in window if residue in BREAKERS) / len(window)
    low_complexity = _dominates(biased_fraction, unique_fraction) or _dominates(pro_gly_fraction, unique_fraction)
    disorder = _attenuate_by_competing_pressure(
        _endogenous_support(biased_fraction, pro_gly_fraction, charged_fraction),
        hydrophobic_fraction,
    )
    membrane = _attenuate_by_competing_pressure(hydrophobic_fraction, charged_fraction, pro_gly_fraction)
    interface = _endogenous_support(hydrophobic_fraction, aromatic_fraction, charged_fraction)
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
        "domain_boundary_tendency": _endogenous_support(pro_gly_fraction, charged_fraction),
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
            "cysteine_density": _avg(1.0 if row["cysteine_mark"] else 0.0 for row in segment_residues),
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
        "cysteine_count": sum(1 for row in residues if row["cysteine_mark"]),
        "cysteine_density": _avg(1.0 if row["cysteine_mark"] else 0.0 for row in residues),
        "histidine_density": _avg(1.0 if row["residue_identity"] == "H" else 0.0 for row in residues),
        "heptad_hydrophobic_periodicity": _periodic_hydrophobic_signal(sequence),
        "repeat_signature": _repeat_signature(sequence),
        "n_terminal_hydrophobic_signal": _n_terminal_hydrophobic_signal(sequence),
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
    return _bounded_span(sequence_length, 1, sequence_length)


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
        early_span = _route_span(length, "full_sequence")
        late_span = _route_span(length, "full_sequence")
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
        early_span = _route_span(length, "alpha_core")
        late_span = _route_span(length, "early_n")
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


def _token_hits(text: str, tokens: list[str]) -> list[str]:
    return sorted({token for token in tokens if token.lower() in text})


def _sequence_signal_label(grammar: str, sequence_field: dict[str, Any]) -> str:
    metrics = sequence_field["global_metrics"]
    segments = sequence_field["segments"]
    local_hydrophobic = max((row["hydrophobic_density"] for row in segments), default=0.0)
    local_interface = max((row["interface_density"] for row in segments), default=0.0)
    local_membrane = max((row["membrane_density"] for row in segments), default=0.0)
    local_beta = max((row["beta_propensity_density"] for row in segments), default=0.0)
    local_disorder = max((row["disorder_density"] + row["low_complexity_density"] for row in segments), default=0.0)
    cysteine_count = int(metrics.get("cysteine_count", 0))
    if grammar in {"coiled_coil_register", "coiled_coil_register_topology"}:
        if _dominates(metrics.get("heptad_hydrophobic_periodicity", 0.0), metrics["hydrophobic_density"]):
            return "heptad_hydrophobic_periodicity_visible"
        return "heptad_signal_not_dominant"
    if grammar == "repeat_solenoid_topology":
        if _present(metrics.get("repeat_signature", 0.0)):
            return "repeat_or_long_solenoid_signature_visible"
        return "repeat_signature_not_dominant"
    if grammar == "disulfide_secretory_redox_context":
        return _cysteine_topology_label(sequence_field)
    if grammar == "secretory_disulfide_redox_topology":
        if cysteine_count >= 2:
            return _cysteine_topology_label(sequence_field)
        return "secretory_evidence_without_cysteine_pair_sequence_support"
    if grammar == "signal_peptide_vs_true_tm_routing":
        return _signal_peptide_routing_label(sequence_field)
    if grammar == "signal_peptide_vs_true_TM":
        if _not_dominated_by(metrics.get("n_terminal_hydrophobic_signal", 0.0), local_membrane):
            return "n_terminal_hydrophobic_signal_visible"
        return "signal_peptide_sequence_signal_not_dominant"
    if grammar == "knotted_topology":
        return "topological_evidence_required_beyond_sequence"
    if grammar == "globular_closure":
        if _not_dominated_by(local_hydrophobic, metrics["hydrophobic_density"]) and _not_dominated_by(metrics["hydrophobic_density"], metrics["mean_disorder"]):
            return "compact_core_sequence_signal_visible"
        return "globular_sequence_signal_weak"
    if grammar in {"intrinsic_disorder_phase_separation", "disorder_boundary_and_fold_upon_binding", "fold_upon_binding_disorder"}:
        if _not_dominated_by(local_disorder, local_hydrophobic):
            return "disorder_or_low_complexity_sequence_signal_visible"
        return "disorder_sequence_signal_weak"
    if grammar == "beta_closure_topology":
        if _not_dominated_by(local_beta, metrics["beta_propensity_density"]):
            return "beta_strand_propensity_signal_visible"
        return "beta_sequence_signal_weak"
    if grammar == "multidomain_allosteric_architecture":
        if _not_dominated_by(local_interface, metrics["mean_interface"]):
            return "long_modular_interface_sequence_signal_visible"
        return "multidomain_sequence_signal_weak"
    if grammar == "membrane_multidomain_folding_proteostasis":
        if _not_dominated_by(local_membrane, metrics["mean_membrane"]):
            return "membrane_segment_sequence_signal_visible"
        return "membrane_sequence_signal_weak"
    if grammar in {"cofactor_ligand_assisted_stabilization", "metal_cluster_and_ligand_locked_basin"}:
        if _not_dominated_by(local_interface, metrics["mean_interface"]) or _present(metrics.get("histidine_density", 0.0)):
            return "pocket_or_coordination_sequence_signal_visible"
        return "cofactor_sequence_signal_weak"
    if grammar in {"assembly_required_folding", "oligomerization_controlled_folding"}:
        if _not_dominated_by(local_interface, metrics["mean_interface"]):
            return "interface_sequence_signal_visible"
        return "assembly_sequence_signal_weak"
    return "sequence_signal_not_specific"


def _sequence_signal_is_supportive(label: str) -> bool:
    return "visible" in label


def _trajectory_supports_grammar(grammar: str, trajectory: dict[str, Any]) -> list[str]:
    final = trajectory["final_state_summary"]
    expected = {
        "globular_closure": ["contact_probability", "segment_compaction"],
        "intrinsic_disorder_phase_separation": ["phase_prone_low_complexity", "disorder_order_balance"],
        "disorder_boundary_and_fold_upon_binding": ["IDR_boundary", "fold_upon_binding_region"],
        "beta_closure_topology": ["closed_beta_topology", "strand_register", "beta_sheet_closure"],
        "multidomain_allosteric_architecture": ["multidomain_allostery", "domain_boundary", "interdomain_lock"],
        "secretory_disulfide_redox_topology": [
            "secretory_redox_context",
            "disulfide_pairing_topology",
            "cysteine_pairing_constraint",
            "extracellular_stabilized_fold",
        ],
        "signal_peptide_vs_true_tm_routing": [
            "signal_peptide_routing_context",
            "cleavage_site_context",
            "n_terminal_secretory_hydrophobic_patch",
            "secretory_lumenal_routing",
        ],
        "coiled_coil_register_topology": [
            "heptad_register_context",
            "hydrophobic_repeat_phase",
            "oligomeric_coiled_coil_core",
        ],
        "repeat_solenoid_topology": [
            "repeat_unit_context",
            "solenoid_axis_context",
            "global_repeat_topology",
        ],
        "knotted_topology": [
            "knot_core_context",
            "threading_loop_context",
            "topological_closure_constraint",
        ],
        "membrane_multidomain_folding_proteostasis": ["proteostasis_routing"],
        "metamorphic_fold_switching": ["state_basin_occupancy"],
        "short_region_host_interface_hijacking": ["interface_readiness"],
        "fold_upon_binding_disorder": ["interface_readiness", "disorder_order_balance"],
        "cofactor_ligand_assisted_stabilization": ["interface_readiness", "contact_probability"],
        "metal_cluster_and_ligand_locked_basin": ["metal_cluster_geometry", "ligand_locked_basin", "coordination_shell_integrity"],
        "assembly_required_folding": ["partner_completed_core", "assembly_required_core"],
        "oligomerization_controlled_folding": ["interface_readiness"],
    }
    supported = []
    for key in expected.get(grammar, []):
        value = final.get(key)
        if isinstance(value, dict) and any(float(item) > 0.0 for item in value.values()):
            supported.append(key)
        elif isinstance(value, (int, float)) and float(value) > 0.0:
            supported.append(key)
    return supported


def _rank_readout(readout: dict[str, Any]) -> tuple[int, int, int, str]:
    return (
        int(readout["internal_view_count"]),
        len(readout["evidence_hits"]),
        0 if readout["grammar_status"] == "candidate_missing_word" else 1,
        readout["grammar"],
    )


def _learned_grammar_readout(
    *,
    grammar: str,
    text: str,
    sequence_field: dict[str, Any],
    selected_mechanism: str,
    operator_field: dict[str, Any],
    trajectory: dict[str, Any],
) -> dict[str, Any]:
    evidence_hits = _token_hits(text, SELF_DECISION_LEARNED_GRAMMAR_FAMILIES.get(grammar, []))
    sequence_label = _sequence_signal_label(grammar, sequence_field)
    trajectory_hits = _trajectory_supports_grammar(grammar, trajectory) if grammar == selected_mechanism else []
    view_sources = []
    if evidence_hits:
        view_sources.append("evidence_family")
    if _sequence_signal_is_supportive(sequence_label):
        view_sources.append("sequence_signature")
    if grammar == selected_mechanism and operator_field["active_operator_count"] > 0:
        view_sources.append("operator_readout")
    if trajectory_hits:
        view_sources.append("trajectory_readout")
    return {
        "grammar": grammar,
        "grammar_status": "learned",
        "evidence_hits": evidence_hits,
        "sequence_support": sequence_label,
        "operator_support": "selected_operator_field_agrees" if grammar == selected_mechanism and operator_field["active_operator_count"] > 0 else "not_selected_operator",
        "trajectory_support": trajectory_hits,
        "view_sources": view_sources,
        "internal_view_count": len(view_sources),
    }


def _candidate_grammar_readout(
    *,
    grammar: str,
    spec: dict[str, Any],
    text: str,
    sequence_field: dict[str, Any],
    selected_mechanism: str,
) -> dict[str, Any]:
    evidence_hits = _token_hits(text, spec["evidence_tokens"])
    sequence_label = _sequence_signal_label(grammar, sequence_field)
    view_sources = []
    if evidence_hits:
        view_sources.append("evidence_family")
    if _sequence_signal_is_supportive(sequence_label):
        view_sources.append("sequence_signature")
    if selected_mechanism in spec.get("competes_with", []):
        view_sources.append("wrong_grammar_competition")
    return {
        "grammar": grammar,
        "grammar_status": spec["grammar_status"],
        "evidence_hits": evidence_hits,
        "sequence_support": sequence_label,
        "operator_support": "no_learned_operator_available",
        "trajectory_support": [],
        "view_sources": view_sources,
        "internal_view_count": len(view_sources),
    }


def _mechanism_competition_readouts(
    *,
    text: str,
    sequence_field: dict[str, Any],
    mechanism: dict[str, Any],
    operator_field: dict[str, Any],
    trajectory: dict[str, Any],
) -> list[dict[str, Any]]:
    selected = mechanism["mechanism_class"]
    readouts = [
        _learned_grammar_readout(
            grammar=grammar,
            text=text,
            sequence_field=sequence_field,
            selected_mechanism=selected,
            operator_field=operator_field,
            trajectory=trajectory,
        )
        for grammar in MECHANISM_CLASSES
        if grammar != "insufficient_evidence_clean_abstain"
    ]
    readouts.extend(
        _candidate_grammar_readout(
            grammar=grammar,
            spec=spec,
            text=text,
            sequence_field=sequence_field,
            selected_mechanism=selected,
        )
        for grammar, spec in SELF_DECISION_CANDIDATE_GRAMMARS.items()
    )
    return sorted(readouts, key=_rank_readout, reverse=True)


def _remove_family_tokens(value: Any, tokens: list[str]) -> Any:
    if isinstance(value, dict):
        return {
            key: _remove_family_tokens(item, tokens)
            for key, item in value.items()
            if not str(key).startswith("withheld")
        }
    if isinstance(value, list):
        return [_remove_family_tokens(item, tokens) for item in value]
    if isinstance(value, str):
        masked = value
        for token in sorted(tokens, key=len, reverse=True):
            masked = masked.replace(token, " ").replace(token.lower(), " ").replace(token.upper(), " ")
        return " ".join(masked.split())
    return value


def _masking_probe(
    *,
    sequence_field: dict[str, Any],
    evidence_manifest: dict[str, Any],
    sources: list[dict[str, Any]],
    top_mechanism: str,
) -> list[dict[str, Any]]:
    rows = []
    family_map = {
        **SELF_DECISION_LEARNED_GRAMMAR_FAMILIES,
        **{grammar: spec["evidence_tokens"] for grammar, spec in SELF_DECISION_CANDIDATE_GRAMMARS.items()},
    }
    text = _allowed_source_text(sources, evidence_manifest)
    for family, tokens in family_map.items():
        hits = _token_hits(text, tokens)
        if not hits:
            continue
        masked_sources = [_remove_family_tokens(source, tokens) for source in sources]
        masked_gate = evidence_boundary_gate(masked_sources)
        masked_mechanism = select_mechanism_grammar(
            sequence_field=sequence_field,
            evidence_manifest=masked_gate,
            sources=masked_sources,
        )
        rows.append({
            "masked_family": family,
            "removed_signal_count": len(hits),
            "mechanism_after_mask": masked_mechanism["mechanism_class"],
            "selected_mechanism_preserved": masked_mechanism["mechanism_class"] == top_mechanism,
            "interpretation": (
                "mechanism_defining_signal_removed"
                if family == top_mechanism
                else "stable_under_nondefining_mask"
                if masked_mechanism["mechanism_class"] == top_mechanism
                else "unstable_under_competing_mask"
            ),
        })
    return rows


def _wrong_grammar_challenge(
    *,
    sequence_field: dict[str, Any],
    evidence_manifest: dict[str, Any],
    sources: list[dict[str, Any]],
    mechanism: dict[str, Any],
    readouts: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    natural = mechanism["natural_mechanism_class"]
    challengers = [
        row["grammar"]
        for row in readouts
        if row["grammar_status"] == "learned" and row["grammar"] not in {natural, "insufficient_evidence_clean_abstain"}
    ][:3]
    rows = []
    for challenger in challengers:
        forced = select_mechanism_grammar(
            sequence_field=sequence_field,
            evidence_manifest=evidence_manifest,
            sources=sources,
            forced_grammar=challenger,
        )
        forced_operator = build_operator_field(
            sequence_field=sequence_field,
            mechanism=forced,
            evidence_manifest=evidence_manifest,
        )
        forced_trajectory = simulate_operator_trajectory(
            sequence_field=sequence_field,
            operator_field=forced_operator,
            mechanism_class=forced["mechanism_class"],
        )
        rows.append({
            "challenger_grammar": challenger,
            "forced_grammar_rejected": forced["forced_grammar_rejected"],
            "forced_result_mechanism": forced["mechanism_class"],
            "challenger_operator_count": forced_operator["active_operator_count"],
            "challenger_trajectory_support": _trajectory_supports_grammar(challenger, forced_trajectory),
            "interpretation": "wrong_grammar_failed" if forced["forced_grammar_rejected"] else "wrong_grammar_competes",
        })
    return rows


def _counterfactual_sequence_controls(sequence_field: dict[str, Any], grammar: str) -> dict[str, Any]:
    sequence = "".join(row["residue_identity"] for row in sequence_field["residues"])
    shuffled = build_sequence_field(shuffled_sequence(sequence))
    masked = build_sequence_field("A" * len(sequence))
    real_label = _sequence_signal_label(grammar, sequence_field)
    shuffled_label = _sequence_signal_label(grammar, shuffled)
    masked_label = _sequence_signal_label(grammar, masked)
    real_support = _sequence_signal_is_supportive(real_label)
    shuffled_support = _sequence_signal_is_supportive(shuffled_label)
    masked_support = _sequence_signal_is_supportive(masked_label)
    if real_support and not (shuffled_support and masked_support):
        interpretation = "counterfactuals_weaken_signal"
    elif real_support:
        interpretation = "counterfactuals_do_not_separate_sequence_signal"
    else:
        interpretation = "context_or_operator_driven_not_sequence_only"
    return {
        "real_sequence_support": real_label,
        "composition_preserved_shuffle_support": shuffled_label,
        "region_masked_sequence_support": masked_label,
        "interpretation": interpretation,
    }


def _mechanism_state_signature(mechanism_class: str, trajectory: dict[str, Any]) -> str:
    final = trajectory["final_state_summary"]
    if mechanism_class == "globular_closure":
        if float(final.get("contact_probability", 0.0)) > 0.0 or float(final.get("segment_compaction", 0.0)) > 0.0:
            return "globular_contact_closure"
        return "globular_uninitialized"
    if mechanism_class == "intrinsic_disorder_phase_separation":
        basin = final.get("state_basin_occupancy", {})
        return max(basin, key=basin.get) if isinstance(basin, dict) and basin else "disorder_no_basin"
    if mechanism_class == "disorder_boundary_and_fold_upon_binding":
        if float(final.get("fold_upon_binding_region", 0.0)) >= float(final.get("IDR_boundary", 0.0)):
            return "fold_upon_binding_boundary"
        return "idr_boundary_ensemble"
    if mechanism_class == "beta_closure_topology":
        return "closed_beta_register" if float(final.get("strand_register", 0.0)) > 0.0 else "beta_unregistered"
    if mechanism_class == "multidomain_allosteric_architecture":
        candidates = {
            "interdomain_lock": final.get("interdomain_lock", 0.0),
            "allosteric_basin_shift": final.get("allosteric_basin_shift", 0.0),
            "domain_reorientation": final.get("domain_reorientation", 0.0),
            "domain_swapping": final.get("domain_swapping", 0.0),
        }
        return max(candidates, key=candidates.get)
    if mechanism_class == "secretory_disulfide_redox_topology":
        if float(final.get("disulfide_pairing_topology", 0.0)) >= float(final.get("redox_mispaired_frustration", 0.0)):
            return "secretory_disulfide_pairing_topology"
        return "redox_mispaired_frustration"
    if mechanism_class == "signal_peptide_vs_true_tm_routing":
        signal = max(
            float(final.get("signal_peptide_routing_context", 0.0)),
            float(final.get("secretory_lumenal_routing", 0.0)),
            float(final.get("cleavage_site_context", 0.0)),
        )
        membrane = float(final.get("membrane_insertion_routing", 0.0))
        anchor = float(final.get("signal_anchor_ambiguity", 0.0))
        if anchor >= signal and anchor >= membrane:
            return "signal_anchor_ambiguity"
        if membrane >= signal:
            return "true_tm_membrane_insertion_route"
        if signal > 0.0:
            return "cleavable_signal_peptide_secretory_route"
        return "signal_peptide_route_unresolved"
    if mechanism_class == "coiled_coil_register_topology":
        if float(final.get("register_shift_frustration", 0.0)) > max(
            float(final.get("heptad_register_context", 0.0)),
            float(final.get("oligomeric_coiled_coil_core", 0.0)),
        ):
            return "coiled_coil_register_shift_conflict"
        if float(final.get("heptad_register_context", 0.0)) > 0.0 and float(final.get("oligomeric_coiled_coil_core", 0.0)) > 0.0:
            return "heptad_registered_coiled_coil_core"
        return "coiled_coil_register_unresolved"
    if mechanism_class == "repeat_solenoid_topology":
        if float(final.get("repeat_boundary_frustration", 0.0)) > float(final.get("global_repeat_topology", 0.0)):
            return "repeat_solenoid_boundary_conflict"
        if float(final.get("repeat_unit_context", 0.0)) > 0.0 and float(final.get("solenoid_axis_context", 0.0)) > 0.0:
            return "phase_aligned_repeat_solenoid_axis"
        return "repeat_solenoid_axis_unresolved"
    if mechanism_class == "knotted_topology":
        if float(final.get("unknotted_decoy_dominance", 0.0)) > float(final.get("topological_closure_constraint", 0.0)):
            return "unknotted_decoy_dominates"
        if float(final.get("threading_loop_context", 0.0)) > 0.0 and float(final.get("topological_closure_constraint", 0.0)) > 0.0:
            return "threaded_knot_topology"
        return "knotted_topology_unresolved"
    if mechanism_class == "membrane_multidomain_folding_proteostasis":
        return "membrane_proteostasis_route" if float(final.get("proteostasis_routing", 0.0)) > 0.0 else "membrane_route_unresolved"
    if mechanism_class == "metal_cluster_and_ligand_locked_basin":
        return "ligand_locked_basin" if float(final.get("ligand_locked_basin", 0.0)) > 0.0 else "metal_ligand_unlocked"
    if mechanism_class == "assembly_required_folding":
        return "partner_completed_core" if float(final.get("partner_completed_core", 0.0)) > 0.0 else "assembly_uncompleted"
    if mechanism_class == "oligomerization_controlled_folding":
        return "oligomer_interface_ready" if float(final.get("interface_readiness", 0.0)) > 0.0 else "oligomer_interface_unready"
    return mechanism_class


def _with_operator_strength_order(operator_field: dict[str, Any], strengths: tuple[float, ...]) -> dict[str, Any]:
    field = deepcopy(operator_field)
    for row, strength in zip(field["operators"], strengths):
        row["activation_strength"] = bounded(strength)
    return field


def _without_operator_index(operator_field: dict[str, Any], index: int) -> dict[str, Any]:
    field = deepcopy(operator_field)
    field["operators"] = [row for row_index, row in enumerate(field["operators"]) if row_index != index]
    field["operator_names"] = [row["operator"] for row in field["operators"]]
    field["active_operator_count"] = len(field["operators"])
    return field


def _operator_basis_stability_probe(
    *,
    sequence_field: dict[str, Any],
    operator_field: dict[str, Any],
    mechanism_class: str,
    baseline_trajectory: dict[str, Any],
) -> dict[str, Any]:
    if mechanism_class == "insufficient_evidence_clean_abstain" or not operator_field["operators"]:
        return {
            "operator_basis_stability": "no_operator_basis_to_probe",
            "coefficient_probe_mode": "no_static_scale_range",
            "baseline_signature": mechanism_class,
            "coefficient_permutation_rows": [],
            "operator_ablation_rows": [],
        }
    baseline_signature = _mechanism_state_signature(mechanism_class, baseline_trajectory)
    observed_strengths = tuple(row["activation_strength"] for row in operator_field["operators"])
    unique_strength_orders = sorted(set(permutations(observed_strengths)))
    declared_grammar_operators = set(GRAMMAR_RULES.get(mechanism_class, {}).get("operators", []))
    permutation_rows = []
    for order_index, strengths in enumerate(unique_strength_orders, start=1):
        perturbed_field = _with_operator_strength_order(operator_field, strengths)
        perturbed_trajectory = simulate_operator_trajectory(
            sequence_field=sequence_field,
            operator_field=perturbed_field,
            mechanism_class=mechanism_class,
        )
        signature = _mechanism_state_signature(mechanism_class, perturbed_trajectory)
        permutation_rows.append({
            "probe_id": f"observed_operator_strength_permutation_{order_index:02d}",
            "probe_kind": "endogenous_observed_operator_coefficient_permutation",
            "external_scale_used": None,
            "selected_signature_preserved": signature == baseline_signature,
            "baseline_signature": baseline_signature,
            "perturbed_signature": signature,
            "operator_strength_order": strengths,
        })
    ablation_rows = []
    for operator_index, operator in enumerate(operator_field["operators"]):
        ablated_field = _without_operator_index(operator_field, operator_index)
        perturbed_trajectory = simulate_operator_trajectory(
            sequence_field=sequence_field,
            operator_field=ablated_field,
            mechanism_class=mechanism_class,
        )
        signature = _mechanism_state_signature(mechanism_class, perturbed_trajectory)
        ablation_rows.append({
            "probe_id": f"remove_operator_{operator_index + 1:02d}_{operator['operator']}",
            "probe_kind": "single_operator_ablation_map",
            "selected_signature_preserved": signature == baseline_signature,
            "baseline_signature": baseline_signature,
            "perturbed_signature": signature,
            "removed_operator": operator["operator"],
            "removed_state_variable": operator["state_variable"],
            "removed_operator_declared_by_selected_grammar": operator["operator"] in declared_grammar_operators,
        })
    permutation_stable = all(row["selected_signature_preserved"] for row in permutation_rows)
    nondefining_ablation_sensitive = [
        row
        for row in ablation_rows
        if not row["removed_operator_declared_by_selected_grammar"] and not row["selected_signature_preserved"]
    ]
    defining_ablation_sensitive = [
        row
        for row in ablation_rows
        if row["removed_operator_declared_by_selected_grammar"] and not row["selected_signature_preserved"]
    ]
    ablation_stable = not nondefining_ablation_sensitive
    if permutation_stable:
        stability = "stable_under_endogenous_operator_basis_probe"
    elif ablation_stable:
        stability = "definition_sensitive_under_semantic_operator_basis_probe"
    else:
        stability = "nondefining_operator_basis_sensitive"
    return {
        "operator_basis_stability": stability,
        "coefficient_probe_mode": "endogenous_observed_operator_permutations_no_static_scale_range",
        "coefficient_scale_values_used": [],
        "baseline_signature": baseline_signature,
        "coefficient_assignment_cross_semantic_role_sensitive": not permutation_stable,
        "definition_sensitive_operator_ablation_observed": bool(defining_ablation_sensitive),
        "nondefining_operator_ablation_sensitive": bool(nondefining_ablation_sensitive),
        "semantic_operator_ablation_stable": ablation_stable,
        "coefficient_permutation_rows": permutation_rows,
        "operator_ablation_rows": ablation_rows,
    }


def _temporal_binding_probe(mechanism_class: str, trajectory: dict[str, Any]) -> dict[str, Any]:
    if mechanism_class == "insufficient_evidence_clean_abstain":
        return {
            "temporal_binding": "no_mechanism_trajectory_to_bind",
            "signature_path": [],
            "observable_trends": [],
        }
    signature_path = [
        {
            "timepoint": row["timepoint"],
            "signature": _mechanism_state_signature(mechanism_class, {"final_state_summary": row}),
        }
        for row in trajectory["timepoints"]
    ]
    support_keys = _trajectory_supports_grammar(mechanism_class, trajectory)
    trends = []
    first = trajectory["timepoints"][0]
    last = trajectory["timepoints"][-1]
    for key in support_keys:
        if not isinstance(first.get(key), (int, float)) or not isinstance(last.get(key), (int, float)):
            continue
        if last[key] > first[key]:
            trend = "increases_across_simulated_path"
        elif last[key] == first[key]:
            trend = "unchanged_across_simulated_path"
        else:
            trend = "decreases_across_simulated_path"
        trends.append({
            "observable": key,
            "first_value": first[key],
            "final_value": last[key],
            "trend": trend,
        })
    signatures = [row["signature"] for row in signature_path]
    if signatures and all(signature == signatures[0] for signature in signatures):
        temporal_binding = "selected_observables_temporally_coherent"
    elif any(row["trend"] == "decreases_across_simulated_path" and last.get(row["observable"], 0.0) <= 0.0 for row in trends):
        temporal_binding = "selected_observable_temporal_conflict"
    elif trends:
        temporal_binding = "selected_observables_temporally_coherent"
    else:
        temporal_binding = "no_selected_numeric_temporal_observable"
    return {
        "temporal_binding": temporal_binding,
        "signature_path": signature_path,
        "observable_trends": trends,
    }


@lru_cache(maxsize=1)
def _physical_backend_probe() -> dict[str, Any]:
    openmm_spec = importlib.util.find_spec("openmm")
    if openmm_spec is None:
        return {
            "kind": "PROTEIN_ESPERANTO_PHYSICAL_BACKEND_PROBE_v0",
            "backend": "openmm",
            "backend_available": False,
            "backend_version": None,
            "forcefield_api_available": False,
        }
    try:
        import openmm  # type: ignore
        import openmm.app as openmm_app  # type: ignore

        version = getattr(openmm, "__version__", "unknown")
        forcefield_api_available = hasattr(openmm_app, "ForceField")
    except Exception as exc:  # pragma: no cover - defensive environment probe.
        return {
            "kind": "PROTEIN_ESPERANTO_PHYSICAL_BACKEND_PROBE_v0",
            "backend": "openmm",
            "backend_available": False,
            "backend_version": None,
            "forcefield_api_available": False,
            "backend_error": type(exc).__name__,
        }
    return {
        "kind": "PROTEIN_ESPERANTO_PHYSICAL_BACKEND_PROBE_v0",
        "backend": "openmm",
        "backend_available": True,
        "backend_version": version,
        "forcefield_api_available": forcefield_api_available,
    }


def _physical_calibration_input_summary(physical_calibration_inputs: dict[str, Any] | None) -> dict[str, Any]:
    if not isinstance(physical_calibration_inputs, dict):
        return {
            "real_physical_calibration_inputs_used": False,
            "physical_calibration_input_status": "no_real_physical_calibration_inputs",
            "real_physical_calibration_kind": None,
            "real_physical_calibration_row_count": 0,
            "real_physical_calibration_hash": None,
            "source_coordinate_database": None,
            "target_native_excluded_from_calibration": None,
            "target_native_contacts_used_before_prediction": None,
            "coordinate_truth_used_as_prediction_input": None,
            "leave_one_target_out_calibration": None,
            "calibration_input_type": None,
            "fold_class_coverage": [],
            "observable_families": [],
        }
    required_truth_boundaries = {
        "target_native_excluded_from_calibration": True,
        "target_native_contacts_used_before_prediction": False,
        "coordinate_truth_used_as_prediction_input": False,
        "leave_one_target_out_calibration": True,
    }
    boundary_ok = all(physical_calibration_inputs.get(key) == value for key, value in required_truth_boundaries.items())
    row_count = int(physical_calibration_inputs.get("row_count") or 0)
    is_real_coordinate = physical_calibration_inputs.get("source_coordinate_database") == "RCSB_PDB"
    inputs_used = boundary_ok and row_count > 0 and is_real_coordinate
    if inputs_used:
        status = "real_coordinate_calibration_inputs_loaded_truth_boundary_preserved"
    elif row_count > 0:
        status = "physical_calibration_inputs_rejected_truth_boundary_or_source_mismatch"
    else:
        status = "physical_calibration_inputs_empty"
    return {
        "real_physical_calibration_inputs_used": inputs_used,
        "physical_calibration_input_status": status,
        "real_physical_calibration_kind": physical_calibration_inputs.get("kind"),
        "real_physical_calibration_row_count": row_count,
        "real_physical_calibration_hash": physical_calibration_inputs.get("calibration_hash"),
        "source_dataset": physical_calibration_inputs.get("source_dataset"),
        "source_dataset_sha256": physical_calibration_inputs.get("source_dataset_sha256"),
        "source_coordinate_database": physical_calibration_inputs.get("source_coordinate_database"),
        "target_native_excluded_from_calibration": physical_calibration_inputs.get("target_native_excluded_from_calibration"),
        "target_native_contacts_used_before_prediction": physical_calibration_inputs.get("target_native_contacts_used_before_prediction"),
        "coordinate_truth_used_as_prediction_input": physical_calibration_inputs.get("coordinate_truth_used_as_prediction_input"),
        "leave_one_target_out_calibration": physical_calibration_inputs.get("leave_one_target_out_calibration"),
        "calibration_input_type": physical_calibration_inputs.get("calibration_input_type"),
        "fold_class_coverage": list(physical_calibration_inputs.get("fold_class_coverage") or []),
        "observable_families": list(physical_calibration_inputs.get("observable_families") or []),
    }


def _physical_grounding_gate(
    *,
    sequence_field: dict[str, Any],
    evidence_manifest: dict[str, Any],
    mechanism: dict[str, Any],
    operator_field: dict[str, Any],
    trajectory: dict[str, Any],
    temporal_probe: dict[str, Any],
    operator_basis_probe: dict[str, Any],
    physical_calibration_inputs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    backend = _physical_backend_probe()
    calibration_summary = _physical_calibration_input_summary(physical_calibration_inputs)
    mechanism_class = mechanism["mechanism_class"]
    grammar = GRAMMAR_RULES.get(mechanism_class, {})
    operator_names = set(operator_field["operator_names"])
    grammar_operators = set(grammar.get("operators", []))
    operator_mapping = [
        {
            "operator": row["operator"],
            "state_variable": row["state_variable"],
            "declared_in_grammar": row["operator"] in grammar_operators,
            "physical_scope": "coarse_pressure_proxy_not_force_field_term",
        }
        for row in operator_field["operators"]
    ]
    if calibration_summary["real_physical_calibration_inputs_used"]:
        missing_physical_inputs = [
            "target_specific_physical_topology_and_environment_for_execution",
            "force_field_or_calibrated_coarse_potential_bound_to_target_execution",
            "calibrated_energy_units_temperature_and_integrator_protocol",
            "independent_dynamic_physical_observable_holdout",
        ]
    else:
        missing_physical_inputs = [
            "real_coordinate_calibration_inputs_with_truth_boundary",
            "atomistic_or_validated_coarse_coordinates",
            "force_field_or_calibrated_coarse_potential",
            "solvent_membrane_ligand_environment_when_relevant",
            "calibrated_energy_units_and_temperature_protocol",
            "independent_physical_observable_holdout",
        ]
    if evidence_manifest["coordinate_derived_source_count_before_prediction"] > 0:
        status = "blocked_coordinate_truth_before_seal"
    elif calibration_summary["real_physical_calibration_inputs_used"] and backend["backend_available"]:
        status = "real_physical_calibration_inputs_loaded_openmm_available_but_target_physical_execution_not_run"
    elif calibration_summary["real_physical_calibration_inputs_used"]:
        status = "real_physical_calibration_inputs_loaded_without_physical_backend"
    elif backend["backend_available"]:
        status = "openmm_available_but_physical_model_not_calibrated"
    else:
        status = "no_physical_backend_available"
    return {
        "kind": "PROTEIN_ESPERANTO_PHYSICAL_GROUNDING_GATE_v0",
        "physical_grounding_status": status,
        "physical_basis_claim_allowed": False,
        "physical_limitation_known_by_engine": True,
        "backend_probe": backend,
        "coarse_operator_layer": True,
        "coarse_operator_decision_scope": "mechanism_language_acceptance_not_physical_fold_solution",
        "coefficient_source": "heuristic_internal_operator_weights_not_physical_force_constants",
        "missing_physical_inputs": missing_physical_inputs,
        "physical_calibration_input_summary": calibration_summary,
        "real_physical_calibration_inputs_used": calibration_summary["real_physical_calibration_inputs_used"],
        "real_physical_calibration_kind": calibration_summary["real_physical_calibration_kind"],
        "real_physical_calibration_row_count": calibration_summary["real_physical_calibration_row_count"],
        "real_physical_calibration_hash": calibration_summary["real_physical_calibration_hash"],
        "target_native_excluded_from_calibration": calibration_summary["target_native_excluded_from_calibration"],
        "target_native_contacts_used_before_prediction": calibration_summary["target_native_contacts_used_before_prediction"],
        "coordinate_truth_used_as_prediction_input": calibration_summary["coordinate_truth_used_as_prediction_input"],
        "leave_one_target_out_calibration": calibration_summary["leave_one_target_out_calibration"],
        "operator_to_physical_proxy_map": operator_mapping,
        "grammar_declared_operator_coverage": sorted(operator_names.intersection(grammar_operators)),
        "grammar_missing_declared_operators": sorted(grammar_operators.difference(operator_names)),
        "temporal_binding": temporal_probe["temporal_binding"],
        "operator_basis_stability": operator_basis_probe["operator_basis_stability"],
        "folding_problem_solved": False,
    }


def _cross_view_binding_probe(selected_readout: dict[str, Any] | None) -> dict[str, Any]:
    if not selected_readout or selected_readout["grammar_status"] != "learned":
        return {
            "cross_view_binding": "no_learned_mechanism_to_bind",
            "bound_view_families": [],
            "self_required_view_families": [],
            "missing_view_families": [],
        }
    bound = set(selected_readout["view_sources"])
    rule = GRAMMAR_RULES.get(selected_readout["grammar"], {})
    required = []
    if selected_readout["evidence_hits"]:
        required.append("evidence_family")
    if _sequence_signal_is_supportive(selected_readout["sequence_support"]):
        required.append("sequence_signature")
    if rule.get("operators"):
        required.append("operator_readout")
    if selected_readout["trajectory_support"]:
        required.append("trajectory_readout")
    missing = [family for family in required if family not in bound]
    if missing:
        binding = "unbound_learned_mechanism"
    elif required:
        binding = "_".join(required) + "_bound"
    else:
        binding = "learned_mechanism_has_no_self_required_views"
    return {
        "cross_view_binding": binding,
        "bound_view_families": selected_readout["view_sources"],
        "self_required_view_families": required,
        "missing_view_families": missing,
        "sequence_support": selected_readout["sequence_support"],
        "operator_support": selected_readout["operator_support"],
        "trajectory_support": selected_readout["trajectory_support"],
    }


def _explicit_dominance_law(
    *,
    mechanism_class: str,
    selected_readout: dict[str, Any] | None,
    top_readout: dict[str, Any] | None,
    missing_word: str | None,
    cross_view_binding: dict[str, Any],
    operator_basis_probe: dict[str, Any],
    temporal_probe: dict[str, Any],
) -> str:
    if missing_word:
        return "candidate_missing_word_competes_with_learned_grammar"
    if mechanism_class == "insufficient_evidence_clean_abstain":
        return "no_dominant_learned_mechanism"
    if not selected_readout or selected_readout["grammar_status"] != "learned":
        return "no_learned_mechanism_readout"
    if not top_readout or selected_readout["grammar"] != top_readout["grammar"]:
        return "selected_mechanism_not_top_internal_readout"
    if cross_view_binding["cross_view_binding"] == "unbound_learned_mechanism":
        return "insufficient_cross_view_binding_for_acceptance"
    if operator_basis_probe["operator_basis_stability"] == "nondefining_operator_basis_sensitive":
        return "operator_basis_sensitive_mechanism"
    if temporal_probe["temporal_binding"] == "selected_observable_temporal_conflict":
        return "temporal_binding_conflict"
    return "single_dominant_learned_mechanism_bound_across_views"


def _internal_contradictions(
    *,
    mechanism: dict[str, Any],
    text: str,
    trajectory: dict[str, Any],
    readouts: list[dict[str, Any]],
    operator_basis_probe: dict[str, Any] | None = None,
    cross_view_binding: dict[str, Any] | None = None,
    temporal_probe: dict[str, Any] | None = None,
) -> list[str]:
    selected = mechanism["mechanism_class"]
    final = trajectory["final_state_summary"]
    contradictions: list[str] = []
    if _missing_word_candidate(readouts) and selected != "insufficient_evidence_clean_abstain":
        contradictions.append("learned_grammar_competes_with_candidate_missing_word")
    if operator_basis_probe and operator_basis_probe["operator_basis_stability"] == "nondefining_operator_basis_sensitive":
        contradictions.append("selected_mechanism_sensitive_to_endogenous_operator_coefficient_assignment")
    if cross_view_binding and cross_view_binding["cross_view_binding"] == "unbound_learned_mechanism":
        contradictions.append("selected_mechanism_lacks_evidence_operator_trajectory_binding")
    if temporal_probe and temporal_probe["temporal_binding"] == "selected_observable_temporal_conflict":
        contradictions.append("selected_mechanism_observable_decreases_across_simulated_timepoints")
    if selected == "globular_closure" and _contains_any(text, STRONG_MEMBRANE_CONTEXT_TOKENS):
        contradictions.append("soluble_globular_plus_strong_membrane_topology")
    if selected in {"intrinsic_disorder_phase_separation", "disorder_boundary_and_fold_upon_binding"} and float(final.get("contact_probability", 0.0)) > float(final.get("disorder_order_balance", 0.0)):
        contradictions.append("disorder_call_with_compact_core_readout")
    if selected == "assembly_required_folding" and _contains_any(text, ["complete soluble monomer", "standalone soluble fold", "soluble_monomeric_core_context"]):
        contradictions.append("assembly_required_plus_complete_monomer_context")
    if selected == "metal_cluster_and_ligand_locked_basin" and not _token_hits(text, METAL_CLUSTER_CONTEXT_TOKENS + LIGAND_LOCKED_CONTEXT_TOKENS):
        contradictions.append("ligand_locked_call_without_ligand_or_metal_evidence")
    if selected == "beta_closure_topology" and not _token_hits(text, BETA_CLOSURE_CONTEXT_TOKENS):
        contradictions.append("beta_closure_call_without_register_or_closure_evidence")
    if selected == "multidomain_allosteric_architecture" and not _token_hits(text, MULTIDOMAIN_ALLOSTERIC_CONTEXT_TOKENS):
        contradictions.append("multidomain_call_without_boundary_hinge_or_modular_signal")
    if selected == "secretory_disulfide_redox_topology":
        if not _token_hits(text, SECRETORY_DISULFIDE_CONTEXT_TOKENS):
            contradictions.append("secretory_disulfide_call_without_secretory_redox_evidence")
        if float(final.get("disulfide_pairing_topology", 0.0)) <= float(final.get("redox_mispaired_frustration", 0.0)):
            contradictions.append("secretory_disulfide_call_with_mispaired_redox_dominance")
    if selected == "signal_peptide_vs_true_tm_routing":
        if not _token_hits(text, SIGNAL_PEPTIDE_ROUTING_CONTEXT_TOKENS):
            contradictions.append("signal_peptide_routing_call_without_signal_context")
        if float(final.get("membrane_insertion_routing", 0.0)) >= float(final.get("secretory_lumenal_routing", 0.0)):
            contradictions.append("signal_peptide_routing_call_with_membrane_insertion_dominance")
        if float(final.get("signal_anchor_ambiguity", 0.0)) >= float(final.get("signal_peptide_routing_context", 0.0)):
            contradictions.append("signal_peptide_routing_call_with_signal_anchor_ambiguity_dominance")
    if selected == "coiled_coil_register_topology":
        if not _token_hits(text, COILED_COIL_REGISTER_CONTEXT_TOKENS):
            contradictions.append("coiled_coil_call_without_register_evidence")
        if float(final.get("register_shift_frustration", 0.0)) > max(
            float(final.get("heptad_register_context", 0.0)),
            float(final.get("oligomeric_coiled_coil_core", 0.0)),
        ):
            contradictions.append("coiled_coil_call_with_register_shift_frustration_dominance")
    if selected == "repeat_solenoid_topology":
        if not _token_hits(text, REPEAT_SOLENOID_CONTEXT_TOKENS):
            contradictions.append("repeat_solenoid_call_without_repeat_axis_evidence")
        if float(final.get("repeat_boundary_frustration", 0.0)) > float(final.get("global_repeat_topology", 0.0)):
            contradictions.append("repeat_solenoid_call_with_boundary_frustration_dominance")
    if selected == "knotted_topology":
        if not _token_hits(text, KNOTTED_TOPOLOGY_CONTEXT_TOKENS):
            contradictions.append("knotted_topology_call_without_explicit_topology_evidence")
        if float(final.get("unknotted_decoy_dominance", 0.0)) > float(final.get("topological_closure_constraint", 0.0)):
            contradictions.append("knotted_topology_call_with_unknotted_decoy_dominance")
    return contradictions


def _missing_word_candidate(readouts: list[dict[str, Any]]) -> str | None:
    learned_best = max(
        (row["internal_view_count"] for row in readouts if row["grammar_status"] == "learned"),
        default=0,
    )
    candidates = [
        row
        for row in readouts
        if row["grammar_status"] == "candidate_missing_word"
        and (
            row["evidence_hits"]
            or (
                "sequence_signature" in row["view_sources"]
                and "wrong_grammar_competition" in row["view_sources"]
                and row["internal_view_count"] > learned_best
            )
        )
    ]
    if not candidates:
        return None
    candidates = sorted(candidates, key=_rank_readout, reverse=True)
    return candidates[0]["grammar"]


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
        cys_his_metal_context = _contains_any(text, CYS_HIS_METAL_COORDINATION_TOKENS)
        metal_cluster_context = _contains_any(text, METAL_CLUSTER_CONTEXT_TOKENS) or cys_his_metal_context
        ligand_locked_context = _contains_any(text, LIGAND_LOCKED_CONTEXT_TOKENS)
        disorder_boundary_context = _contains_any(text, DISORDER_BOUNDARY_CONTEXT_TOKENS)
        beta_topology_word = _beta_topology_word_from_text(text)
        beta_closure_context = beta_topology_word is not None or _contains_any(text, BETA_CLOSURE_CONTEXT_TOKENS)
        beta_ambiguous_context = _contains_any(text, BETA_AMBIGUOUS_CONTEXT_TOKENS)
        multidomain_word = _multidomain_allosteric_word_from_text(text)
        multidomain_context = multidomain_word is not None or _contains_any(text, MULTIDOMAIN_ALLOSTERIC_CONTEXT_TOKENS)
        secretory_disulfide_context = _contains_any(text, SECRETORY_DISULFIDE_CONTEXT_TOKENS)
        signal_peptide_only_context = _contains_any(text, SIGNAL_PEPTIDE_ONLY_CONTEXT_TOKENS)
        signal_peptide_context = signal_peptide_only_context or _contains_any(text, SIGNAL_PEPTIDE_ROUTING_CONTEXT_TOKENS)
        explicit_cleavable_signal_context = _contains_any(
            text,
            [
                "signal_peptide_routing_context",
                "cleavage_site_context",
                "secretory_lumenal_routing",
                "cleaved signal peptide",
                "n-terminal signal peptide",
            ],
        )
        true_tm_routing_context = _contains_any(text, TRUE_TM_ROUTING_CONTEXT_TOKENS)
        signal_anchor_context = _contains_any(text, SIGNAL_ANCHOR_AMBIGUITY_TOKENS)
        coiled_coil_context = _contains_any(text, COILED_COIL_REGISTER_CONTEXT_TOKENS)
        repeat_solenoid_context = _contains_any(text, REPEAT_SOLENOID_CONTEXT_TOKENS)
        knotted_context = _contains_any(text, KNOTTED_TOPOLOGY_CONTEXT_TOKENS)
        signal_peptide_sequence_topology = _has_signal_peptide_routing_topology(sequence_field)
        paired_cysteine_topology = _has_paired_cysteine_topology(sequence_field)
        cysteine_pairing_conflict = _cysteine_topology_label(sequence_field) == "odd_cysteine_pairing_conflict"
        coiled_or_repeat_candidate_context = coiled_coil_context or repeat_solenoid_context
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
            natural = "membrane_multidomain_folding_proteostasis"
            reason = "strong_membrane_topology_context_prioritized_over_incidental_ligand_or_assembly_context"
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
        elif signal_peptide_context and signal_anchor_context:
            natural = "insufficient_evidence_clean_abstain"
            reason = "signal_anchor_ambiguity_requires_boundary_abstention"
        elif signal_peptide_context and true_tm_routing_context and not signal_peptide_only_context:
            natural = "insufficient_evidence_clean_abstain"
            reason = "true_tm_context_competes_with_signal_peptide_routing"
        elif (
            signal_peptide_context
            and signal_peptide_only_context
            and explicit_cleavable_signal_context
            and not true_tm_routing_context
            and not signal_anchor_context
        ):
            natural = "signal_peptide_vs_true_tm_routing"
            reason = "explicit_cleavable_signal_peptide_context_without_true_tm_or_anchor_competitor"
        elif (
            signal_peptide_context
            and (not secretory_disulfide_context or signal_peptide_only_context)
            and signal_peptide_sequence_topology
        ):
            natural = "signal_peptide_vs_true_tm_routing"
            reason = "cleavable_signal_peptide_routing_context_separated_from_true_tm"
        elif signal_peptide_context and not secretory_disulfide_context:
            natural = "insufficient_evidence_clean_abstain"
            reason = "signal_peptide_text_without_n_terminal_routing_topology_requires_abstention"
        elif secretory_disulfide_context and coiled_or_repeat_candidate_context:
            natural = "insufficient_evidence_clean_abstain"
            reason = "coiled_coil_or_repeat_candidate_competes_with_secretory_disulfide_text"
        elif secretory_disulfide_context and signal_peptide_only_context:
            natural = "insufficient_evidence_clean_abstain"
            reason = "signal_peptide_only_without_supported_disulfide_topology_requires_abstention"
        elif secretory_disulfide_context and cysteine_pairing_conflict:
            natural = "insufficient_evidence_clean_abstain"
            reason = "odd_or_conflicting_cysteine_pairing_topology_requires_abstention"
        elif secretory_disulfide_context and paired_cysteine_topology:
            natural = "secretory_disulfide_redox_topology"
            reason = "secretory_disulfide_redox_context_with_paired_cysteine_topology"
        elif secretory_disulfide_context:
            natural = "insufficient_evidence_clean_abstain"
            reason = "secretory_disulfide_text_without_cysteine_pairing_topology_requires_abstention"
        elif knotted_context:
            natural = "knotted_topology"
            reason = "explicit_knot_or_slipknot_threading_topology_context"
        elif repeat_solenoid_context:
            natural = "repeat_solenoid_topology"
            reason = "explicit_repeat_unit_phase_and_solenoid_axis_context"
        elif coiled_coil_context:
            natural = "coiled_coil_register_topology"
            reason = "explicit_heptad_register_coiled_coil_context"
        elif multidomain_context:
            natural = "multidomain_allosteric_architecture"
            reason = "explicit_multidomain_allosteric_architecture_context"
        elif cofactor_context:
            natural = "cofactor_ligand_assisted_stabilization"
            reason = "explicit_ligand_cofactor_or_metal_context"
        elif soluble_monomeric_core_context and _not_dominated_by(metrics["hydrophobic_density"], metrics["mean_disorder"]):
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
                or _dominates(metrics["mean_membrane"], metrics["mean_disorder"])
            ):
                natural = "membrane_multidomain_folding_proteostasis"
                reason = "membrane_multidomain_mutation_interface_evidence"
            else:
                natural = "insufficient_evidence_clean_abstain"
                reason = "generic_membrane_annotation_without_specific_operator"
        elif any(token in text for token in ["low complexity", "disprot", "disordered", "phase", "prion", "lcd", "fus", "tdp-43", "tdp43"]):
            if _dominates(_endogenous_support(metrics["low_complexity_density"], metrics["mean_disorder"]), metrics["hydrophobic_density"]):
                natural = "intrinsic_disorder_phase_separation"
                reason = "low_complexity_disorder_phase_evidence"
            else:
                natural = "insufficient_evidence_clean_abstain"
                reason = "disorder_text_without_sequence_support"
        elif any(token in text for token in ["binding", "partner", "motif"]) and _dominates(metrics["mean_disorder"], metrics["hydrophobic_density"]):
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
        ]) and _not_dominated_by(metrics["hydrophobic_density"], metrics["mean_disorder"]):
            natural = "globular_closure"
            reason = "noncoordinate_process_metadata_and_sequence_support_globular_closure"
        elif any(token in text for token in ["mini-protein", "miniprotein", "trp-cage", "caged aromatic core"]):
            natural = "globular_closure"
            reason = "small_fast_folder_metadata_supports_globular_closure"
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
        "selected_secretory_disulfide_word": "disulfide_secretory_redox_context"
        if mechanism_class == "secretory_disulfide_redox_topology"
        else None,
        "selected_signal_peptide_word": "signal_peptide_vs_true_TM"
        if mechanism_class == "signal_peptide_vs_true_tm_routing"
        else None,
        "selected_coiled_coil_word": "coiled_coil_register"
        if mechanism_class == "coiled_coil_register_topology"
        else None,
        "selected_repeat_solenoid_word": "repeat_solenoid_topology"
        if mechanism_class == "repeat_solenoid_topology"
        else None,
        "selected_knotted_topology_word": "knotted_topology"
        if mechanism_class == "knotted_topology"
        else None,
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


def _endogenous_mean(*values: float) -> float:
    positives = [float(value) for value in values if float(value) > 0.0]
    return bounded(sum(positives) / len(positives)) if positives else 0.0


def _endogenous_support(*values: float) -> float:
    return _endogenous_mean(*values)


def _attenuate_by_competing_pressure(value: float, *pressures: float) -> float:
    pressure = _endogenous_mean(*pressures)
    value = bounded(value)
    total = value + pressure
    if not _present(total):
        return 0.0
    return bounded(value * (value / total))


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
    cysteine_segments = _strong_segments(sequence_field, "cysteine_density")
    n_terminal_segments = sequence_field["segments"][:2]
    if mechanism_class == "globular_closure":
        operators.append(_operator(
            "closure_operator",
            ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
            _endogenous_support(metrics["hydrophobic_density"], metrics["aromatic_density"], metrics["mean_interface"]),
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
                _endogenous_support(metrics["mean_disorder"], metrics["low_complexity_density"]),
                evidence,
                "expanded dynamic ensemble with no stable single fold",
                "charge/proline/glycine changes tune expansion",
                "stable globular fold dominates without condition dependence",
                "disorder_order_balance",
            ),
            _operator(
                "phase_operator",
                ", ".join(_span_from_segment(row) for row in phase_segments),
                _endogenous_support(metrics["low_complexity_density"], metrics["aromatic_density"]),
                evidence,
                "weak multivalent attraction and phase-prone condensate threshold",
                "aromatic/charge/salt/RNA perturbations shift phase threshold",
                "phase behavior unaffected by predicted stickers or conditions",
                "state_basin_occupancy",
            ),
            _operator(
                "repulsion_operator",
                "charged low-complexity windows",
                _endogenous_support(abs(metrics["net_charge_per_residue"]), metrics["low_complexity_density"]),
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
                _endogenous_support(metrics["mean_disorder"], metrics["low_complexity_density"]),
                evidence,
                "persistent IDR boundary and structured-domain/tail separation",
                "boundary truncation or charge/proline/glycine edits shift disorder persistence",
                "generic oligomer interface explains the holdout without an IDR boundary",
                "IDR_boundary",
            ),
            _operator(
                "interface_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments),
                _endogenous_support(metrics["mean_interface"], metrics["mean_disorder"]),
                evidence,
                "local fold-upon-binding motif readiness without whole-chain compaction",
                "partner or motif removal weakens local ordering",
                "partner context has no directional effect on the predicted local motif",
                "fold_upon_binding_region",
            ),
            _operator(
                "phase_operator",
                ", ".join(_span_from_segment(row) for row in phase_segments),
                _endogenous_support(metrics["low_complexity_density"], metrics["aromatic_density"]),
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
                _endogenous_support(metrics["hydrophobic_density"], metrics["beta_propensity_density"]),
                evidence,
                f"{beta_word} strand-register closure rather than generic beta propensity",
                "register-shift or blade/repeat perturbation weakens closed beta topology",
                "open beta sheet or wrong beta subtype explains the holdout better",
                beta_word,
            ),
            _operator(
                "interface_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments),
                _endogenous_support(metrics["mean_interface"], metrics["beta_propensity_density"]),
                evidence,
                "inter-strand edge pairing and closure interface readiness",
                "edge/interface perturbation opens the beta topology",
                "beta-rich region stays open without closure contacts",
                "beta_sheet_closure",
            ),
            _operator(
                "frustration_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                _endogenous_support(metrics["aromatic_density"], metrics["mean_interface"]),
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
                _endogenous_support(metrics["mean_interface"], metrics["hydrophobic_density"]),
                evidence,
                "interdomain lock and domain-boundary interface readiness",
                "interface or lock perturbation weakens interdomain coupling",
                "single-domain closure explains the holdout without interdomain dependence",
                "interdomain_lock",
            ),
            _operator(
                "dual_basin_switch_operator",
                "hinge/allosteric domain-coupling axis",
                _endogenous_support(metrics["mean_interface"], metrics["aromatic_density"]),
                evidence,
                "domain reorientation shifts an allosteric basin rather than one averaged fold",
                "hinge or allosteric push changes basin occupancy and orientation",
                "domain orientation is condition-independent",
                "allosteric_basin_shift",
            ),
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                _endogenous_support(metrics["hydrophobic_density"], metrics["mean_interface"]),
                evidence,
                "modular architecture closes through coupled domain packing",
                "domain-boundary or packing damage lowers modular compaction",
                "a standalone compact domain explains the full architecture",
                "modular_architecture",
            ),
            _operator(
                "frustration_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments[:2]),
                _endogenous_support(metrics["mean_interface"], metrics["aromatic_density"]),
                evidence,
                f"{multidomain_word} separates hinge, lock, reorientation, and swapped-domain topology from generic domains",
                "wrong domain-boundary operator raises allosteric or swapping conflict",
                "generic domain text alone predicts all states equally",
                multidomain_word,
            ),
        ])
    elif mechanism_class == "secretory_disulfide_redox_topology":
        cysteine_span = ", ".join(_span_from_segment(row) for row in cysteine_segments)
        interface_span = ", ".join(_span_from_segment(row) for row in interface_segments)
        operators.extend([
            _operator(
                "disulfide_pairing_operator",
                cysteine_span,
                _endogenous_support(metrics["cysteine_density"], metrics["mean_interface"]),
                evidence,
                "paired cysteine topology closes as a secretory disulfide network",
                "cysteine-pairing damage weakens disulfide topology and extracellular stabilization",
                "ordinary cysteine noise or metal coordination explains the evidence better",
                "disulfide_pairing_topology",
            ),
            _operator(
                "secretory_redox_operator",
                cysteine_span,
                _endogenous_support(metrics["cysteine_density"], metrics["n_terminal_hydrophobic_signal"]),
                evidence,
                "secretory redox environment supports disulfide maturation after signal peptide removal",
                "redox perturbation or signal-context masking lowers secretory redox support",
                "signal peptide text alone explains the evidence without redox topology",
                "secretory_redox_context",
            ),
            _operator(
                "frustration_operator",
                cysteine_span,
                _endogenous_support(metrics["cysteine_density"], metrics["histidine_density"]),
                evidence,
                "mispaired redox frustration separates disulfide topology from Cys-His metal grammar",
                "mispaired cysteine or wrong-redox perturbation raises frustration",
                "all cysteine-rich proteins behave as the same globular fold",
                "redox_mispaired_frustration",
            ),
            _operator(
                "proteostasis_operator",
                "secretory quality-control route",
                _endogenous_support(metrics["cysteine_density"], metrics["mean_interface"]),
                evidence,
                "secretory quality control distinguishes extracellular maturation from cytosolic cysteine noise",
                "quality-control stress weakens extracellular stabilized fold",
                "secretory route has no directional effect",
                "secretory_quality_control",
            ),
            _operator(
                "closure_operator",
                interface_span,
                _endogenous_support(metrics["hydrophobic_density"], metrics["cysteine_density"]),
                evidence,
                "extracellular stabilized fold closes through disulfide-constrained topology",
                "disulfide reduction lowers closure more than matched controls",
                "ordinary globular closure explains the target without disulfide/redox state",
                "extracellular_stabilized_fold",
            ),
        ])
    elif mechanism_class == "signal_peptide_vs_true_tm_routing":
        n_terminal_span = ", ".join(_span_from_segment(row) for row in n_terminal_segments)
        membrane_span = ", ".join(_span_from_segment(row) for row in membrane_segments)
        operators.extend([
            _operator(
                "signal_peptide_routing_operator",
                n_terminal_span,
                _endogenous_support(metrics["n_terminal_hydrophobic_signal"], metrics["mean_interface"]),
                evidence,
                "cleavable N-terminal signal peptide routes to secretory entry",
                "N-terminal hydrophobic masking weakens signal-peptide routing",
                "internal true TM span explains the hydrophobic signal better",
                "signal_peptide_routing_context",
            ),
            _operator(
                "cleavage_context_operator",
                n_terminal_span,
                _endogenous_support(metrics["n_terminal_hydrophobic_signal"], abs(metrics["net_charge_per_residue"])),
                evidence,
                "cleavage-site context separates signal peptide from signal anchor",
                "cleavage-context masking weakens lumenal routing",
                "uncleaved signal-anchor decoy explains the route",
                "cleavage_site_context",
            ),
            _operator(
                "secretory_routing_operator",
                n_terminal_span,
                _endogenous_support(metrics["n_terminal_hydrophobic_signal"], metrics["mean_interface"]),
                evidence,
                "secretory lumenal routing follows cleavable signal context",
                "signal peptide perturbation lowers secretory lumenal route",
                "membrane insertion route dominates instead",
                "secretory_lumenal_routing",
            ),
            _operator(
                "tm_insertion_operator",
                membrane_span,
                _endogenous_support(metrics["mean_membrane"], metrics["hydrophobic_density"]),
                evidence,
                "true TM insertion remains a competing boundary readout",
                "TM decoy context must not validate signal-peptide routing",
                "signal peptide grammar steals multi-pass membrane topology",
                "true_transmembrane_span_context",
            ),
            _operator(
                "membrane_pressure_operator",
                membrane_span,
                _endogenous_support(metrics["mean_membrane"], metrics["n_terminal_hydrophobic_signal"]),
                evidence,
                "membrane pressure is reported separately from cleaved secretory routing",
                "membrane pressure increase shifts toward insertion/anchor ambiguity",
                "membrane pressure has no boundary effect",
                "membrane_insertion_routing",
            ),
            _operator(
                "frustration_operator",
                n_terminal_span,
                _endogenous_support(metrics["mean_membrane"], metrics["n_terminal_hydrophobic_signal"]),
                evidence,
                "signal-anchor ambiguity is exposed instead of accepted",
                "signal-anchor decoy raises boundary ambiguity",
                "signal anchor and cleaved signal peptide are indistinguishable",
                "signal_anchor_ambiguity",
            ),
        ])
    elif mechanism_class == "coiled_coil_register_topology":
        heptad_span = ", ".join(_span_from_segment(row) for row in hydrophobic_segments)
        interface_span = ", ".join(_span_from_segment(row) for row in interface_segments)
        operators.extend([
            _operator(
                "heptad_register_operator",
                heptad_span,
                _endogenous_support(
                    metrics.get("heptad_hydrophobic_periodicity", 0.0),
                    metrics["hydrophobic_density"],
                    metrics["mean_interface"],
                ),
                evidence,
                "heptad phase locks the coiled-coil register instead of generic helix closure",
                "heptad shuffle or register shift lowers register alignment",
                "ordinary helix bundle or assembly interface explains the holdout better",
                "heptad_register_context",
            ),
            _operator(
                "coiled_coil_interface_operator",
                interface_span or heptad_span,
                _endogenous_support(metrics["mean_interface"], metrics["hydrophobic_density"]),
                evidence,
                "hydrophobic seam packs as an oligomeric coiled-coil interface",
                "hydrophobic phase shift weakens the coiled-coil interface",
                "globular hydrophobic core remains equally coherent after phase shift",
                "hydrophobic_repeat_phase",
            ),
            _operator(
                "oligomeric_register_operator",
                interface_span or heptad_span,
                _endogenous_support(metrics["mean_interface"], metrics.get("heptad_hydrophobic_periodicity", 0.0)),
                evidence,
                "parallel/antiparallel register and partner-copy core are coupled",
                "partner-register perturbation weakens oligomeric core alignment",
                "generic assembly without heptad phase explains the evidence",
                "oligomeric_coiled_coil_core",
            ),
            _operator(
                "register_shift_frustration_operator",
                heptad_span,
                _endogenous_support(metrics["aromatic_density"], metrics["mean_interface"]),
                evidence,
                "register-shift frustration is exposed as a falsifier",
                "forcing wrong register raises frustration",
                "register shifts have no directional effect",
                "register_shift_frustration",
            ),
            _operator(
                "closure_operator",
                heptad_span,
                _endogenous_support(metrics["hydrophobic_density"], metrics["mean_interface"]),
                evidence,
                "local closure is subordinated to heptad register",
                "register damage lowers closure support",
                "generic closure beats registered coiled-coil state",
                "coiled_coil_assembly_dependency",
            ),
        ])
    elif mechanism_class == "repeat_solenoid_topology":
        repeat_segments = sorted(
            sequence_field["segments"],
            key=lambda row: row["interface_density"] + row["beta_propensity_density"] + row["pro_gly_density"],
            reverse=True,
        )[:3]
        repeat_span = ", ".join(_span_from_segment(row) for row in repeat_segments or interface_segments)
        operators.extend([
            _operator(
                "repeat_phase_operator",
                repeat_span,
                _endogenous_support(metrics.get("repeat_signature", 0.0), metrics["mean_interface"]),
                evidence,
                "repeat units align in phase rather than as unrelated domains",
                "repeat-order shuffle lowers phase alignment",
                "multidomain or beta closure grammar explains the target better",
                "repeat_phase_alignment",
            ),
            _operator(
                "solenoid_axis_operator",
                repeat_span,
                _endogenous_support(metrics.get("repeat_signature", 0.0), metrics["hydrophobic_density"]),
                evidence,
                "local repeat units curve around a coherent solenoid axis",
                "axis or boundary masking weakens solenoid continuity",
                "generic long protein context has no repeat-axis operator",
                "solenoid_axis_context",
            ),
            _operator(
                "local_repeat_closure_operator",
                repeat_span,
                _endogenous_support(metrics["hydrophobic_density"], metrics["beta_propensity_density"]),
                evidence,
                "local repeat closure propagates through adjacent units",
                "repeat boundary damage lowers local closure",
                "single closed beta or globular unit explains all contacts",
                "local_repeat_closure",
            ),
            _operator(
                "global_repeat_stack_operator",
                repeat_span,
                _endogenous_support(metrics.get("repeat_signature", 0.0), metrics["mean_interface"]),
                evidence,
                "global topology is a phase-aligned repeat stack",
                "repeat-order shuffle damages global repeat topology",
                "a generic multidomain chain preserves the same topology",
                "global_repeat_topology",
            ),
            _operator(
                "repeat_boundary_frustration_operator",
                repeat_span,
                _endogenous_support(metrics["aromatic_density"], metrics["mean_interface"]),
                evidence,
                "boundary frustration distinguishes true solenoid axis from accidental repeats",
                "boundary masking raises repeat conflict",
                "boundary edits are inert",
                "repeat_boundary_frustration",
            ),
        ])
    elif mechanism_class == "knotted_topology":
        interface_span = ", ".join(_span_from_segment(row) for row in interface_segments)
        hydrophobic_span = ", ".join(_span_from_segment(row) for row in hydrophobic_segments)
        operators.extend([
            _operator(
                "threading_operator",
                interface_span or hydrophobic_span,
                _endogenous_support(metrics["mean_interface"], metrics["hydrophobic_density"]),
                evidence,
                "threading loop forms a topology-constrained path rather than generic closure",
                "threading damage lowers knot-core and loop readouts",
                "unknotted globular closure explains the target better",
                "threading_loop_context",
            ),
            _operator(
                "topological_closure_operator",
                hydrophobic_span or interface_span,
                _endogenous_support(metrics["hydrophobic_density"], metrics["mean_interface"]),
                evidence,
                "topological closure preserves knot core after threading",
                "topology mask lowers closure constraint",
                "closed beta/repeat grammar captures the same topology",
                "topological_closure_constraint",
            ),
            _operator(
                "long_range_threading_operator",
                f"{_span_from_segment(sequence_field['segments'][0])}-{_span_from_segment(sequence_field['segments'][-1])}",
                _endogenous_support(metrics["mean_interface"], metrics["aromatic_density"]),
                evidence,
                "long-range order is required for the knot/slipknot route",
                "long-range shuffle weakens threading dependency",
                "local closure alone predicts the state",
                "long_range_threading_dependency",
            ),
            _operator(
                "slipknot_intermediate_operator",
                interface_span or hydrophobic_span,
                _endogenous_support(metrics["mean_interface"], metrics["hydrophobic_density"]),
                evidence,
                "slipknot intermediate is retained as a mechanistic waypoint",
                "slipknot masking lowers the intermediate readout",
                "no intermediate is needed for the final topology",
                "slipknot_intermediate_context",
            ),
            _operator(
                "knotting_frustration_operator",
                hydrophobic_span or interface_span,
                _endogenous_support(metrics["aromatic_density"], metrics["mean_interface"]),
                evidence,
                "unknotted decoy and knotting frustration are explicit falsifiers",
                "wrong topology raises unknotted decoy dominance",
                "topology perturbations are inert",
                "knotting_frustration",
            ),
        ])
    elif mechanism_class == "membrane_multidomain_folding_proteostasis":
        f508 = next((row for row in focus_regions if "F508" in str(row.get("name", "")) or row.get("position") == 508), None)
        f508_span = "508" if f508 is None else str(f508.get("span", f508.get("position", "508")))
        operators.extend([
            _operator(
                "membrane_pressure_operator",
                ", ".join(_span_from_segment(row) for row in membrane_segments),
                _endogenous_support(metrics["mean_membrane"], metrics["hydrophobic_density"]),
                evidence,
                "membrane-buried routing pressure and topology context",
                "membrane/context disruption weakens maturation route",
                "membrane context irrelevant to holdout rescue logic",
                "proteostasis_routing",
            ),
            _operator(
                "closure_operator",
                f"NBD1/local deletion focus {f508_span}",
                _endogenous_support(metrics["hydrophobic_density"], metrics["mean_interface"]),
                evidence,
                "NBD1 local stability and partial domain closure",
                "F508del weakens; NBD1-only correction partially rescues",
                "NBD1 stability evidence fails to track F508del",
                "segment_compaction",
            ),
            _operator(
                "interface_operator",
                "NBD1-MSD/interdomain correction axis",
                _endogenous_support(metrics["mean_interface"], metrics["mean_membrane"]),
                evidence,
                "interdomain interface readiness",
                "interface correction strengthens rescue beyond NBD1-only correction",
                "interface/proteostasis correction has no added effect",
                "interface_readiness",
            ),
            _operator(
                "proteostasis_operator",
                "folding quality-control and trafficking route",
                _endogenous_support(metrics["mean_membrane"], metrics["mean_interface"]),
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
                _endogenous_support(metrics["hydrophobic_density"], metrics["aromatic_density"], metrics["mean_interface"]),
                evidence,
                "two incompatible alpha/beta state basins",
                "partner/release perturbations shift alpha/beta occupancy",
                "only one stable state under all contexts",
                "state_basin_occupancy",
            ),
            _operator(
                "interface_operator",
                "partner-context interface",
                _endogenous_support(metrics["mean_interface"], metrics["aromatic_density"]),
                evidence,
                "context-conditioned state selection",
                "partner removal or stabilization shifts state occupancy",
                "partner context has no directional effect",
                "interface_readiness",
            ),
            _operator(
                "frustration_operator",
                "secondary-structure conflict region",
                _endogenous_support(metrics["aromatic_density"], metrics["mean_interface"]),
                evidence,
                "state conflict rather than averaged consensus fold",
                "forcing one fold fails",
                "one averaged state explains all holdouts",
                "frustration",
            ),
        ])
    elif mechanism_class == "short_region_host_interface_hijacking":
        cterm = _span_from_segment(sequence_field["segments"][-1])
        operators.extend([
            _operator(
                "host_hijack_operator",
                f"C-terminal host-interface region {cterm}",
                _endogenous_support(metrics["mean_interface"], metrics["mean_disorder"], metrics["aromatic_density"]),
                evidence,
                "host-interface capture without global fold requirement",
                "C-terminal disruption weakens host binding and transport/IFN consequences",
                "wrong-region perturbations dominate while C-terminus is inert",
                "interface_readiness",
            ),
            _operator(
                "interface_operator",
                f"short linear motif/interface window {cterm}",
                _endogenous_support(metrics["mean_interface"], metrics["aromatic_density"]),
                evidence,
                "RAE1/NUP98 or host-surface readiness",
                "host partner removal weakens the operator",
                "host partner evidence does not localize to predicted region",
                "contact_probability",
            ),
            _operator(
                "disorder_operator",
                "short exposed accessory protein context",
                _endogenous_support(metrics["mean_disorder"], metrics["low_complexity_density"]),
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
                _endogenous_support(metrics["mean_interface"], metrics["aromatic_density"]),
                evidence,
                "ligand/cofactor pocket readiness and local stabilization",
                "cofactor removal or pocket disruption weakens the stabilized basin",
                "apo and ligand-bound contexts are indistinguishable",
                "interface_readiness",
            ),
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                _endogenous_support(metrics["hydrophobic_density"], metrics["mean_interface"]),
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
                _endogenous_support(metrics["mean_interface"], metrics["aromatic_density"], metrics["histidine_density"]),
                evidence,
                "coordination shell and ligand pocket complete the locked basin",
                "metal/cofactor removal or coordinating-side-chain disruption unlocks the basin",
                "generic cofactor annotation explains the holdout without a geometry-locked basin",
                "coordination_shell_integrity",
            ),
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                _endogenous_support(metrics["hydrophobic_density"], metrics["mean_interface"]),
                evidence,
                "ligand-locked compaction rather than free apo closure",
                "apo conversion lowers basin occupancy and compaction",
                "ligand-free closure remains equally stable",
                "ligand_locked_basin",
            ),
            _operator(
                "frustration_operator",
                ", ".join(_span_from_segment(row) for row in interface_segments[:2]),
                _endogenous_support(metrics["mean_interface"], metrics["aromatic_density"]),
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
                _endogenous_support(metrics["mean_interface"], metrics["hydrophobic_density"]),
                evidence,
                "partner-completed core and biological assembly readiness",
                "partner/interface disruption exposes incomplete monomer topology",
                "complete monomer, membrane, or ligand grammar explains the holdout without assembly",
                "partner_completed_core",
            ),
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                _endogenous_support(metrics["hydrophobic_density"], metrics["mean_interface"]),
                evidence,
                "hydrophobic core closure only after assembly context resolves the surface",
                "assembly-interface mutation lowers partner-completed closure",
                "monomer-only closure remains stable and complete",
                "assembly_required_core",
            ),
            _operator(
                "frustration_operator",
                ", ".join(_span_from_segment(row) for row in membrane_segments),
                _endogenous_support(metrics["mean_membrane"], metrics["mean_interface"]),
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
                _endogenous_support(metrics["mean_interface"], metrics["hydrophobic_density"]),
                evidence,
                "partner-copy interface readiness and assembly-stabilized folding",
                "interface or concentration perturbation weakens the assembled basin",
                "monomer-only grammar explains the holdout equally well",
                "interface_readiness",
            ),
            _operator(
                "closure_operator",
                ", ".join(_span_from_segment(row) for row in hydrophobic_segments),
                _endogenous_support(metrics["hydrophobic_density"], metrics["mean_interface"]),
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
            _endogenous_support(metrics["mean_interface"], metrics["hydrophobic_density"]),
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
    return 0.0


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
    trajectory_steps = list(range(4))
    for step_index, timepoint in enumerate(trajectory_steps):
        progress = step_index / max(1, len(trajectory_steps) - 1)
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
            expanded = _attenuate_by_competing_pressure(_endogenous_support(disorder, phase), closure * progress)
            phase_dynamic = _endogenous_support(phase * progress, disorder)
            compact = _attenuate_by_competing_pressure(closure * progress, disorder, phase)
            basin = {
                "expanded_disordered": expanded,
                "phase_prone_dynamic": phase_dynamic,
                "compact_single_fold": compact,
            }
            segment_compaction = _endogenous_support(phase * progress, closure * progress, noise)
            contact_probability = _endogenous_support(phase * progress, noise)
            interface_readiness = _endogenous_support(interface * progress, phase_dynamic)
            proteostasis_routing = 0.0
        elif mechanism_class == "disorder_boundary_and_fold_upon_binding":
            partner_loss = float((perturbation or {}).get("partner_loss", 0.0))
            motif_damage = float((perturbation or {}).get("motif_damage", 0.0))
            idr_boundary = _attenuate_by_competing_pressure(_endogenous_support(disorder, phase), motif_damage)
            local_order = _attenuate_by_competing_pressure(
                _endogenous_support(interface * progress, disorder),
                partner_loss,
                motif_damage,
            )
            low_complexity_basin = _endogenous_support(phase * progress, disorder)
            basin = {
                "disorder_boundary_ensemble": idr_boundary,
                "fold_upon_binding_basin": local_order,
                "phase_prone_low_complexity": low_complexity_basin,
                "compact_single_fold": _attenuate_by_competing_pressure(closure * progress, disorder, phase),
            }
            segment_compaction = _endogenous_support(local_order, closure * progress, noise)
            contact_probability = _endogenous_support(local_order, phase * progress)
            interface_readiness = local_order
            proteostasis_routing = 0.0
        elif mechanism_class == "beta_closure_topology":
            register_damage = float((perturbation or {}).get("register_damage", 0.0))
            closure_conflict = float((perturbation or {}).get("closure_conflict", 0.0))
            register = _attenuate_by_competing_pressure(_endogenous_support(closure * progress, interface), register_damage)
            closure_state = _attenuate_by_competing_pressure(
                _endogenous_support(closure * progress, interface, metrics["beta_propensity_density"]),
                closure_conflict,
            )
            conflict = _attenuate_by_competing_pressure(
                _endogenous_support(frustration_strength, closure_conflict, register_damage),
                closure_state,
            )
            basin = {
                "closed_beta_topology": closure_state,
                "strand_register": register,
                "open_beta_sheet_ambiguous": _attenuate_by_competing_pressure(conflict, closure_state),
                "wrong_beta_topology": conflict,
            }
            segment_compaction = _endogenous_support(closure_state, closure * progress, noise)
            contact_probability = _endogenous_support(closure_state, register)
            interface_readiness = _endogenous_support(interface * progress, register)
            proteostasis_routing = 0.0
        elif mechanism_class == "multidomain_allosteric_architecture":
            hinge_damage = float((perturbation or {}).get("hinge_damage", 0.0))
            lock_damage = float((perturbation or {}).get("lock_damage", 0.0))
            allosteric_push = float((perturbation or {}).get("allosteric_push", 0.0))
            domain_boundary_state = _attenuate_by_competing_pressure(
                _endogenous_support(interface, frustration_strength, metrics["mean_interface"], progress),
                hinge_damage,
            )
            interdomain_lock = _attenuate_by_competing_pressure(
                _endogenous_support(interface * progress, closure, domain_boundary_state),
                lock_damage,
            )
            allosteric_shift = _attenuate_by_competing_pressure(
                _endogenous_support(switch * progress, frustration_strength, allosteric_push),
                lock_damage,
            )
            reorientation = _attenuate_by_competing_pressure(
                _endogenous_support(switch * progress, interface, allosteric_shift),
                hinge_damage,
            )
            hinge = _attenuate_by_competing_pressure(
                _endogenous_support(frustration_strength, domain_boundary_state, switch),
                hinge_damage,
            )
            modular = _endogenous_support(closure * progress, interface, domain_boundary_state)
            swapped = interdomain_lock if multidomain_subtype == "domain_swapping" else _attenuate_by_competing_pressure(interdomain_lock, modular)
            basin = {
                "multidomain_allosteric_basin": bounded(max(allosteric_shift, reorientation, interdomain_lock)),
                "domain_boundary": domain_boundary_state,
                "hinge_region": hinge,
                "interdomain_locked_basin": interdomain_lock,
                "allosteric_basin_shift": allosteric_shift,
                "domain_reorientation_basin": reorientation,
                "modular_architecture": modular,
                "domain_swapping": swapped,
                "single_domain_shortcut": _attenuate_by_competing_pressure(lock_damage, modular, interdomain_lock),
            }
            segment_compaction = _endogenous_support(modular, closure * progress, noise)
            contact_probability = _endogenous_support(interdomain_lock, allosteric_shift, closure * progress)
            interface_readiness = interdomain_lock
            proteostasis_routing = 0.0
        elif mechanism_class == "secretory_disulfide_redox_topology":
            disulfide_damage = float((perturbation or {}).get("disulfide_damage", 0.0))
            redox_shift = float((perturbation or {}).get("redox_shift", 0.0))
            glycosylation_loss = float((perturbation or {}).get("glycosylation_loss", 0.0))
            quality_control_stress = float((perturbation or {}).get("quality_control_stress", 0.0))
            disulfide_pairing = strengths["disulfide_pairing_operator"]
            secretory_redox = strengths["secretory_redox_operator"]
            pairing = _attenuate_by_competing_pressure(
                _endogenous_support(disulfide_pairing * progress, secretory_redox, closure),
                disulfide_damage,
                redox_shift,
            )
            redox_context = _attenuate_by_competing_pressure(
                _endogenous_support(secretory_redox * progress, proteostasis),
                redox_shift,
                quality_control_stress,
            )
            extracellular = _attenuate_by_competing_pressure(
                _endogenous_support(pairing, redox_context, closure * progress),
                glycosylation_loss,
            )
            glyco = _attenuate_by_competing_pressure(
                _endogenous_support(secretory_redox, proteostasis, metrics["cysteine_density"]),
                glycosylation_loss,
            )
            quality = _attenuate_by_competing_pressure(
                _endogenous_support(proteostasis * progress, redox_context),
                quality_control_stress,
            )
            mispaired = _attenuate_by_competing_pressure(
                _endogenous_support(frustration_strength, disulfide_damage, redox_shift),
                pairing,
            )
            signal_removed = _endogenous_support(
                secretory_redox,
                metrics["n_terminal_hydrophobic_signal"],
            )
            basin = {
                "secretory_redox_context": redox_context,
                "disulfide_pairing_topology": pairing,
                "extracellular_stabilized_fold": extracellular,
                "glycosylation_context": glyco,
                "secretory_quality_control": quality,
                "redox_mispaired_frustration": mispaired,
                "signal_peptide_removed_context": signal_removed,
                "ordinary_cysteine_globular_noise": _attenuate_by_competing_pressure(disulfide_damage, pairing),
            }
            segment_compaction = _endogenous_support(extracellular, closure * progress, noise)
            contact_probability = _endogenous_support(pairing, extracellular)
            interface_readiness = _endogenous_support(pairing, redox_context)
            proteostasis_routing = quality
        elif mechanism_class == "signal_peptide_vs_true_tm_routing":
            n_terminal_mask = float((perturbation or {}).get("n_terminal_mask", 0.0))
            cleavage_loss = float((perturbation or {}).get("cleavage_loss", 0.0))
            true_tm_decoy = float((perturbation or {}).get("true_tm_decoy", 0.0))
            signal_anchor_decoy = float((perturbation or {}).get("signal_anchor_decoy", 0.0))
            signal_route = strengths["signal_peptide_routing_operator"]
            cleavage = strengths["cleavage_context_operator"]
            secretory_route = strengths["secretory_routing_operator"]
            tm_insert = strengths["tm_insertion_operator"]
            n_patch = _attenuate_by_competing_pressure(
                _endogenous_support(metrics["n_terminal_hydrophobic_signal"], signal_route),
                n_terminal_mask,
            )
            cleavage_state = _attenuate_by_competing_pressure(
                _endogenous_support(cleavage * progress, signal_route),
                cleavage_loss,
                signal_anchor_decoy,
            )
            signal_context = _attenuate_by_competing_pressure(
                _endogenous_support(signal_route * progress, cleavage_state, secretory_route),
                n_terminal_mask,
                true_tm_decoy,
                signal_anchor_decoy,
            )
            tm_context = _attenuate_by_competing_pressure(
                _endogenous_support(tm_insert * progress, membrane, true_tm_decoy, signal_anchor_decoy),
                signal_context,
            )
            secretory_lumenal = _attenuate_by_competing_pressure(
                _endogenous_support(signal_context, cleavage_state, secretory_route * progress),
                true_tm_decoy,
                signal_anchor_decoy,
            )
            membrane_insert = _attenuate_by_competing_pressure(
                _endogenous_support(tm_context, membrane, true_tm_decoy, signal_anchor_decoy),
                secretory_lumenal,
            )
            single_conflict = _endogenous_support(max(0.0, tm_context - signal_context), signal_anchor_decoy)
            multi_conflict = _attenuate_by_competing_pressure(
                _endogenous_support(true_tm_decoy, membrane),
                signal_context,
            )
            anchor_ambiguity = _endogenous_support(
                min(signal_context, tm_context),
                signal_anchor_decoy,
                frustration_strength,
            )
            basin = {
                "signal_peptide_routing_context": signal_context,
                "cleavage_site_context": cleavage_state,
                "n_terminal_secretory_hydrophobic_patch": n_patch,
                "true_transmembrane_span_context": tm_context,
                "single_pass_tm_conflict": single_conflict,
                "multi_pass_tm_conflict": multi_conflict,
                "secretory_lumenal_routing": secretory_lumenal,
                "membrane_insertion_routing": membrane_insert,
                "signal_anchor_ambiguity": anchor_ambiguity,
            }
            segment_compaction = _endogenous_support(secretory_lumenal, membrane_insert, closure * progress, noise)
            contact_probability = _endogenous_support(signal_context, secretory_lumenal, membrane_insert)
            interface_readiness = _endogenous_support(secretory_lumenal, signal_context)
            proteostasis_routing = _endogenous_support(secretory_lumenal, membrane_insert)
        elif mechanism_class == "coiled_coil_register_topology":
            register_shift = float((perturbation or {}).get("register_shift", 0.0))
            heptad_shuffle = float((perturbation or {}).get("heptad_shuffle", 0.0))
            phase_shift = float((perturbation or {}).get("phase_shift", 0.0))
            register_operator = strengths["heptad_register_operator"]
            coil_interface = strengths["coiled_coil_interface_operator"]
            oligomeric_register = strengths["oligomeric_register_operator"]
            register_frustration = strengths["register_shift_frustration_operator"]
            heptad_state = _attenuate_by_competing_pressure(
                _endogenous_support(register_operator * progress, metrics.get("heptad_hydrophobic_periodicity", 0.0)),
                heptad_shuffle,
                register_shift,
            )
            phase_state = _attenuate_by_competing_pressure(
                _endogenous_support(coil_interface * progress, heptad_state, metrics["hydrophobic_density"]),
                phase_shift,
                heptad_shuffle,
            )
            register_pairing = _attenuate_by_competing_pressure(
                _endogenous_support(oligomeric_register * progress, heptad_state, phase_state),
                register_shift,
            )
            core = _endogenous_support(min(heptad_state, phase_state), register_pairing, closure * progress)
            frustration = _attenuate_by_competing_pressure(
                _endogenous_support(register_frustration, register_shift, heptad_shuffle, phase_shift),
                heptad_state,
                core,
            )
            assembly_dependency = _endogenous_support(core, interface)
            zipper = _endogenous_support(heptad_state, metrics.get("heptad_hydrophobic_periodicity", 0.0))
            basin = {
                "heptad_registered_core": core,
                "heptad_register_context": heptad_state,
                "hydrophobic_repeat_phase": phase_state,
                "parallel_antiparallel_register": register_pairing,
                "oligomeric_coiled_coil_core": core,
                "coiled_coil_assembly_dependency": assembly_dependency,
                "leucine_zipper_context": zipper,
                "register_shift_frustration": frustration,
                "generic_helix_bundle_decoy": _attenuate_by_competing_pressure(
                    _endogenous_support(register_shift, heptad_shuffle, phase_shift, closure),
                    core,
                ),
            }
            segment_compaction = _endogenous_support(core, closure * progress, noise)
            contact_probability = _endogenous_support(core, heptad_state, phase_state)
            interface_readiness = _endogenous_support(core, register_pairing)
            proteostasis_routing = 0.0
        elif mechanism_class == "repeat_solenoid_topology":
            repeat_phase_damage = float((perturbation or {}).get("repeat_phase_damage", 0.0))
            repeat_boundary_mask = float((perturbation or {}).get("repeat_boundary_mask", 0.0))
            repeat_order_shuffle = float((perturbation or {}).get("repeat_order_shuffle", 0.0))
            repeat_phase = strengths["repeat_phase_operator"]
            solenoid_axis = strengths["solenoid_axis_operator"]
            local_repeat = strengths["local_repeat_closure_operator"]
            global_stack = strengths["global_repeat_stack_operator"]
            boundary_frustration = strengths["repeat_boundary_frustration_operator"]
            repeat_unit = _attenuate_by_competing_pressure(
                _endogenous_support(repeat_phase * progress, metrics.get("repeat_signature", 0.0)),
                repeat_phase_damage,
                repeat_order_shuffle,
            )
            axis = _attenuate_by_competing_pressure(
                _endogenous_support(solenoid_axis * progress, repeat_unit),
                repeat_order_shuffle,
                repeat_boundary_mask,
            )
            local_closure = _attenuate_by_competing_pressure(
                _endogenous_support(local_repeat * progress, repeat_unit),
                repeat_boundary_mask,
            )
            phase_alignment = _attenuate_by_competing_pressure(
                _endogenous_support(repeat_phase * progress, min(repeat_unit, axis)),
                repeat_phase_damage,
                repeat_order_shuffle,
            )
            global_topology = _attenuate_by_competing_pressure(
                _endogenous_support(global_stack * progress, axis, local_closure),
                repeat_order_shuffle,
                repeat_boundary_mask,
            )
            boundary_conflict = _attenuate_by_competing_pressure(
                _endogenous_support(boundary_frustration, repeat_boundary_mask, repeat_order_shuffle),
                global_topology,
            )
            lineage_context = _endogenous_support(repeat_unit, axis)
            basin = {
                "repeat_unit_context": repeat_unit,
                "solenoid_axis_context": axis,
                "curved_repeat_stack": _endogenous_support(axis, global_topology),
                "local_repeat_closure": local_closure,
                "global_repeat_topology": global_topology,
                "repeat_phase_alignment": phase_alignment,
                "repeat_boundary_frustration": boundary_conflict,
                "ankyrin_armadillo_tpr_lrr_context": lineage_context,
                "generic_multidomain_decoy": _attenuate_by_competing_pressure(
                    _endogenous_support(repeat_order_shuffle, repeat_boundary_mask, local_repeat),
                    global_topology,
                ),
            }
            segment_compaction = _endogenous_support(local_closure, global_topology, closure * progress, noise)
            contact_probability = _endogenous_support(local_closure, global_topology, axis)
            interface_readiness = _endogenous_support(axis, global_topology)
            proteostasis_routing = 0.0
        elif mechanism_class == "knotted_topology":
            threading_damage = float((perturbation or {}).get("threading_damage", 0.0))
            topology_mask = float((perturbation or {}).get("topology_mask", 0.0))
            long_range_shuffle = float((perturbation or {}).get("long_range_shuffle", 0.0))
            threading = strengths["threading_operator"]
            topology = strengths["topological_closure_operator"]
            long_range = strengths["long_range_threading_operator"]
            slipknot = strengths["slipknot_intermediate_operator"]
            knot_frustration = strengths["knotting_frustration_operator"]
            threading_loop = _attenuate_by_competing_pressure(
                _endogenous_support(threading * progress, interface),
                threading_damage,
                long_range_shuffle,
            )
            long_range_state = _attenuate_by_competing_pressure(
                _endogenous_support(long_range * progress, threading_loop),
                long_range_shuffle,
                topology_mask,
            )
            closure_constraint = _attenuate_by_competing_pressure(
                _endogenous_support(topology * progress, min(threading_loop, long_range_state)),
                topology_mask,
                threading_damage,
            )
            slipknot_state = _attenuate_by_competing_pressure(
                _endogenous_support(slipknot * progress, threading_loop),
                threading_damage,
            )
            knot_core = _endogenous_support(min(threading_loop, closure_constraint), long_range_state, slipknot_state)
            frustration = _attenuate_by_competing_pressure(
                _endogenous_support(knot_frustration, threading_damage, topology_mask, long_range_shuffle),
                knot_core,
            )
            unknotted = _attenuate_by_competing_pressure(
                _endogenous_support(topology_mask, long_range_shuffle, threading_damage, closure),
                knot_core,
            )
            basin = {
                "knot_core_context": knot_core,
                "threading_loop_context": threading_loop,
                "slipknot_intermediate_context": slipknot_state,
                "topological_closure_constraint": closure_constraint,
                "long_range_threading_dependency": long_range_state,
                "knotting_frustration": frustration,
                "unknotted_decoy_dominance": unknotted,
            }
            segment_compaction = _endogenous_support(knot_core, closure_constraint, closure * progress, noise)
            contact_probability = _endogenous_support(knot_core, threading_loop, long_range_state)
            interface_readiness = _endogenous_support(threading_loop, closure_constraint)
            proteostasis_routing = 0.0
        elif mechanism_class == "membrane_multidomain_folding_proteostasis":
            damage = float((perturbation or {}).get("damage", 0.0))
            rescue = float((perturbation or {}).get("rescue", 0.0))
            stability = _attenuate_by_competing_pressure(_endogenous_support(closure * progress, rescue), damage)
            interface_ready = _attenuate_by_competing_pressure(_endogenous_support(interface * progress, rescue), damage)
            routing = _attenuate_by_competing_pressure(_endogenous_support(proteostasis * progress, membrane, rescue), damage)
            basin = {
                "mature_membrane_routed": routing,
                "qc_retained_misfolded": _attenuate_by_competing_pressure(damage, routing, rescue),
                "partial_nbd1_rescue": stability,
            }
            segment_compaction = stability
            contact_probability = _endogenous_support(interface_ready, stability)
            interface_readiness = interface_ready
            proteostasis_routing = routing
        elif mechanism_class == "metamorphic_fold_switching":
            release = float((perturbation or {}).get("release", 0.0))
            alpha_bias = float((perturbation or {}).get("alpha_bias", 0.0))
            beta_bias = float((perturbation or {}).get("beta_bias", 0.0))
            alpha = _attenuate_by_competing_pressure(_endogenous_support(switch, alpha_bias), release, beta_bias, progress)
            beta = _endogenous_support(switch * progress, release, beta_bias)
            basin = {
                "alpha_context_basin": alpha,
                "beta_released_basin": beta,
                "averaged_single_fold": _attenuate_by_competing_pressure(1.0 - switch, alpha, beta),
            }
            segment_compaction = _endogenous_support(switch, closure * progress)
            contact_probability = _endogenous_support(switch, interface)
            interface_readiness = _endogenous_support(interface, switch)
            proteostasis_routing = 0.0
        elif mechanism_class == "short_region_host_interface_hijacking":
            disruption = float((perturbation or {}).get("interface_disruption", 0.0))
            host_ready = _attenuate_by_competing_pressure(_endogenous_support(host * progress, interface), disruption)
            basin = {
                "host_interface_engaged": host_ready,
                "exposed_short_region": _attenuate_by_competing_pressure(disorder, disruption),
                "compact_single_fold": _attenuate_by_competing_pressure(closure * progress, disorder),
            }
            segment_compaction = _endogenous_support(closure * progress)
            contact_probability = _endogenous_support(host_ready)
            interface_readiness = host_ready
            proteostasis_routing = 0.0
        elif mechanism_class == "globular_closure":
            basin = {
                "compact_folded": _endogenous_support(closure * progress),
                "expanded_unfolded": _attenuate_by_competing_pressure(disorder, closure * progress),
            }
            segment_compaction = _endogenous_support(closure * progress)
            contact_probability = _endogenous_support(closure * progress)
            interface_readiness = _endogenous_support(interface * progress)
            proteostasis_routing = 0.0
        elif mechanism_class == "cofactor_ligand_assisted_stabilization":
            cofactor_loss = float((perturbation or {}).get("cofactor_loss", 0.0))
            rescue = float((perturbation or {}).get("rescue", 0.0))
            pocket_ready = _attenuate_by_competing_pressure(_endogenous_support(interface * progress, closure, rescue), cofactor_loss)
            compact = _attenuate_by_competing_pressure(_endogenous_support(closure * progress, pocket_ready), cofactor_loss)
            basin = {
                "ligand_stabilized_basin": pocket_ready,
                "apo_weak_basin": _attenuate_by_competing_pressure(cofactor_loss, pocket_ready, rescue),
                "generic_compact_basin": _attenuate_by_competing_pressure(closure * progress, pocket_ready),
            }
            segment_compaction = compact
            contact_probability = _endogenous_support(compact, pocket_ready)
            interface_readiness = pocket_ready
            proteostasis_routing = 0.0
        elif mechanism_class == "metal_cluster_and_ligand_locked_basin":
            cofactor_loss = float((perturbation or {}).get("cofactor_loss", 0.0))
            coordination_damage = float((perturbation or {}).get("coordination_damage", 0.0))
            rescue = float((perturbation or {}).get("rescue", 0.0))
            coordination_shell = _attenuate_by_competing_pressure(
                _endogenous_support(interface * progress, frustration_strength, rescue),
                coordination_damage,
                cofactor_loss,
            )
            locked_basin = _attenuate_by_competing_pressure(
                _endogenous_support(closure * progress, coordination_shell, rescue),
                cofactor_loss,
            )
            basin = {
                "metal_cluster_locked_basin": coordination_shell,
                "ligand_locked_basin": locked_basin,
                "apo_unlocked_basin": _attenuate_by_competing_pressure(
                    _endogenous_support(cofactor_loss, coordination_damage),
                    locked_basin,
                ),
            }
            segment_compaction = _endogenous_support(closure * progress, locked_basin)
            contact_probability = _endogenous_support(segment_compaction, coordination_shell)
            interface_readiness = coordination_shell
            proteostasis_routing = 0.0
        elif mechanism_class == "assembly_required_folding":
            interface_disruption = float((perturbation or {}).get("interface_disruption", 0.0))
            concentration_rescue = float((perturbation or {}).get("concentration_rescue", 0.0))
            partner_completion = _attenuate_by_competing_pressure(
                _endogenous_support(interface * progress, closure, frustration_strength, concentration_rescue),
                interface_disruption,
            )
            monomer_gap = _attenuate_by_competing_pressure(
                _endogenous_support(frustration_strength, interface_disruption),
                partner_completion,
            )
            basin = {
                "assembly_required_basin": partner_completion,
                "monomer_incomplete_topology": monomer_gap,
                "assembly_ambiguous_basin": _attenuate_by_competing_pressure(interface_disruption, partner_completion),
            }
            segment_compaction = _endogenous_support(closure * progress, partner_completion)
            contact_probability = _endogenous_support(partner_completion, closure)
            interface_readiness = partner_completion
            proteostasis_routing = 0.0
        elif mechanism_class == "oligomerization_controlled_folding":
            interface_disruption = float((perturbation or {}).get("interface_disruption", 0.0))
            concentration_rescue = float((perturbation or {}).get("concentration_rescue", 0.0))
            assembly_ready = _attenuate_by_competing_pressure(
                _endogenous_support(interface * progress, closure, concentration_rescue),
                interface_disruption,
            )
            basin = {
                "assembly_stabilized_basin": assembly_ready,
                "monomer_partial_order": _attenuate_by_competing_pressure(closure, assembly_ready),
                "interface_rejected_basin": interface_disruption,
            }
            segment_compaction = _endogenous_support(closure * progress, assembly_ready)
            contact_probability = _endogenous_support(assembly_ready, closure)
            interface_readiness = assembly_ready
            proteostasis_routing = 0.0
        else:
            basin = {"clean_abstain": 1.0}
            segment_compaction = 0.0
            contact_probability = 0.0
            interface_readiness = 0.0
            proteostasis_routing = 0.0
        exposure = _attenuate_by_competing_pressure(_endogenous_support(metrics["mean_disorder"], disorder), segment_compaction)
        disorder_order = _attenuate_by_competing_pressure(_endogenous_support(disorder, phase), closure, segment_compaction)
        timepoints.append({
            "timepoint": timepoint,
            "residue_exposure": exposure,
            "segment_compaction": segment_compaction,
            "contact_probability": contact_probability,
            "operator_activation": bounded(sum(strengths.values()) / max(1, len(operator_field["operators"]))),
            "frustration": _attenuate_by_competing_pressure(strengths["frustration_operator"], contact_probability),
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
                0.0
                if mechanism_class == "disorder_boundary_and_fold_upon_binding"
                else _attenuate_by_competing_pressure(metrics["mean_interface"], metrics["mean_disorder"])
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
                _attenuate_by_competing_pressure(
                    basin.get("closed_beta_topology", 0.0),
                    basin.get("wrong_beta_topology", 0.0),
                )
                if mechanism_class == "beta_closure_topology"
                else 0.0
            ),
            "closed_beta_ambiguous": bounded(
                max(basin.get("open_beta_sheet_ambiguous", 0.0), basin.get("wrong_beta_topology", 0.0))
                if mechanism_class == "beta_closure_topology"
                else 0.0
            ),
            "strand_register_insufficient": bounded(
                _attenuate_by_competing_pressure(
                    basin.get("open_beta_sheet_ambiguous", 0.0),
                    basin.get("strand_register", 0.0),
                )
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
            "secretory_redox_context": bounded(
                basin.get("secretory_redox_context", 0.0)
                if mechanism_class == "secretory_disulfide_redox_topology"
                else 0.0
            ),
            "disulfide_pairing_topology": bounded(
                basin.get("disulfide_pairing_topology", 0.0)
                if mechanism_class == "secretory_disulfide_redox_topology"
                else 0.0
            ),
            "cysteine_pairing_constraint": bounded(
                min(basin.get("secretory_redox_context", 0.0), basin.get("disulfide_pairing_topology", 0.0))
                if mechanism_class == "secretory_disulfide_redox_topology"
                else 0.0
            ),
            "extracellular_stabilized_fold": bounded(
                basin.get("extracellular_stabilized_fold", 0.0)
                if mechanism_class == "secretory_disulfide_redox_topology"
                else 0.0
            ),
            "glycosylation_context": bounded(
                basin.get("glycosylation_context", 0.0)
                if mechanism_class == "secretory_disulfide_redox_topology"
                else 0.0
            ),
            "redox_mispaired_frustration": bounded(
                basin.get("redox_mispaired_frustration", 0.0)
                if mechanism_class == "secretory_disulfide_redox_topology"
                else 0.0
            ),
            "signal_peptide_removed_context": bounded(
                basin.get("signal_peptide_removed_context", 0.0)
                if mechanism_class == "secretory_disulfide_redox_topology"
                else 0.0
            ),
            "secretory_quality_control": bounded(
                basin.get("secretory_quality_control", 0.0)
                if mechanism_class == "secretory_disulfide_redox_topology"
                else 0.0
            ),
            "signal_peptide_routing_context": bounded(
                basin.get("signal_peptide_routing_context", 0.0)
                if mechanism_class == "signal_peptide_vs_true_tm_routing"
                else 0.0
            ),
            "cleavage_site_context": bounded(
                basin.get("cleavage_site_context", 0.0)
                if mechanism_class == "signal_peptide_vs_true_tm_routing"
                else 0.0
            ),
            "n_terminal_secretory_hydrophobic_patch": bounded(
                basin.get("n_terminal_secretory_hydrophobic_patch", 0.0)
                if mechanism_class == "signal_peptide_vs_true_tm_routing"
                else 0.0
            ),
            "true_transmembrane_span_context": bounded(
                basin.get("true_transmembrane_span_context", 0.0)
                if mechanism_class == "signal_peptide_vs_true_tm_routing"
                else 0.0
            ),
            "single_pass_tm_conflict": bounded(
                basin.get("single_pass_tm_conflict", 0.0)
                if mechanism_class == "signal_peptide_vs_true_tm_routing"
                else 0.0
            ),
            "multi_pass_tm_conflict": bounded(
                basin.get("multi_pass_tm_conflict", 0.0)
                if mechanism_class == "signal_peptide_vs_true_tm_routing"
                else 0.0
            ),
            "secretory_lumenal_routing": bounded(
                basin.get("secretory_lumenal_routing", 0.0)
                if mechanism_class == "signal_peptide_vs_true_tm_routing"
                else 0.0
            ),
            "membrane_insertion_routing": bounded(
                basin.get("membrane_insertion_routing", 0.0)
                if mechanism_class == "signal_peptide_vs_true_tm_routing"
                else 0.0
            ),
            "signal_anchor_ambiguity": bounded(
                basin.get("signal_anchor_ambiguity", 0.0)
                if mechanism_class == "signal_peptide_vs_true_tm_routing"
                else 0.0
            ),
            "heptad_register_context": bounded(
                basin.get("heptad_register_context", 0.0)
                if mechanism_class == "coiled_coil_register_topology"
                else 0.0
            ),
            "hydrophobic_repeat_phase": bounded(
                basin.get("hydrophobic_repeat_phase", 0.0)
                if mechanism_class == "coiled_coil_register_topology"
                else 0.0
            ),
            "parallel_antiparallel_register": bounded(
                basin.get("parallel_antiparallel_register", 0.0)
                if mechanism_class == "coiled_coil_register_topology"
                else 0.0
            ),
            "oligomeric_coiled_coil_core": bounded(
                basin.get("oligomeric_coiled_coil_core", 0.0)
                if mechanism_class == "coiled_coil_register_topology"
                else 0.0
            ),
            "register_shift_frustration": bounded(
                basin.get("register_shift_frustration", 0.0)
                if mechanism_class == "coiled_coil_register_topology"
                else 0.0
            ),
            "coiled_coil_assembly_dependency": bounded(
                basin.get("coiled_coil_assembly_dependency", 0.0)
                if mechanism_class == "coiled_coil_register_topology"
                else 0.0
            ),
            "leucine_zipper_context": bounded(
                basin.get("leucine_zipper_context", 0.0)
                if mechanism_class == "coiled_coil_register_topology"
                else 0.0
            ),
            "repeat_unit_context": bounded(
                basin.get("repeat_unit_context", 0.0)
                if mechanism_class == "repeat_solenoid_topology"
                else 0.0
            ),
            "solenoid_axis_context": bounded(
                basin.get("solenoid_axis_context", 0.0)
                if mechanism_class == "repeat_solenoid_topology"
                else 0.0
            ),
            "curved_repeat_stack": bounded(
                basin.get("curved_repeat_stack", 0.0)
                if mechanism_class == "repeat_solenoid_topology"
                else 0.0
            ),
            "local_repeat_closure": bounded(
                basin.get("local_repeat_closure", 0.0)
                if mechanism_class == "repeat_solenoid_topology"
                else 0.0
            ),
            "global_repeat_topology": bounded(
                basin.get("global_repeat_topology", 0.0)
                if mechanism_class == "repeat_solenoid_topology"
                else 0.0
            ),
            "repeat_phase_alignment": bounded(
                basin.get("repeat_phase_alignment", 0.0)
                if mechanism_class == "repeat_solenoid_topology"
                else 0.0
            ),
            "repeat_boundary_frustration": bounded(
                basin.get("repeat_boundary_frustration", 0.0)
                if mechanism_class == "repeat_solenoid_topology"
                else 0.0
            ),
            "ankyrin_armadillo_tpr_lrr_context": bounded(
                basin.get("ankyrin_armadillo_tpr_lrr_context", 0.0)
                if mechanism_class == "repeat_solenoid_topology"
                else 0.0
            ),
            "knot_core_context": bounded(
                basin.get("knot_core_context", 0.0)
                if mechanism_class == "knotted_topology"
                else 0.0
            ),
            "threading_loop_context": bounded(
                basin.get("threading_loop_context", 0.0)
                if mechanism_class == "knotted_topology"
                else 0.0
            ),
            "slipknot_intermediate_context": bounded(
                basin.get("slipknot_intermediate_context", 0.0)
                if mechanism_class == "knotted_topology"
                else 0.0
            ),
            "topological_closure_constraint": bounded(
                basin.get("topological_closure_constraint", 0.0)
                if mechanism_class == "knotted_topology"
                else 0.0
            ),
            "long_range_threading_dependency": bounded(
                basin.get("long_range_threading_dependency", 0.0)
                if mechanism_class == "knotted_topology"
                else 0.0
            ),
            "knotting_frustration": bounded(
                basin.get("knotting_frustration", 0.0)
                if mechanism_class == "knotted_topology"
                else 0.0
            ),
            "unknotted_decoy_dominance": bounded(
                basin.get("unknotted_decoy_dominance", 0.0)
                if mechanism_class == "knotted_topology"
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
                    "probability": _endogenous_support(left["aromatic_density"], right["aromatic_density"], _operator_strength(operator_field, "phase_operator")),
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
                    "probability": _endogenous_support(left["disorder_density"], right["interface_density"], _operator_strength(operator_field, "interface_operator")),
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
                    "probability": _endogenous_support(left["beta_propensity_density"], right["beta_propensity_density"], _operator_strength(operator_field, "closure_operator")),
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
                    "probability": _endogenous_support(
                        left["interface_density"],
                        right["interface_density"],
                        _operator_strength(operator_field, "interface_operator"),
                        _operator_strength(operator_field, "dual_basin_switch_operator"),
                    ),
                    "interaction_type": "interdomain_allosteric_lock",
                })
    elif mechanism_class == "secretory_disulfide_redox_topology":
        cysteine_top = _strong_segments(sequence_field, "cysteine_density", limit=6)
        if len(cysteine_top) < 2:
            cysteine_top = _strong_segments(sequence_field, "interface_density", limit=4)
        for index in range(0, max(0, len(cysteine_top) - 1), 2):
            left = cysteine_top[index]
            right = cysteine_top[index + 1]
            if left["segment_id"] == right["segment_id"]:
                continue
            pairs.append({
                "segment_a": left["segment_id"],
                "segment_b": right["segment_id"],
                "probability": _endogenous_support(
                    left["cysteine_density"],
                    right["cysteine_density"],
                    _operator_strength(operator_field, "disulfide_pairing_operator"),
                    _operator_strength(operator_field, "secretory_redox_operator"),
                ),
                "interaction_type": "secretory_disulfide_pairing_contact",
            })
    elif mechanism_class == "signal_peptide_vs_true_tm_routing":
        n_terminal = segments[0]
        membrane_top = _strong_segments(sequence_field, "membrane_density", limit=2)
        pairs.append({
            "segment_a": n_terminal["segment_id"],
            "segment_b": "signal_peptidase_translocon_route",
            "probability": _endogenous_support(
                n_terminal["membrane_density"],
                _operator_strength(operator_field, "signal_peptide_routing_operator"),
                _operator_strength(operator_field, "cleavage_context_operator"),
            ),
            "interaction_type": "signal_peptide_secretory_routing_contact",
        })
        for segment in membrane_top:
            if segment["segment_id"] == n_terminal["segment_id"]:
                continue
            pairs.append({
                "segment_a": n_terminal["segment_id"],
                "segment_b": segment["segment_id"],
                "probability": _endogenous_support(
                    n_terminal["membrane_density"],
                    segment["membrane_density"],
                    _operator_strength(operator_field, "tm_insertion_operator"),
                ),
                "interaction_type": "signal_peptide_true_tm_boundary_probe",
            })
    elif mechanism_class == "coiled_coil_register_topology":
        hydrophobic_top = _strong_segments(sequence_field, "hydrophobic_density", limit=4)
        if len(hydrophobic_top) < 2:
            hydrophobic_top = _strong_segments(sequence_field, "interface_density", limit=4)
        for left in hydrophobic_top[:2]:
            for right in hydrophobic_top[2:]:
                if left["segment_id"] == right["segment_id"]:
                    continue
                pairs.append({
                    "segment_a": left["segment_id"],
                    "segment_b": right["segment_id"],
                    "probability": _endogenous_support(
                        _operator_strength(operator_field, "heptad_register_operator"),
                        _operator_strength(operator_field, "coiled_coil_interface_operator"),
                        left["hydrophobic_density"],
                        right["hydrophobic_density"],
                    ),
                    "interaction_type": "heptad_registered_coiled_coil_contact",
                })
    elif mechanism_class == "repeat_solenoid_topology":
        repeat_top = sorted(
            segments,
            key=lambda row: row["interface_density"] + row["beta_propensity_density"] + row["pro_gly_density"],
            reverse=True,
        )[:4]
        for left, right in zip(repeat_top, repeat_top[1:]):
            if left["segment_id"] == right["segment_id"]:
                continue
            pairs.append({
                "segment_a": left["segment_id"],
                "segment_b": right["segment_id"],
                "probability": _endogenous_support(
                    _operator_strength(operator_field, "repeat_phase_operator"),
                    _operator_strength(operator_field, "solenoid_axis_operator"),
                    left["interface_density"],
                    right["interface_density"],
                ),
                "interaction_type": "phase_aligned_repeat_solenoid_contact",
            })
    elif mechanism_class == "knotted_topology":
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
                    "probability": _endogenous_support(
                        _operator_strength(operator_field, "threading_operator"),
                        _operator_strength(operator_field, "topological_closure_operator"),
                        left["interface_density"],
                        right["interface_density"],
                    ),
                    "interaction_type": "threaded_knot_topology_contact",
                })
    elif mechanism_class == "short_region_host_interface_hijacking":
        cterm = segments[-1]
        pairs.append({
            "segment_a": cterm["segment_id"],
            "segment_b": "host_RAE1_NUP98_surface",
            "probability": _endogenous_support(cterm["interface_density"], _operator_strength(operator_field, "host_hijack_operator")),
            "interaction_type": "host_interface_capture",
        })
    elif mechanism_class == "metamorphic_fold_switching":
        pairs.append({
            "segment_a": segments[0]["segment_id"],
            "segment_b": segments[-1]["segment_id"],
            "probability": _endogenous_support(
                segments[0]["interface_density"],
                segments[-1]["interface_density"],
                _operator_strength(operator_field, "dual_basin_switch_operator"),
            ),
            "interaction_type": "state_basin_contact_rewrite",
        })
    elif mechanism_class == "membrane_multidomain_folding_proteostasis":
        top_membrane = _strong_segments(sequence_field, "membrane_density", limit=2)
        focus = max(segments, key=lambda row: _endogenous_support(row["interface_density"], row["membrane_density"]))
        for segment in top_membrane:
            pairs.append({
                "segment_a": focus["segment_id"],
                "segment_b": segment["segment_id"],
                "probability": _endogenous_support(focus["interface_density"], segment["membrane_density"], _operator_strength(operator_field, "interface_operator")),
                "interaction_type": "interdomain_membrane_interface",
            })
    elif mechanism_class == "cofactor_ligand_assisted_stabilization":
        top = _strong_segments(sequence_field, "interface_density", limit=3)
        for segment in top:
            pairs.append({
                "segment_a": segment["segment_id"],
                "segment_b": "ligand_or_cofactor_pocket",
                "probability": _endogenous_support(segment["interface_density"], _operator_strength(operator_field, "interface_operator")),
                "interaction_type": "cofactor_stabilized_interface",
            })
    elif mechanism_class == "metal_cluster_and_ligand_locked_basin":
        top = _strong_segments(sequence_field, "interface_density", limit=3)
        for segment in top:
            pairs.append({
                "segment_a": segment["segment_id"],
                "segment_b": "metal_cluster_or_ligand_locked_pocket",
                "probability": _endogenous_support(segment["interface_density"], segment["aromatic_density"], _operator_strength(operator_field, "interface_operator")),
                "interaction_type": "coordination_shell_locked_interface",
            })
    elif mechanism_class == "assembly_required_folding":
        top = _strong_segments(sequence_field, "interface_density", limit=3)
        for segment in top:
            pairs.append({
                "segment_a": segment["segment_id"],
                "segment_b": "partner_completed_core_interface",
                "probability": _endogenous_support(segment["interface_density"], _operator_strength(operator_field, "interface_operator")),
                "interaction_type": "assembly_required_partner_completion",
            })
    elif mechanism_class == "oligomerization_controlled_folding":
        top = _strong_segments(sequence_field, "interface_density", limit=3)
        for segment in top:
            pairs.append({
                "segment_a": segment["segment_id"],
                "segment_b": "partner_copy_interface",
                "probability": _endogenous_support(segment["interface_density"], _operator_strength(operator_field, "interface_operator")),
                "interaction_type": "assembly_stabilized_interface",
            })
    else:
        top = _strong_segments(sequence_field, "hydrophobic_density", limit=4)
        for left in top[:2]:
            for right in top[2:]:
                pairs.append({
                    "segment_a": left["segment_id"],
                    "segment_b": right["segment_id"],
                    "probability": _endogenous_support(left["hydrophobic_density"], right["hydrophobic_density"], _operator_strength(operator_field, "closure_operator")),
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
    if delta > 0.0:
        return "increase"
    if delta < 0.0:
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


def _legacy_acceptance_from_self_decision(judge: dict[str, Any]) -> dict[str, Any]:
    final = judge["final_self_decision"]
    if final in {"accepted", "accepted_with_caution"}:
        decision = "accepted"
    elif final == "blocked_for_leakage":
        decision = "blocked_for_leakage"
    else:
        decision = "abstain_recommended"
    return {
        "kind": "PROTEIN_ESPERANTO_SELF_DECISION_COMPAT_ACCEPTANCE_VIEW_v0",
        "zero_failed_accepted_required": True,
        "acceptance_decision": decision,
        "firewall_reason": judge["self_decision_reason"],
        "blocked_reasons": judge["blocked_reasons"],
        "unknown_word_signals": [judge["missing_word_candidate"]] if judge.get("missing_word_candidate") else [],
        "missing_esperanto_word": judge.get("missing_word_candidate"),
        "known_mechanism_class": judge["top_mechanism"],
        "known_multidomain_word": judge.get("known_multidomain_word"),
        "known_beta_topology_word": judge.get("known_beta_topology_word"),
        "known_secretory_disulfide_word": judge.get("known_secretory_disulfide_word"),
        "known_signal_peptide_word": judge.get("known_signal_peptide_word"),
        "known_coiled_coil_word": judge.get("known_coiled_coil_word"),
        "known_repeat_solenoid_word": judge.get("known_repeat_solenoid_word"),
        "known_knotted_topology_word": judge.get("known_knotted_topology_word"),
        "operator_activation": judge.get("operator_activation", 0.0),
        "coordinate_truth_used_before_prediction": judge["coordinate_truth_used_before_prediction"],
        "folding_problem_solved": False,
    }


def self_decision_judge(
    *,
    sequence_field: dict[str, Any],
    evidence_manifest: dict[str, Any],
    sources: list[dict[str, Any]],
    mechanism: dict[str, Any],
    operator_field: dict[str, Any],
    trajectory: dict[str, Any],
    physical_calibration_inputs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    text = _allowed_source_text(sources, evidence_manifest)
    mechanism_class = mechanism["mechanism_class"]
    blocked_reasons: list[str] = []
    if evidence_manifest["coordinate_derived_source_count_before_prediction"] > 0:
        blocked_reasons.append("coordinate_derived_source_before_prediction")
    if evidence_manifest["internal_runtime_source_count_for_prediction"] > 0:
        blocked_reasons.append("internal_runtime_source_before_prediction")
    if evidence_manifest["holdout_opened_before_seal"]:
        blocked_reasons.append("holdout_opened_before_seal")
    readouts = _mechanism_competition_readouts(
        text=text,
        sequence_field=sequence_field,
        mechanism=mechanism,
        operator_field=operator_field,
        trajectory=trajectory,
    )
    top_readout = readouts[0] if readouts else None
    selected_readout = next((row for row in readouts if row["grammar"] == mechanism_class), None)
    masking_rows = _masking_probe(
        sequence_field=sequence_field,
        evidence_manifest=evidence_manifest,
        sources=sources,
        top_mechanism=mechanism_class,
    )
    wrong_grammar_rows = _wrong_grammar_challenge(
        sequence_field=sequence_field,
        evidence_manifest=evidence_manifest,
        sources=sources,
        mechanism=mechanism,
        readouts=readouts,
    )
    wrong_grammar_separation = (
        "wrong_grammars_fail"
        if wrong_grammar_rows and all(row["forced_grammar_rejected"] for row in wrong_grammar_rows)
        else "no_wrong_grammar_needed"
        if not wrong_grammar_rows
        else "wrong_grammar_competes"
    )
    counterfactuals = _counterfactual_sequence_controls(sequence_field, mechanism_class)
    operator_basis_probe = _operator_basis_stability_probe(
        sequence_field=sequence_field,
        operator_field=operator_field,
        mechanism_class=mechanism_class,
        baseline_trajectory=trajectory,
    )
    temporal_probe = _temporal_binding_probe(mechanism_class, trajectory)
    physical_gate = _physical_grounding_gate(
        sequence_field=sequence_field,
        evidence_manifest=evidence_manifest,
        mechanism=mechanism,
        operator_field=operator_field,
        trajectory=trajectory,
        temporal_probe=temporal_probe,
        operator_basis_probe=operator_basis_probe,
        physical_calibration_inputs=physical_calibration_inputs,
    )
    cross_view_binding = _cross_view_binding_probe(selected_readout)
    contradictions = _internal_contradictions(
        mechanism=mechanism,
        text=text,
        trajectory=trajectory,
        readouts=readouts,
        operator_basis_probe=operator_basis_probe,
        cross_view_binding=cross_view_binding,
        temporal_probe=temporal_probe,
    )
    missing_word = _missing_word_candidate(readouts)
    selected_views = selected_readout["view_sources"] if selected_readout else []
    internal_consensus = _explicit_dominance_law(
        mechanism_class=mechanism_class,
        selected_readout=selected_readout,
        top_readout=top_readout,
        missing_word=missing_word,
        cross_view_binding=cross_view_binding,
        operator_basis_probe=operator_basis_probe,
        temporal_probe=temporal_probe,
    )
    if masking_rows:
        nondefining_flips = [
            row
            for row in masking_rows
            if row["masked_family"] != mechanism_class and not row["selected_mechanism_preserved"]
        ]
        masking_stability = "unstable_under_nondefining_mask" if nondefining_flips else "stable_or_definition_sensitive_under_masking"
    else:
        masking_stability = "no_evidence_family_to_mask"
    if blocked_reasons:
        final_decision = "blocked_for_leakage"
        reason = "blocked_prediction_boundary_violation"
    elif missing_word:
        final_decision = "clean_abstain_missing_word"
        reason = "self_decision_candidate_word_competes_without_learned_grammar"
    elif mechanism_class == "insufficient_evidence_clean_abstain":
        final_decision = "clean_abstain_low_internal_consensus"
        reason = "self_decision_no_dominant_learned_mechanism"
    elif contradictions:
        final_decision = "clean_abstain_conflict"
        reason = "self_decision_unresolved_internal_contradiction"
    elif wrong_grammar_separation == "wrong_grammar_competes":
        final_decision = "clean_abstain_low_internal_consensus"
        reason = "self_decision_wrong_grammar_competes"
    elif internal_consensus == "single_dominant_learned_mechanism_bound_across_views":
        final_decision = "accepted"
        reason = "self_decision_internal_views_agree"
    else:
        final_decision = "clean_abstain_low_internal_consensus"
        reason = "self_decision_internal_views_do_not_converge"
    final_state = trajectory["final_state_summary"]
    operator_activation = final_state.get("operator_activation", 0.0)
    judge = {
        "kind": "PROTEIN_ESPERANTO_SELF_DECISION_JUDGE_v0",
        "zero_failed_accepted_required": True,
        "top_mechanism": mechanism_class,
        "natural_mechanism_class": mechanism["natural_mechanism_class"],
        "runner_up_mechanisms": [row["grammar"] for row in readouts if row["grammar"] != mechanism_class][:3],
        "mechanism_competition": readouts[:8],
        "internal_consensus": internal_consensus,
        "dominance_law": internal_consensus,
        "selected_mechanism_views": selected_views,
        "cross_view_binding": cross_view_binding["cross_view_binding"],
        "cross_view_binding_probe": cross_view_binding,
        "masking_stability": masking_stability,
        "evidence_masking": masking_rows,
        "wrong_grammar_separation": wrong_grammar_separation,
        "wrong_grammar_challenge": wrong_grammar_rows,
        "counterfactual_separation": counterfactuals["interpretation"],
        "counterfactual_sequence_controls": counterfactuals,
        "operator_basis_stability": operator_basis_probe["operator_basis_stability"],
        "operator_basis_stability_probe": operator_basis_probe,
        "coefficient_perturbation_stability": operator_basis_probe["operator_basis_stability"],
        "coefficient_perturbation_probe": operator_basis_probe,
        "coefficient_probe_mode": operator_basis_probe["coefficient_probe_mode"],
        "temporal_binding": temporal_probe["temporal_binding"],
        "temporal_binding_probe": temporal_probe,
        "physics_grounding_status": "coarse_operator_heuristic_not_atomistic_physics",
        "coefficient_source": "heuristic_internal_operator_weights_not_physical_force_constants",
        "physical_basis_claim_allowed": False,
        "physical_grounding_gate": physical_gate,
        "physical_grounding_status": physical_gate["physical_grounding_status"],
        "physical_backend_available": physical_gate["backend_probe"]["backend_available"],
        "real_physical_calibration_inputs_used": physical_gate["real_physical_calibration_inputs_used"],
        "real_physical_calibration_kind": physical_gate["real_physical_calibration_kind"],
        "real_physical_calibration_row_count": physical_gate["real_physical_calibration_row_count"],
        "real_physical_calibration_hash": physical_gate["real_physical_calibration_hash"],
        "contradiction_count": len(contradictions),
        "contradictions": contradictions,
        "missing_word_candidate": missing_word,
        "candidate_grammar_scope": [
            {
                "grammar": grammar,
                "grammar_status": spec["grammar_status"],
                "acceptance_role": "clean_abstain_until_revision_implements_grammar",
                "sequence_signal": spec["sequence_signal"],
            }
            for grammar, spec in SELF_DECISION_CANDIDATE_GRAMMARS.items()
        ],
        "final_self_decision": final_decision,
        "self_decision_reason": reason,
        "blocked_reasons": blocked_reasons,
        "known_multidomain_word": mechanism.get("selected_multidomain_word"),
        "known_beta_topology_word": mechanism.get("selected_beta_topology_word"),
        "known_secretory_disulfide_word": mechanism.get("selected_secretory_disulfide_word"),
        "known_signal_peptide_word": mechanism.get("selected_signal_peptide_word"),
        "known_coiled_coil_word": mechanism.get("selected_coiled_coil_word"),
        "known_repeat_solenoid_word": mechanism.get("selected_repeat_solenoid_word"),
        "known_knotted_topology_word": mechanism.get("selected_knotted_topology_word"),
        "operator_activation": operator_activation,
        "coordinate_truth_used_before_prediction": evidence_manifest["coordinate_truth_used_before_prediction"],
        "folding_problem_solved": False,
    }
    legacy = _legacy_acceptance_from_self_decision(judge)
    judge["legacy_acceptance_view"] = legacy
    for key, value in legacy.items():
        if key != "kind":
            judge[key] = value
    return judge


def acceptance_firewall(
    *,
    sequence_field: dict[str, Any],
    evidence_manifest: dict[str, Any],
    sources: list[dict[str, Any]],
    mechanism: dict[str, Any],
    trajectory: dict[str, Any],
    physical_calibration_inputs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    operator_field = build_operator_field(
        sequence_field=sequence_field,
        mechanism=mechanism,
        evidence_manifest=evidence_manifest,
    )
    judge = self_decision_judge(
        sequence_field=sequence_field,
        evidence_manifest=evidence_manifest,
        sources=sources,
        mechanism=mechanism,
        operator_field=operator_field,
        trajectory=trajectory,
        physical_calibration_inputs=physical_calibration_inputs,
    )
    return _legacy_acceptance_from_self_decision(judge)


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
        {
            "operator": "disulfide_pairing_operator",
            "custom_force_term": "coarse cysteine-pair proximity bias for secretory disulfide grammar",
            "guard": "does not encode native disulfide geometry; only compares biased versus unbiased sealed observables",
        },
        {
            "operator": "secretory_redox_operator",
            "custom_force_term": "secretory redox-state observable bias outside atomistic chemistry",
            "guard": "redox state is reported as a coarse proxy, not a quantum or force-field disulfide formation claim",
        },
        {
            "operator": "signal_peptide_routing_operator",
            "custom_force_term": "coarse N-terminal secretory-routing bias for cleavable signal peptide grammar",
            "guard": "does not encode a native membrane path; compares signal-routing observables against matched TM decoys",
        },
        {
            "operator": "tm_insertion_operator",
            "custom_force_term": "coarse competing membrane-insertion observable for signal/TM boundary checks",
            "guard": "true TM insertion remains a falsifier for signal-peptide acceptance, not a hidden label",
        },
        {
            "operator": "cleavage_context_operator",
            "custom_force_term": "coarse cleavage-context observable for signal peptidase routing",
            "guard": "cleavage context must come from sealed non-coordinate evidence and matched controls",
        },
        {
            "operator": "secretory_routing_operator",
            "custom_force_term": "coarse secretory-lumenal routing observable",
            "guard": "secretory route is a language observable unless independent physical holdouts earn a physical claim",
        },
        {
            "operator": "heptad_register_operator",
            "custom_force_term": "coarse heptad/register alignment observable derived from sealed sequence and evidence",
            "guard": "register bias is compared with heptad-shuffled controls and must not encode native coordinates",
        },
        {
            "operator": "coiled_coil_interface_operator",
            "custom_force_term": "coarse coiled-coil interface observable",
            "guard": "ordinary helix bundle and generic assembly controls must fail the coiled-coil claim",
        },
        {
            "operator": "oligomeric_register_operator",
            "custom_force_term": "coarse partner-register observable for coiled-coil cores",
            "guard": "partner-copy context remains sealed non-coordinate evidence",
        },
        {
            "operator": "register_shift_frustration_operator",
            "custom_force_term": "coarse wrong-register conflict observable",
            "guard": "accepted rows must beat register-shift controls without static observable cutoffs",
        },
        {
            "operator": "repeat_phase_operator",
            "custom_force_term": "coarse repeat-unit phase observable",
            "guard": "repeat order controls must reduce phase support without native contact leakage",
        },
        {
            "operator": "solenoid_axis_operator",
            "custom_force_term": "coarse solenoid-axis continuity observable",
            "guard": "generic multidomain and beta controls remain falsifiers",
        },
        {
            "operator": "local_repeat_closure_operator",
            "custom_force_term": "coarse adjacent-repeat closure observable",
            "guard": "local closure cannot by itself claim a physical fold",
        },
        {
            "operator": "global_repeat_stack_operator",
            "custom_force_term": "coarse global repeat-stack observable",
            "guard": "global stack is scored by paired controls, not fixed thresholds",
        },
        {
            "operator": "repeat_boundary_frustration_operator",
            "custom_force_term": "coarse repeat-boundary conflict observable",
            "guard": "boundary masking should expose conflict instead of forcing acceptance",
        },
        {
            "operator": "threading_operator",
            "custom_force_term": "coarse threading-loop observable for knot/slipknot grammar",
            "guard": "explicit non-coordinate topology context is required; sequence alone cannot validate the knot",
        },
        {
            "operator": "topological_closure_operator",
            "custom_force_term": "coarse topological-closure observable",
            "guard": "topology-masked controls must fail before any accepted knot grammar is supported",
        },
        {
            "operator": "long_range_threading_operator",
            "custom_force_term": "coarse long-range threading-dependency observable",
            "guard": "long-range order shuffles are paired falsifiers",
        },
        {
            "operator": "slipknot_intermediate_operator",
            "custom_force_term": "coarse slipknot-intermediate observable",
            "guard": "intermediate support is a language observable, not an atomistic path claim",
        },
        {
            "operator": "knotting_frustration_operator",
            "custom_force_term": "coarse unknotted-decoy/frustration observable",
            "guard": "accepted knot rows must beat unknotted topology controls and physical claims stay blocked",
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
    physical_calibration_inputs: dict[str, Any] | None = None,
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
    judge = self_decision_judge(
        sequence_field=sequence_field,
        evidence_manifest=evidence_manifest,
        sources=sources,
        mechanism=mechanism,
        operator_field=operator_field,
        trajectory=trajectory,
        physical_calibration_inputs=physical_calibration_inputs,
    )
    acceptance_view = judge["legacy_acceptance_view"]
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
        "epistemological_status": protein_esperanto_epistemological_status(),
        "input_evidence_manifest": {
            "source_count": len(sources),
            "source_ids": [str(source.get("source_id", source.get("accession", "unknown_source"))) for source in sources],
        },
        "evidence_manifest": evidence_manifest,
        "selected_mechanism_grammar": mechanism,
        "self_decision_judge": judge,
        "acceptance_firewall": acceptance_view,
        "physical_grounding_gate": judge["physical_grounding_gate"],
        "physical_calibration_input_summary": judge["physical_grounding_gate"]["physical_calibration_input_summary"],
        "operator_field": operator_field,
        "initial_sequence_field_map": sequence_field,
        "trajectory_summary": trajectory,
        "operator_state_propagation_summary": trajectory,
        "predicted_contact_interaction_probability_map": trajectory["predicted_contact_interaction_probability_map"],
        "hypothesized_interaction_language_map": trajectory["predicted_contact_interaction_probability_map"],
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
                _endogenous_support(row["low_complexity_density"], row["aromatic_density"], row["disorder_density"])
                for row in segments
            ),
            default=0.0,
        )
        support = _endogenous_support(field["low_complexity_density"], field["aromatic_density"], field["mean_disorder"], local)
    elif mechanism == "disorder_boundary_and_fold_upon_binding":
        local_disorder = max((_endogenous_support(row["low_complexity_density"], row["disorder_density"]) for row in segments), default=0.0)
        local_interface = max((row["interface_density"] for row in segments), default=0.0)
        support = _endogenous_support(field["low_complexity_density"], field["mean_disorder"], local_disorder, local_interface, activation)
    elif mechanism == "beta_closure_topology":
        local_beta = max((row["beta_propensity_density"] for row in segments), default=0.0)
        local_interface = max((row["interface_density"] for row in segments), default=0.0)
        local_aromatic = max((row["aromatic_density"] for row in segments), default=0.0)
        support = _endogenous_support(
            field["hydrophobic_density"],
            field["beta_propensity_density"],
            field["mean_interface"],
            local_beta,
            max(local_interface, local_aromatic),
            activation,
        )
    elif mechanism == "multidomain_allosteric_architecture":
        local_interface = max((row["interface_density"] for row in segments), default=0.0)
        local_hydrophobic = max((row["hydrophobic_density"] for row in segments), default=0.0)
        local_aromatic = max((row["aromatic_density"] for row in segments), default=0.0)
        support = _endogenous_support(
            field["hydrophobic_density"],
            field["mean_interface"],
            local_interface,
            max(local_hydrophobic, local_aromatic),
            activation,
        )
    elif mechanism == "secretory_disulfide_redox_topology":
        local_cysteine = max((row["cysteine_density"] for row in segments), default=0.0)
        local_interface = max((row["interface_density"] for row in segments), default=0.0)
        support = _endogenous_support(
            field["cysteine_density"],
            field["n_terminal_hydrophobic_signal"],
            local_cysteine,
            local_interface,
            activation,
        )
    elif mechanism == "signal_peptide_vs_true_tm_routing":
        n_terminal_segments = segments[:2]
        internal_segments = segments[2:]
        n_terminal_membrane = max((row["membrane_density"] for row in n_terminal_segments), default=0.0)
        internal_membrane = max((row["membrane_density"] for row in internal_segments), default=0.0)
        n_terminal_interface = max((row["interface_density"] for row in n_terminal_segments), default=0.0)
        route_separation = max(0.0, max(field["n_terminal_hydrophobic_signal"], n_terminal_membrane) - internal_membrane)
        support = _endogenous_support(
            field["n_terminal_hydrophobic_signal"],
            n_terminal_membrane,
            route_separation,
            n_terminal_interface,
            activation,
        )
    elif mechanism == "coiled_coil_register_topology":
        local_hydrophobic = max((row["hydrophobic_density"] for row in segments), default=0.0)
        local_interface = max((row["interface_density"] for row in segments), default=0.0)
        support = _endogenous_support(
            field.get("heptad_hydrophobic_periodicity", 0.0),
            field["hydrophobic_density"],
            field["mean_interface"],
            local_hydrophobic,
            local_interface,
            activation,
        )
    elif mechanism == "repeat_solenoid_topology":
        local_repeat_like = max(
            (
                _endogenous_support(row["interface_density"], row["beta_propensity_density"], row["pro_gly_density"])
                for row in segments
            ),
            default=0.0,
        )
        support = _endogenous_support(
            field.get("repeat_signature", 0.0),
            field["mean_interface"],
            field["beta_propensity_density"],
            local_repeat_like,
            activation,
        )
    elif mechanism == "knotted_topology":
        local_interface = max((row["interface_density"] for row in segments), default=0.0)
        long_range_span = _endogenous_support(segments[0]["interface_density"], segments[-1]["interface_density"]) if segments else 0.0
        support = _endogenous_support(
            field["hydrophobic_density"],
            field["mean_interface"],
            field["aromatic_density"],
            local_interface,
            long_range_span,
            activation,
        )
    elif mechanism == "membrane_multidomain_folding_proteostasis":
        local = max((row["membrane_density"] for row in segments), default=0.0)
        support = _endogenous_support(field["mean_membrane"], local, activation)
    elif mechanism == "metamorphic_fold_switching":
        support = activation
    elif mechanism == "short_region_host_interface_hijacking":
        cterminal = segments[-2:] if len(segments) >= 2 else segments
        local = _avg(row["interface_density"] for row in cterminal)
        support = _endogenous_support(field["mean_interface"], local, activation)
    elif mechanism == "cofactor_ligand_assisted_stabilization":
        local = max((row["interface_density"] for row in segments), default=0.0)
        support = _endogenous_support(field["hydrophobic_density"], field["mean_interface"], local, activation)
    elif mechanism == "metal_cluster_and_ligand_locked_basin":
        local = max((_endogenous_support(row["interface_density"], row["aromatic_density"]) for row in segments), default=0.0)
        support = _endogenous_support(field["hydrophobic_density"], field["mean_interface"], local, activation)
    elif mechanism == "assembly_required_folding":
        local = max((row["interface_density"] for row in segments), default=0.0)
        membrane_like = max((row["membrane_density"] for row in segments), default=0.0)
        support = _endogenous_support(field["hydrophobic_density"], field["mean_interface"], local, membrane_like, activation)
    elif mechanism == "oligomerization_controlled_folding":
        local = max((row["interface_density"] for row in segments), default=0.0)
        support = _endogenous_support(field["hydrophobic_density"], field["mean_interface"], local, activation)
    elif mechanism == "globular_closure":
        local = max((row["hydrophobic_density"] for row in segments), default=0.0)
        support = _endogenous_support(field["hydrophobic_density"], local, activation)
    else:
        support = 0.0
    return _endogenous_support(activation, support)
