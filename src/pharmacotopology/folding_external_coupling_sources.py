from __future__ import annotations

from dataclasses import dataclass


EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID = (
    "EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_V0"
)

EXTERNAL_COUPLING_TRACE_LOOP_REPORT_KIND = (
    "external_evolutionary_coupling_trace_loop_benchmark_v0"
)
EXTERNAL_COUPLING_TRACE_LOOP_CERTIFICATE_KIND = (
    "external_evolutionary_coupling_trace_loop_certificate_v0"
)

ACCEPTED_EXTERNAL_COUPLING_SOURCE_KINDS = frozenset(
    {
        "external_msa_dca_plmc_v1",
        "external_msa_dca_ccmpred_v1",
        "external_evcouplings_sequence_covariation_v1",
        "external_pfam_msa_dca_v1",
        "external_uniref_msa_dca_v1",
    }
)

REJECTED_EXTERNAL_COUPLING_SOURCE_KINDS = frozenset(
    {
        "coordinate_native_contacts",
        "pdb_contact_map",
        "alphafold_distance_map",
        "supervised_structure_contact_predictor",
        "manual_contact_annotation",
        "oracle_from_benchmark_coordinates",
        "oracle_from_native_contacts",
    }
)

EXTERNAL_COUPLING_ROW_STATUSES = frozenset(
    {
        "external_couplings_available",
        "external_couplings_rejected_low_depth",
        "external_couplings_rejected_low_coverage",
        "external_couplings_rejected_mapping_ambiguous",
        "external_couplings_rejected_coordinate_taint",
    }
)


@dataclass(frozen=True)
class ExternalCouplingQualityPolicy:
    target_coverage_min: float
    focus_sequence_mapping_confidence_min: float
    effective_sequence_count_over_length_min: float
    require_top_l_couplings: bool = True


EXPLORATORY_EXTERNAL_COUPLING_POLICY = ExternalCouplingQualityPolicy(
    target_coverage_min=0.70,
    focus_sequence_mapping_confidence_min=0.98,
    effective_sequence_count_over_length_min=1.0,
)

SERIOUS_EXTERNAL_COUPLING_POLICY = ExternalCouplingQualityPolicy(
    target_coverage_min=0.70,
    focus_sequence_mapping_confidence_min=0.98,
    effective_sequence_count_over_length_min=5.0,
)


def accepted_external_source_kind(source_kind: str) -> bool:
    return source_kind in ACCEPTED_EXTERNAL_COUPLING_SOURCE_KINDS


def rejected_external_source_kind(source_kind: str) -> bool:
    return source_kind in REJECTED_EXTERNAL_COUPLING_SOURCE_KINDS
