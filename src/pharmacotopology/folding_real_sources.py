from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Sequence


@dataclass(frozen=True)
class FoldingBenchmarkSource:
    source_id: str
    display_name: str
    role: str
    accepted_prefixes: tuple[str, ...]
    source_url: str
    use_boundary: str
    clinical_use_allowed: bool = False
    drug_design_created: bool = False
    molecule_generated: bool = False
    protein_sequence_design_created: bool = False
    folding_solution_claim_created: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class FoldingBenchmarkSourceReview:
    valid: bool
    sources_reviewed: int
    source_ids: tuple[str, ...]
    accepted_prefixes: tuple[str, ...]
    violations: tuple[str, ...]
    clinical_use_allowed: bool = False
    drug_design_created: bool = False
    molecule_generated: bool = False
    protein_sequence_design_created: bool = False
    folding_solution_claim_created: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


REAL_FOLDING_BENCHMARK_SOURCES: tuple[FoldingBenchmarkSource, ...] = (
    FoldingBenchmarkSource(
        source_id="RCSB_PDB",
        display_name="RCSB Protein Data Bank",
        role="experimental_structure_and_sequence_reference",
        accepted_prefixes=("pdb:",),
        source_url="https://data.rcsb.org/",
        use_boundary="external_structure_reference_not_prediction_truth_claim",
    ),
    FoldingBenchmarkSource(
        source_id="CATH",
        display_name="CATH Protein Structure Classification",
        role="class_architecture_topology_homology_labels",
        accepted_prefixes=("cath:",),
        source_url="https://www.cathdb.info/",
        use_boundary="external_fold_class_label_reference",
    ),
    FoldingBenchmarkSource(
        source_id="SCOPe",
        display_name="SCOPe Structural Classification of Proteins",
        role="structural_fold_domain_label_reference",
        accepted_prefixes=("scop:",),
        source_url="https://scop.berkeley.edu/",
        use_boundary="external_fold_class_label_reference",
    ),
    FoldingBenchmarkSource(
        source_id="DisProt",
        display_name="DisProt",
        role="intrinsic_disorder_reference",
        accepted_prefixes=("disprot:",),
        source_url="https://disprot.org/",
        use_boundary="external_disorder_label_reference",
    ),
    FoldingBenchmarkSource(
        source_id="CASP",
        display_name="CASP-style blind assessment",
        role="blind_benchmark_protocol_reference",
        accepted_prefixes=("casp:",),
        source_url="https://predictioncenter.org/",
        use_boundary="benchmark_protocol_reference_not_training_data",
    ),
    FoldingBenchmarkSource(
        source_id="AlphaFold_DB",
        display_name="AlphaFold Database",
        role="comparison_baseline_not_final_truth",
        accepted_prefixes=("afdb:", "alphafold:"),
        source_url="https://alphafold.ebi.ac.uk/",
        use_boundary="predicted_structure_baseline_not_experimental_proof",
    ),
)

STRATIFIED_500_TARGETS: dict[str, int] = {
    "mostly_alpha_or_compact": 100,
    "mostly_beta_or_long_range": 100,
    "alpha_beta_mixed": 100,
    "multidomain_or_boundary_sensitive": 100,
    "disordered_or_flexible": 100,
}


def source_catalog_to_dict(
    sources: Sequence[FoldingBenchmarkSource] = REAL_FOLDING_BENCHMARK_SOURCES,
) -> list[dict[str, object]]:
    return [source.to_dict() for source in sources]


def source_ids(
    sources: Sequence[FoldingBenchmarkSource] = REAL_FOLDING_BENCHMARK_SOURCES,
) -> tuple[str, ...]:
    return tuple(source.source_id for source in sources)


def accepted_reference_prefixes(
    sources: Sequence[FoldingBenchmarkSource] = REAL_FOLDING_BENCHMARK_SOURCES,
) -> tuple[str, ...]:
    prefixes: list[str] = []
    for source in sources:
        prefixes.extend(source.accepted_prefixes)
    return tuple(sorted(set(prefixes)))


def validate_source_catalog(
    sources: Sequence[FoldingBenchmarkSource] = REAL_FOLDING_BENCHMARK_SOURCES,
) -> FoldingBenchmarkSourceReview:
    violations: list[str] = []
    seen: set[str] = set()
    for source in sources:
        if not source.source_id:
            violations.append("source_id_empty")
        if source.source_id in seen:
            violations.append(f"duplicate_source_id:{source.source_id}")
        seen.add(source.source_id)
        if not source.accepted_prefixes:
            violations.append(f"{source.source_id}.accepted_prefixes_empty")
        if not source.source_url.startswith("https://"):
            violations.append(f"{source.source_id}.source_url_not_https")
        if any(
            (
                source.clinical_use_allowed,
                source.drug_design_created,
                source.molecule_generated,
                source.protein_sequence_design_created,
                source.folding_solution_claim_created,
            )
        ):
            violations.append(f"{source.source_id}.boundary_opened")

    return FoldingBenchmarkSourceReview(
        valid=not violations,
        sources_reviewed=len(sources),
        source_ids=source_ids(sources),
        accepted_prefixes=accepted_reference_prefixes(sources),
        violations=tuple(violations),
    )
