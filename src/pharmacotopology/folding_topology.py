from __future__ import annotations

from dataclasses import asdict, dataclass
from math import log
from typing import Mapping, Sequence


FOLDING_TOPOLOGY_DIMENSIONS: tuple[str, ...] = (
    "sequence_complexity",
    "secondary_structure_balance",
    "contact_map_closure",
    "hydrophobic_core_closure",
    "loop_disorder_pressure",
    "domain_boundary_stability",
    "long_range_contact_order",
    "conformational_flexibility",
    "knot_or_entanglement_signature",
    "uncertainty_radius",
)

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC_AMINO_ACIDS = frozenset("AVILMFWYC")
AROMATIC_AMINO_ACIDS = frozenset("FWY")
CHARGED_AMINO_ACIDS = frozenset("DEKRH")
POLAR_AMINO_ACIDS = frozenset("STNQ")
DISORDER_PROMOTING_AMINO_ACIDS = frozenset("PGSQEKR")

CONTACT_PROXY_DIMENSIONS: tuple[str, ...] = (
    "contact_map_closure",
    "hydrophobic_core_closure",
    "long_range_contact_order",
    "knot_or_entanglement_signature",
)

BROAD_FOLD_CLASSES: tuple[str, ...] = (
    "alpha_rich",
    "beta_rich",
    "alpha_beta_mixed",
    "multidomain_boundary",
    "disordered_flexible",
)

INTERNAL_TO_BROAD_FOLD_CLASS: dict[str, str] = {
    "compact_folded": "alpha_rich",
    "long_range_contact_rich": "beta_rich",
    "flexible_or_disordered": "disordered_flexible",
    "boundary_sensitive": "multidomain_boundary",
    "entanglement_rich": "multidomain_boundary",
    "mixed_uncertain": "alpha_beta_mixed",
}

KNOWN_FOLD_CLASSES: tuple[str, ...] = BROAD_FOLD_CLASSES + tuple(
    INTERNAL_TO_BROAD_FOLD_CLASS
)

PLACEHOLDER_REFERENCE_SOURCE_MARKERS: tuple[str, ...] = (
    "placeholder",
    "example",
    "replace_me",
    "template",
)


@dataclass(frozen=True)
class FoldingTopologySignature:
    sequence_complexity: float
    secondary_structure_balance: float
    contact_map_closure: float
    hydrophobic_core_closure: float
    loop_disorder_pressure: float
    domain_boundary_stability: float
    long_range_contact_order: float
    conformational_flexibility: float
    knot_or_entanglement_signature: float
    uncertainty_radius: float


@dataclass(frozen=True)
class FoldingReferenceExample:
    protein_id: str
    sequence: str
    reference_structure_source: str
    reference_fold_class: str
    reference_topology_signature: FoldingTopologySignature
    source_database: str = ""
    source_accession: str = ""
    sequence_length: int = 0
    reference_label_source: str = ""
    is_external_reference: bool = False
    curation_notes: tuple[str, ...] = ()
    reference_topology_signature_kind: str = "class_prototype"


@dataclass(frozen=True)
class FoldingTopologyComparison:
    protein_id: str
    sequence_length: int
    source_database: str
    source_accession: str
    reference_structure_source: str
    reference_label_source: str
    predicted_topology_signature: FoldingTopologySignature
    reference_topology_signature: FoldingTopologySignature
    predicted_internal_fold_class: str
    reference_source_fold_class: str
    predicted_fold_class: str
    reference_fold_class: str
    contact_map_similarity: float
    fold_class_match: bool
    uncertainty_radius: float
    evidence_readiness: str
    failure_reason: str
    is_external_reference: bool = False
    curation_notes: tuple[str, ...] = ()
    reference_topology_signature_kind: str = "class_prototype"
    clinical_use_allowed: bool = False
    drug_design_created: bool = False
    molecule_generated: bool = False
    protein_sequence_design_created: bool = False
    folding_solution_claim_created: bool = False


DEFAULT_FOLDING_BENCHMARKS: tuple[FoldingReferenceExample, ...] = (
    FoldingReferenceExample(
        protein_id="alpha_compact_reference",
        sequence=(
            "MALAELAKQVAQLKDELAALQKELGAALEQAKQALEAAGKAVQLAEEG"
            "LAQLKQEVAALEKQLA"
        ),
        reference_structure_source="benchmark_placeholder:compact_alpha_summary",
        reference_fold_class="compact_folded",
        reference_topology_signature=FoldingTopologySignature(
            sequence_complexity=0.61,
            secondary_structure_balance=0.74,
            contact_map_closure=0.72,
            hydrophobic_core_closure=0.65,
            loop_disorder_pressure=0.18,
            domain_boundary_stability=0.76,
            long_range_contact_order=0.48,
            conformational_flexibility=0.28,
            knot_or_entanglement_signature=0.18,
            uncertainty_radius=0.26,
        ),
    ),
    FoldingReferenceExample(
        protein_id="long_range_beta_reference",
        sequence=(
            "TVVTIENEGKTVTFTDGEVTVNGEKTVTFEGQAVTVKDGEVTVNGETVT"
            "FEGQA"
        ),
        reference_structure_source="benchmark_placeholder:long_range_contact_summary",
        reference_fold_class="long_range_contact_rich",
        reference_topology_signature=FoldingTopologySignature(
            sequence_complexity=0.56,
            secondary_structure_balance=0.69,
            contact_map_closure=0.64,
            hydrophobic_core_closure=0.50,
            loop_disorder_pressure=0.26,
            domain_boundary_stability=0.62,
            long_range_contact_order=0.62,
            conformational_flexibility=0.34,
            knot_or_entanglement_signature=0.30,
            uncertainty_radius=0.30,
        ),
    ),
    FoldingReferenceExample(
        protein_id="flexible_loop_reference",
        sequence=(
            "MGGSGQPGSQEKRGGSGQPGSQEKRGGSGQPGSQEKRGGSGQPGSQEKR"
        ),
        reference_structure_source="benchmark_placeholder:flexible_loop_summary",
        reference_fold_class="flexible_or_disordered",
        reference_topology_signature=FoldingTopologySignature(
            sequence_complexity=0.38,
            secondary_structure_balance=0.42,
            contact_map_closure=0.22,
            hydrophobic_core_closure=0.16,
            loop_disorder_pressure=0.78,
            domain_boundary_stability=0.34,
            long_range_contact_order=0.18,
            conformational_flexibility=0.82,
            knot_or_entanglement_signature=0.06,
            uncertainty_radius=0.58,
        ),
    ),
    FoldingReferenceExample(
        protein_id="boundary_sensitive_reference",
        sequence=(
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFY"
            "TPKTQLEKRGIVEQCCTSICSLYQLENYCN"
        ),
        reference_structure_source="benchmark_placeholder:domain_boundary_summary",
        reference_fold_class="boundary_sensitive",
        reference_topology_signature=FoldingTopologySignature(
            sequence_complexity=0.70,
            secondary_structure_balance=0.58,
            contact_map_closure=0.54,
            hydrophobic_core_closure=0.56,
            loop_disorder_pressure=0.38,
            domain_boundary_stability=0.44,
            long_range_contact_order=0.40,
            conformational_flexibility=0.46,
            knot_or_entanglement_signature=0.16,
            uncertainty_radius=0.40,
        ),
    ),
)


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def _rounded(value: float) -> float:
    return round(_clamp(value), 6)


def normalize_sequence(sequence: str) -> str:
    cleaned = "".join(sequence.upper().split()).replace("-", "")
    if not cleaned:
        raise ValueError("Protein sequence must not be empty.")
    invalid = sorted(set(cleaned).difference(AMINO_ACIDS))
    if invalid:
        raise ValueError(
            "Protein sequence contains unsupported amino-acid symbols: "
            + ", ".join(invalid)
        )
    return cleaned


def _fraction(sequence: str, alphabet: frozenset[str]) -> float:
    if not sequence:
        return 0.0
    return sum(1 for residue in sequence if residue in alphabet) / len(sequence)


def sequence_features(sequence: str) -> dict[str, float]:
    normalized = normalize_sequence(sequence)
    counts = {
        residue: normalized.count(residue)
        for residue in sorted(set(normalized))
    }
    entropy = 0.0
    for count in counts.values():
        probability = count / len(normalized)
        entropy -= probability * (log(probability) / log(20))

    length_norm = _clamp((len(normalized) - 30) / 170)
    hydrophobic = _fraction(normalized, HYDROPHOBIC_AMINO_ACIDS)
    aromatic = _fraction(normalized, AROMATIC_AMINO_ACIDS)
    charged = _fraction(normalized, CHARGED_AMINO_ACIDS)
    polar = _fraction(normalized, POLAR_AMINO_ACIDS)
    disorder = _fraction(normalized, DISORDER_PROMOTING_AMINO_ACIDS)
    proline_glycine = _fraction(normalized, frozenset("PG"))
    cysteine = _fraction(normalized, frozenset("C"))

    return {
        "length": float(len(normalized)),
        "length_norm": _rounded(length_norm),
        "sequence_complexity": _rounded(entropy),
        "hydrophobic_fraction": _rounded(hydrophobic),
        "aromatic_fraction": _rounded(aromatic),
        "charged_fraction": _rounded(charged),
        "polar_fraction": _rounded(polar),
        "disorder_promoting_fraction": _rounded(disorder),
        "proline_glycine_fraction": _rounded(proline_glycine),
        "cysteine_fraction": _rounded(cysteine),
    }


def predict_topology_signature(sequence: str) -> FoldingTopologySignature:
    features = sequence_features(sequence)
    complexity = features["sequence_complexity"]
    length_norm = features["length_norm"]
    hydrophobic = features["hydrophobic_fraction"]
    aromatic = features["aromatic_fraction"]
    charged = features["charged_fraction"]
    polar = features["polar_fraction"]
    disorder = features["disorder_promoting_fraction"]
    proline_glycine = features["proline_glycine_fraction"]
    cysteine = features["cysteine_fraction"]

    loop_disorder_pressure = (
        disorder * 0.80
        + proline_glycine * 0.45
        + charged * 0.18
        + (1.0 - complexity) * 0.12
    )
    hydrophobic_core_closure = (
        hydrophobic * 1.28
        + aromatic * 0.30
        + cysteine * 0.20
        + length_norm * 0.10
        - proline_glycine * 0.50
        - charged * 0.12
    )
    secondary_structure_balance = 1.0 - abs(
        (hydrophobic + polar + aromatic) - (disorder + charged)
    )
    contact_map_closure = (
        hydrophobic_core_closure * 0.62
        + complexity * 0.20
        + length_norm * 0.20
        - loop_disorder_pressure * 0.16
    )
    domain_boundary_stability = (
        0.82
        - loop_disorder_pressure * 0.45
        - charged * 0.20
        + cysteine * 0.18
        + complexity * 0.12
    )
    long_range_contact_order = (
        contact_map_closure * 0.46
        + length_norm * 0.34
        + hydrophobic_core_closure * 0.18
        - loop_disorder_pressure * 0.10
    )
    conformational_flexibility = (
        loop_disorder_pressure * 0.68
        + charged * 0.24
        + polar * 0.18
        + (1.0 - hydrophobic_core_closure) * 0.10
    )
    knot_or_entanglement_signature = (
        long_range_contact_order * 0.42
        + length_norm * 0.28
        + hydrophobic_core_closure * 0.12
        - conformational_flexibility * 0.20
    )
    uncertainty_radius = (
        0.22
        + loop_disorder_pressure * 0.26
        + (1.0 - complexity) * 0.18
        + (1.0 - length_norm) * 0.08
    )

    return FoldingTopologySignature(
        sequence_complexity=_rounded(complexity),
        secondary_structure_balance=_rounded(secondary_structure_balance),
        contact_map_closure=_rounded(contact_map_closure),
        hydrophobic_core_closure=_rounded(hydrophobic_core_closure),
        loop_disorder_pressure=_rounded(loop_disorder_pressure),
        domain_boundary_stability=_rounded(domain_boundary_stability),
        long_range_contact_order=_rounded(long_range_contact_order),
        conformational_flexibility=_rounded(conformational_flexibility),
        knot_or_entanglement_signature=_rounded(knot_or_entanglement_signature),
        uncertainty_radius=_rounded(uncertainty_radius),
    )


def signature_to_dict(signature: FoldingTopologySignature) -> dict[str, float]:
    values = asdict(signature)
    return {
        dimension: float(values[dimension])
        for dimension in FOLDING_TOPOLOGY_DIMENSIONS
    }


def classify_fold_topology(signature: FoldingTopologySignature) -> str:
    values = signature_to_dict(signature)
    if (
        values["loop_disorder_pressure"] >= 0.56
        and values["conformational_flexibility"] >= 0.50
    ):
        return "flexible_or_disordered"
    if values["domain_boundary_stability"] <= 0.48:
        return "boundary_sensitive"
    if (
        values["knot_or_entanglement_signature"] >= 0.34
        and values["long_range_contact_order"] >= 0.50
    ):
        return "entanglement_rich"
    if values["long_range_contact_order"] >= 0.52:
        return "long_range_contact_rich"
    if (
        values["contact_map_closure"] >= 0.50
        and values["hydrophobic_core_closure"] >= 0.45
    ):
        return "compact_folded"
    return "mixed_uncertain"


def normalize_fold_class(fold_class: str) -> str:
    if fold_class in BROAD_FOLD_CLASSES:
        return fold_class
    return INTERNAL_TO_BROAD_FOLD_CLASS.get(fold_class, "alpha_beta_mixed")


def classify_broad_fold_topology(signature: FoldingTopologySignature) -> str:
    values = signature_to_dict(signature)
    if (
        values["loop_disorder_pressure"] >= 0.56
        and values["conformational_flexibility"] >= 0.50
    ):
        return "disordered_flexible"
    if (
        values["domain_boundary_stability"] <= 0.48
        or values["knot_or_entanglement_signature"] >= 0.34
    ):
        return "multidomain_boundary"
    if (
        values["long_range_contact_order"] >= 0.52
        and values["hydrophobic_core_closure"] < 0.56
    ):
        return "beta_rich"
    if (
        values["contact_map_closure"] >= 0.50
        and values["hydrophobic_core_closure"] >= 0.45
        and values["conformational_flexibility"] < 0.45
    ):
        return "alpha_rich"
    return "alpha_beta_mixed"


def contact_map_proxy_similarity(
    predicted: FoldingTopologySignature,
    reference: FoldingTopologySignature,
    dimensions: Sequence[str] = CONTACT_PROXY_DIMENSIONS,
) -> float:
    predicted_values = signature_to_dict(predicted)
    reference_values = signature_to_dict(reference)
    distance = sum(
        abs(predicted_values[dimension] - reference_values[dimension])
        for dimension in dimensions
    ) / len(dimensions)
    return _rounded(1.0 - distance)


def evidence_readiness_label(
    reference_structure_source: str,
    contact_map_similarity: float,
    fold_class_match: bool,
    uncertainty_radius: float,
) -> str:
    if reference_source_is_placeholder(reference_structure_source):
        return "benchmark_shell_only"
    return "external_benchmark_attached"


def reference_source_is_placeholder(source: str) -> bool:
    normalized = source.strip().lower()
    return any(marker in normalized for marker in PLACEHOLDER_REFERENCE_SOURCE_MARKERS)


def compare_to_reference(reference: FoldingReferenceExample) -> FoldingTopologyComparison:
    normalized = normalize_sequence(reference.sequence)
    predicted = predict_topology_signature(normalized)
    predicted_internal_class = classify_fold_topology(predicted)
    predicted_class = classify_broad_fold_topology(predicted)
    reference_class = normalize_fold_class(reference.reference_fold_class)
    contact_similarity = contact_map_proxy_similarity(
        predicted,
        reference.reference_topology_signature,
    )
    fold_class_match = predicted_class == reference_class
    uncertainty = _rounded(
        max(
            predicted.uncertainty_radius,
            reference.reference_topology_signature.uncertainty_radius,
        )
    )
    readiness = evidence_readiness_label(
        reference.reference_structure_source,
        contact_similarity,
        fold_class_match,
        uncertainty,
    )
    failure_reason = ""
    if readiness == "benchmark_shell_only":
        failure_reason = "external_structure_benchmark_not_attached"
    elif not fold_class_match:
        failure_reason = "fold_class_proxy_mismatch"
    elif contact_similarity < 0.70:
        failure_reason = "low_topology_similarity"

    return FoldingTopologyComparison(
        protein_id=reference.protein_id,
        sequence_length=len(normalized),
        source_database=reference.source_database,
        source_accession=reference.source_accession,
        reference_structure_source=reference.reference_structure_source,
        reference_label_source=reference.reference_label_source,
        predicted_topology_signature=predicted,
        reference_topology_signature=reference.reference_topology_signature,
        predicted_internal_fold_class=predicted_internal_class,
        reference_source_fold_class=reference.reference_fold_class,
        predicted_fold_class=predicted_class,
        reference_fold_class=reference_class,
        contact_map_similarity=contact_similarity,
        fold_class_match=fold_class_match,
        uncertainty_radius=uncertainty,
        evidence_readiness=readiness,
        failure_reason=failure_reason,
        is_external_reference=reference.is_external_reference,
        curation_notes=reference.curation_notes,
        reference_topology_signature_kind=reference.reference_topology_signature_kind,
    )


def run_folding_topology_benchmark(
    references: Sequence[FoldingReferenceExample] = DEFAULT_FOLDING_BENCHMARKS,
) -> tuple[FoldingTopologyComparison, ...]:
    return tuple(compare_to_reference(reference) for reference in references)


def comparison_to_dict(
    comparison: FoldingTopologyComparison,
) -> dict[str, object]:
    data = asdict(comparison)
    data["predicted_topology_signature"] = signature_to_dict(
        comparison.predicted_topology_signature
    )
    data["reference_topology_signature"] = signature_to_dict(
        comparison.reference_topology_signature
    )
    return data


def reference_from_mapping(data: Mapping[str, object]) -> FoldingReferenceExample:
    signature_data = data.get("reference_topology_signature")
    if not isinstance(signature_data, Mapping):
        raise ValueError("reference_topology_signature must be a mapping.")
    signature = FoldingTopologySignature(
        **{
            dimension: float(signature_data[dimension])
            for dimension in FOLDING_TOPOLOGY_DIMENSIONS
        }
    )
    return FoldingReferenceExample(
        protein_id=str(data["protein_id"]),
        sequence=str(data["sequence"]),
        reference_structure_source=str(data["reference_structure_source"]),
        reference_fold_class=str(data["reference_fold_class"]),
        reference_topology_signature=signature,
        source_database=str(data.get("source_database", "")),
        source_accession=str(data.get("source_accession", "")),
        sequence_length=int(data.get("sequence_length", 0) or 0),
        reference_label_source=str(data.get("reference_label_source", "")),
        is_external_reference=bool(data.get("is_external_reference", False)),
        curation_notes=tuple(
            str(note) for note in data.get("curation_notes", ()) or ()
        ),
        reference_topology_signature_kind=str(
            data.get("reference_topology_signature_kind", "class_prototype")
        ),
    )
