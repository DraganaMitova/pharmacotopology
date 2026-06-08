from __future__ import annotations

from dataclasses import dataclass, asdict
import os
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
from math import sqrt
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_native_contact_eval import (
    ContactPair,
    evaluate_contact_prediction,
    normalized_contact_pairs,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import RealCoordinateVisualRow
from pharmacotopology.folding_structure_benchmark import AA3_TO_1

TEMPLATE_DOCKING_KIND = "template_docking_v0"
TEMPLATE_DOCKING_SOURCE_KIND = "template_docking_interface_v0"
TEMPLATE_DOCKING_SOURCE_FAMILY = "template_docking"
TEMPLATE_DOCKING_KIND_CONTACT_REPORT = "template_docking_contact_vote_v0"
TEMPLATE_DOCKING_THRESHOLD_KIND = "template_docking_v0_internal_gap_vote"

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
DEFAULT_CONTACT_CUTOFF_ANGSTROM = 8.0
DEFAULT_ALIGNMENT_MATCH_SCORE = 1.0
DEFAULT_ALIGNMENT_MISMATCH_SCORE = -1.0
DEFAULT_ALIGNMENT_GAP_PENALTY = -1.0
DEFAULT_MIN_SEQUENCE_SEPARATION = 3
DEFAULT_DOMAIN_SPLIT_MIN = 20
DEFAULT_DOMAINS_MIN_CONTACTS = 12
CHEMICAL_HYDROPHOBIC = set("AVILMFWY")
CHEMICAL_POSITIVE = set("KRH")
CHEMICAL_NEGATIVE = set("DE")
CHEMICAL_POLAR = set("STNQ")
CHEMICAL_THRESHOLD_DEFAULT = 0.6
CHEMICAL_THRESHOLD_AUTO = None


@dataclass(frozen=True)
class TemplateResidue:
    index: int
    aa: str
    x: float
    y: float
    z: float
    residue_number: int
    insertion_code: str
    residue_number_text: str


@dataclass(frozen=True)
class TemplateStructure:
    template_id: str
    sequence: str
    residues: tuple[TemplateResidue, ...]
    source_accession: str
    source_id: str
    chain_id: str | None = None
    source_path: str | None = None
    source_kind: str = "benchmark_row"


@dataclass(frozen=True)
class TemplateInterfaceResult:
    template_id: str
    source_accession: str
    source_id: str
    source_kind: str
    domain_split: int
    interface_contact_count: int
    mapped_contact_count: int
    chem_kept_contact_count: int
    alignment_score: float
    alignment_identity: float
    mapped_pairs: tuple[ContactPair, ...]
    voted_pairs: tuple[ContactPair, ...]
    mapping_mode: str = "full_alignment"
    source_path: str | None = None


@dataclass(frozen=True)
class TemplateDockingReport:
    kind: str
    source_mode: str
    target_source_accession: str
    target_row_id: str
    target_sequence_length: int
    template_count: int
    templates_used: int
    templates_with_mapped_pairs: int
    winner_template_id: str
    winner_template_source: str
    winner_vote_count: int
    min_template_support: int
    predicted_pair_count: int
    voted_pair_count: int
    winner_pair_count: int
    long_range_recall: float
    native_contact_precision: float
    native_contact_recall: float
    contact_map_f1: float
    false_contact_rate: float
    target_accession: str
    certificate_path: str | None

    template_rows: tuple[dict[str, object], ...]
    voted_pairs: tuple[dict[str, object], ...]
    winner_pairs: tuple[dict[str, object], ...]
    predicted_pairs: tuple[dict[str, object], ...]
    metric_after_native_audit: Mapping[str, object]
    long_range_threshold: int

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["metric_after_native_audit"] = dict(self.metric_after_native_audit)
        payload["template_rows"] = list(self.template_rows)
        payload["voted_pairs"] = list(self.voted_pairs)
        payload["winner_pairs"] = list(self.winner_pairs)
        payload["predicted_pairs"] = list(self.predicted_pairs)
        return payload


@dataclass(frozen=True)
class TemplateConstraintRecord:
    template_accession: str
    template_source_id: str
    confidence: float
    support_count: int
    total_template_pairs: int
    pair: ContactPair


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _distance(left: TemplateResidue, right: TemplateResidue) -> float:
    return sqrt(
        (left.x - right.x) ** 2 + (left.y - right.y) ** 2 + (left.z - right.z) ** 2
    )


def chemical_score(aa_i: str, aa_j: str) -> float:
    if aa_i not in AA_ALPHABET or aa_j not in AA_ALPHABET:
        return 0.1
    if aa_i in CHEMICAL_HYDROPHOBIC and aa_j in CHEMICAL_HYDROPHOBIC:
        return 1.0
    if (aa_i in CHEMICAL_POSITIVE and aa_j in CHEMICAL_NEGATIVE) or (
        aa_i in CHEMICAL_NEGATIVE and aa_j in CHEMICAL_POSITIVE
    ):
        return 0.9
    if aa_i in CHEMICAL_POLAR and aa_j in CHEMICAL_POLAR:
        return 0.5
    if (aa_i in CHEMICAL_POLAR and aa_j in CHEMICAL_HYDROPHOBIC) or (
        aa_i in CHEMICAL_HYDROPHOBIC and aa_j in CHEMICAL_POLAR
    ):
        return 0.3
    return 0.1


def _aa3_to_aa1(res_name: str) -> str:
    return AA3_TO_1.get(res_name.upper(), "X")


def parse_benchmark_templates(rows: Sequence[RealCoordinateVisualRow], target_accession: str) -> tuple[TemplateStructure, ...]:
    templates: list[TemplateStructure] = []
    for row in rows:
        if row.source_accession == target_accession:
            continue
        residues: list[TemplateResidue] = []
        for point in row.coordinate_points:
            aa = row.sequence[point.sequence_index - 1]
            residues.append(
                TemplateResidue(
                    index=point.sequence_index,
                    aa=aa,
                    x=point.x,
                    y=point.y,
                    z=point.z,
                    residue_number=point.residue_number,
                    insertion_code=point.insertion_code,
                    residue_number_text=f"{point.residue_number}{point.insertion_code}",
                )
            )
        templates.append(
            TemplateStructure(
                template_id=row.source_accession,
                sequence=row.sequence,
                residues=tuple(residues),
                source_accession=row.source_accession,
                source_id=row.source_id,
                chain_id="A",
                source_path=None,
                source_kind="benchmark_row",
            )
        )
    return tuple(templates)


def parse_pdb_template(path: Path, chain_id: str | None = None) -> tuple[TemplateStructure, ...]:
    if not path.exists():
        raise FileNotFoundError(f"PDB template path does not exist: {path}")

    residues_by_chain: dict[str, list[TemplateResidue]] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        chain = line[21].strip() or "A"
        try:
            residue_number = int(line[22:26].strip())
            insertion_code = line[26].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue
        residues_by_chain.setdefault(chain, [])
        residues = residues_by_chain[chain]
        aa = _aa3_to_aa1(line[17:20].strip())
        residues.append(
            TemplateResidue(
                index=len(residues) + 1,
                aa=aa,
                x=x,
                y=y,
                z=z,
                residue_number=residue_number,
                insertion_code=insertion_code,
                residue_number_text=f"{residue_number}{insertion_code}",
            )
        )

    if chain_id is None:
        if not residues_by_chain:
            return tuple()
        # Keep the chain with most residues as the best default.
        chain_id = max(
            residues_by_chain.items(), key=lambda item: (len(item[1]), item[0])
        )[0]

    residues = residues_by_chain.get(chain_id)
    if not residues:
        return tuple()

    sequence = "".join(residue.aa for residue in residues)
    template_id = f"{path.name}:{chain_id}"
    return (
        TemplateStructure(
            template_id=template_id,
            sequence=sequence,
            residues=tuple(residues),
            source_accession=template_id,
            source_id=template_id,
            chain_id=chain_id,
            source_path=str(path),
            source_kind="pdb_file",
        ),
    )


def collect_pdb_templates(
    template_paths: Sequence[Path],
    chain_id: str | None = None,
    max_template_count: int | None = None,
) -> tuple[TemplateStructure, ...]:
    templates: list[TemplateStructure] = []
    for path in template_paths:
        for template in parse_pdb_template(path, chain_id=chain_id):
            templates.append(template)
            if max_template_count and len(templates) >= max_template_count:
                return tuple(templates)
    return tuple(templates)


def _collect_contacts(
    residues: Sequence[TemplateResidue],
    *,
    contact_cutoff_angstrom: float = DEFAULT_CONTACT_CUTOFF_ANGSTROM,
    minimum_sequence_separation: int = DEFAULT_MIN_SEQUENCE_SEPARATION,
) -> tuple[ContactPair, ...]:
    if not residues:
        return tuple()
    contacts: list[ContactPair] = []
    cutoff_sq = contact_cutoff_angstrom * contact_cutoff_angstrom
    ordered = residues
    for left_position, left in enumerate(ordered[:-1], start=1):
        left_index = left.index
        for right in ordered[left_position:]:
            right_index = right.index
            if right_index - left_index < minimum_sequence_separation:
                continue
            dx = left.x - right.x
            dy = left.y - right.y
            dz = left.z - right.z
            if dx * dx + dy * dy + dz * dz <= cutoff_sq:
                contacts.append((left.index, right_index))
    return normalized_contact_pairs(contacts)


def _best_domain_split(
    contacts: Sequence[ContactPair],
    *,
    sequence_length: int,
    minimum_domain_size: int = DEFAULT_DOMAIN_SPLIT_MIN,
) -> int:
    if sequence_length <= minimum_domain_size * 2:
        return max(1, sequence_length // 2)

    diff = [0] * (sequence_length + 2)
    for left, right in contacts:
        if right <= left:
            left, right = right, left
        start = left
        end = right
        if end - start < 2:
            continue
        diff[start] += 1
        diff[end] -= 1
    running = 0
    best_split = sequence_length // 2
    best_count = -1
    for split in range(minimum_domain_size, sequence_length - minimum_domain_size + 1):
        running += diff[split]
        if running > best_count:
            best_count = running
            best_split = split
    if best_count <= 0:
        return sequence_length // 2
    return best_split


def _interface_contacts(
    contacts: Sequence[ContactPair],
    split: int,
) -> tuple[ContactPair, ...]:
    result = []
    for left, right in contacts:
        if left <= split < right:
            result.append((left, right))
    return normalized_contact_pairs(result)


def _pair_score_left_align(
    template: str,
    target: str,
    *,
    match: float = DEFAULT_ALIGNMENT_MATCH_SCORE,
    mismatch: float = DEFAULT_ALIGNMENT_MISMATCH_SCORE,
    gap: float = DEFAULT_ALIGNMENT_GAP_PENALTY,
) -> tuple[dict[int, int], float, float]:
    n = len(template)
    m = len(target)
    if n == 0:
        return ({}, 0.0, 0.0)

    scores = [[0.0] * (m + 1) for _ in range(n + 1)]
    moves: list[list[int]] = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        scores[i][0] = scores[i - 1][0] + gap
        moves[i][0] = 1
    for j in range(1, m + 1):
        scores[0][j] = scores[0][j - 1] + gap
        moves[0][j] = 2

    for i, aa_t in enumerate(template, start=1):
        for j, aa_s in enumerate(target, start=1):
            match_score = match if aa_t == aa_s else mismatch
            diag = scores[i - 1][j - 1] + match_score
            up = scores[i - 1][j] + gap
            left = scores[i][j - 1] + gap
            best = max(diag, up, left)
            scores[i][j] = best
            if best == diag:
                moves[i][j] = 0
            elif best == up:
                moves[i][j] = 1
            else:
                moves[i][j] = 2

    aligned_template: list[str] = []
    aligned_target: list[str] = []
    i = n
    j = m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and moves[i][j] == 0:
            aligned_template.append(template[i - 1])
            aligned_target.append(target[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and moves[i][j] == 1:
            aligned_template.append(template[i - 1])
            aligned_target.append("-")
            i -= 1
        else:
            aligned_template.append("-")
            aligned_target.append(target[j - 1])
            j -= 1

    aligned_template.reverse()
    aligned_target.reverse()

    template_to_target: dict[int, int] = {}
    t_pos = 0
    s_pos = 0
    exact_matches = 0
    aligned_positions = 0
    for a_t, a_s in zip(aligned_template, aligned_target):
        if a_t != "-":
            t_pos += 1
        if a_s != "-":
            s_pos += 1
        if a_t != "-" and a_s != "-":
            aligned_positions += 1
            if a_t == a_s:
                exact_matches += 1
            template_to_target[t_pos] = s_pos

    identity = 0.0
    if aligned_positions > 0:
        identity = exact_matches / aligned_positions
    return template_to_target, scores[n][m], _rounded(identity)


def _domain_aware_template_mapping(
    template: str,
    template_split: int,
    target: str,
    interface_pairs: Sequence[ContactPair],
    *,
    minimum_sequence_separation: int,
    chemical_threshold: float,
    minimum_domain_size: int = DEFAULT_DOMAIN_SPLIT_MIN,
    domain_scan_step: int = 1,
) -> tuple[tuple[ContactPair, ...], float]:
    if not template or not target:
        return tuple(), 0.0

    left_template = template[:template_split]
    right_template = template[template_split:]
    target_length = len(target)
    if not left_template or not right_template:
        return tuple(), 0.0

    step = max(1, domain_scan_step)
    last_split = max(minimum_domain_size + 1, target_length - minimum_domain_size + 1)
    target_window = list(range(minimum_domain_size, last_split, step))
    if not target_window:
        target_window = [minimum_domain_size]
    elif target_window[-1] != last_split - 1:
        target_window.append(last_split - 1)

    best_mapped: tuple[ContactPair, ...] = tuple()
    best_score = -1.0
    for target_split in target_window:
        left_target = target[:target_split]
        right_target = target[target_split:]
        if not left_target or not right_target:
            continue

        left_map, left_score, left_identity = _pair_score_left_align(
            left_template, left_target
        )
        right_map, right_score, right_identity = _pair_score_left_align(
            right_template, right_target
        )

        mapped_pairs: list[ContactPair] = []
        for left_t, right_t in interface_pairs:
            left_source = left_t
            right_source = right_t
            mapped_left: int | None = None
            mapped_right: int | None = None

            if left_source <= template_split:
                mapped_left = left_map.get(left_source)
            else:
                mapped_left = right_map.get(left_source - template_split)
                if mapped_left is not None:
                    mapped_left += target_split

            if right_source <= template_split:
                mapped_right = left_map.get(right_source)
            else:
                mapped_right = right_map.get(right_source - template_split)
                if mapped_right is not None:
                    mapped_right += target_split

            if mapped_left is None or mapped_right is None:
                continue

            i = min(mapped_left, mapped_right)
            j = max(mapped_left, mapped_right)
            if j - i < minimum_sequence_separation:
                continue
            if i < 1 or j > target_length:
                continue
            if chemical_score(target[i - 1], target[j - 1]) < chemical_threshold:
                continue

            mapped_pairs.append((i, j))

        mapped_pairs_t = normalized_contact_pairs(mapped_pairs)
        if not mapped_pairs_t:
            continue

        score = (
            1000.0 * len(mapped_pairs_t)
            + 50.0 * (left_identity + right_identity)
            + left_score
            + right_score
        )
        if score > best_score:
            best_score = score
            best_mapped = tuple(mapped_pairs_t)

    return best_mapped, best_score


def map_template_interface_to_target(
    template: TemplateStructure,
    target_sequence: str,
    *,
    contact_cutoff_angstrom: float,
    minimum_sequence_separation: int,
    chemical_threshold: float | None = CHEMICAL_THRESHOLD_AUTO,
    domain_scan_step: int = 1,
) -> tuple[tuple[ContactPair, ...], TemplateInterfaceResult]:
    residues = template.residues
    if not residues:
        return (), TemplateInterfaceResult(
            template_id=template.template_id,
            source_accession=template.source_accession,
            source_id=template.source_id,
            source_kind=template.source_kind,
            domain_split=0,
            interface_contact_count=0,
            mapped_contact_count=0,
            chem_kept_contact_count=0,
            mapping_mode="full_alignment",
            alignment_score=0.0,
            alignment_identity=0.0,
            mapped_pairs=tuple(),
            voted_pairs=tuple(),
            source_path=template.source_path,
        )

    contacts = _collect_contacts(
        residues,
        contact_cutoff_angstrom=contact_cutoff_angstrom,
        minimum_sequence_separation=minimum_sequence_separation,
    )
    if not contacts:
        return (), TemplateInterfaceResult(
            template_id=template.template_id,
            source_accession=template.source_accession,
            source_id=template.source_id,
            source_kind=template.source_kind,
            domain_split=0,
            interface_contact_count=0,
            mapped_contact_count=0,
            chem_kept_contact_count=0,
            mapping_mode="full_alignment",
            alignment_score=0.0,
            alignment_identity=0.0,
            mapped_pairs=tuple(),
            voted_pairs=tuple(),
            source_path=template.source_path,
        )

    split = _best_domain_split(contacts, sequence_length=len(template.sequence))
    interface = _interface_contacts(contacts, split=split)
    if len(interface) < DEFAULT_DOMAINS_MIN_CONTACTS:
        split = max(1, len(template.sequence) // 2)
        interface = _interface_contacts(contacts, split=split)

    if not interface:
        return (), TemplateInterfaceResult(
            template_id=template.template_id,
            source_accession=template.source_accession,
            source_id=template.source_id,
            source_kind=template.source_kind,
            domain_split=split,
            interface_contact_count=len(contacts),
            mapped_contact_count=0,
            chem_kept_contact_count=0,
            mapping_mode="full_alignment",
            alignment_score=0.0,
            alignment_identity=0.0,
            mapped_pairs=tuple(),
            voted_pairs=tuple(),
            source_path=template.source_path,
        )

    mapping, alignment_score, alignment_identity = _pair_score_left_align(
        template.sequence,
        target_sequence,
    )
    candidate_pairs: list[tuple[ContactPair, float]] = []
    for left_t, right_t in interface:
        left_s = mapping.get(left_t)
        right_s = mapping.get(right_t)
        if left_s is None or right_s is None:
            continue
        i = min(left_s, right_s)
        j = max(left_s, right_s)
        if j - i < minimum_sequence_separation:
            continue
        if i < 1 or j > len(target_sequence):
            continue
        chem = chemical_score(target_sequence[i - 1], target_sequence[j - 1])
        candidate_pairs.append(((i, j), chem))

    candidate_chem_scores = [score for _, score in candidate_pairs]
    if chemical_threshold is None:
        resolved_chemical_threshold = _internal_gap_threshold(candidate_chem_scores)
    else:
        resolved_chemical_threshold = _rounded(chemical_threshold)

    mapped_pairs = [
        pair for pair, score in candidate_pairs if score >= resolved_chemical_threshold
    ]
    mapped_pairs_t = normalized_contact_pairs(mapped_pairs)
    voted_pairs = tuple(sorted(mapped_pairs_t))

    domain_mapped_pairs, domain_score = _domain_aware_template_mapping(
        template.sequence,
        template_split=split,
        target=target_sequence,
        interface_pairs=interface,
        minimum_sequence_separation=minimum_sequence_separation,
        chemical_threshold=resolved_chemical_threshold,
        domain_scan_step=domain_scan_step,
    )

    mapping_mode = "full_alignment"
    chosen_pairs = voted_pairs
    if len(domain_mapped_pairs) > len(voted_pairs):
        chosen_pairs = tuple(domain_mapped_pairs)
        mapping_mode = "domain_alignment"
    elif len(domain_mapped_pairs) == len(voted_pairs) and domain_score > 0 and alignment_score < domain_score:
        chosen_pairs = tuple(domain_mapped_pairs)
        mapping_mode = "domain_alignment"

    chosen_pairs = normalized_contact_pairs(chosen_pairs)
    aligned_score = alignment_score
    if mapping_mode == "domain_alignment" and domain_score > alignment_score:
        aligned_score = domain_score

    return chosen_pairs, TemplateInterfaceResult(
        template_id=template.template_id,
        source_accession=template.source_accession,
        source_id=template.source_id,
        source_kind=template.source_kind,
        domain_split=split,
        interface_contact_count=len(interface),
        mapped_contact_count=len(chosen_pairs),
        chem_kept_contact_count=len(chosen_pairs),
        mapping_mode=mapping_mode,
        alignment_score=aligned_score,
        alignment_identity=alignment_identity,
        mapped_pairs=tuple(chosen_pairs),
        voted_pairs=tuple(chosen_pairs),
        source_path=template.source_path,
    )


def _internal_gap_threshold(scores: Sequence[float]) -> float:
    if not scores:
        return 0.0
    sorted_scores = sorted(scores)
    max_gap = -1.0
    threshold = sorted_scores[0]
    for left, right in zip(sorted_scores[:-1], sorted_scores[1:]):
        gap = right - left
        if gap > max_gap:
            max_gap = gap
            threshold = left
    return _rounded(threshold + max_gap)


def rank_template_pairs(
    per_template: Mapping[str, tuple[ContactPair, ...]],
    template_weights: Mapping[str, float],
    *,
    long_range_threshold: int = 24,
) -> tuple[tuple[ContactPair, float, int], ...]:
    votes: dict[ContactPair, list[float]] = {}
    for template_id, pairs in per_template.items():
        weight = max(0.0, template_weights.get(template_id, 1.0))
        for pair in pairs:
            votes.setdefault(pair, []).append(weight)

    scored: list[tuple[ContactPair, float, int]] = []
    for pair, weights in votes.items():
        score = sum(weights)
        scored.append((pair, _rounded(score), len(weights)))

    if not scored:
        return tuple()

    # Internal-gap threshold on raw scores is a self-deciding gate.
    score_values = [item[1] for item in scored]
    threshold = _internal_gap_threshold(score_values)
    accepted = [item for item in scored if item[1] >= threshold]
    if not accepted:
        accepted = scored

    def _long_range(pair: ContactPair) -> int:
        return 1 if pair[1] - pair[0] >= long_range_threshold else 0

    accepted.sort(key=lambda item: (-(item[1]), -_long_range(item[0]), item[0][0], item[0][1]))
    return tuple(accepted)


def _build_constraint_payload(
    target_row: RealCoordinateVisualRow,
    pairs: Sequence[ContactPair],
    *,
    pair_support: Mapping[ContactPair, int],
    pair_templates: Mapping[ContactPair, tuple[str, ...]],
    template_count: int,
) -> dict[str, object]:
    constraints: list[dict[str, object]] = []
    for i, j in pairs:
        support = pair_support.get((i, j), 1)
        templates = pair_templates.get((i, j), tuple())
        confidence = _rounded(0.08 + (0.18 * min(5, support)))
        constraints.append(
            CouplingConstraint(
                row_id=target_row.row_id,
                source_accession=target_row.source_accession,
                constraint_id=f"template_dock_{target_row.source_accession.replace(':', '_')}_{i}_{j}",
                i=i,
                j=j,
                sequence_separation=j - i,
                normalized_separation=_rounded((j - i) / max(target_row.sequence_length, 1)),
                confidence=confidence,
                constraint_class="template_docking_interface_contact",
                source_kind=TEMPLATE_DOCKING_SOURCE_KIND,
                coordinate_truth_used_to_build_constraint=False,
                native_truth_used_before_coupling_selection=False,
                structure_model_used=False,
                raw_sequence_exposed=False,
                raw_score=confidence,
                apc_corrected_score=confidence,
                rank=0,
                rank_fraction=0.0,
                msa_source_kind=TEMPLATE_DOCKING_KIND,
                msa_sha256="",
                msa_depth=0,
                effective_sequence_count=0.0,
                effective_sequence_count_over_length=0.0,
                target_coverage=0.0,
                focus_sequence_mapping_confidence=min(1.0, support / max(1, template_count)),
            ).to_dict()
        )
    return {
        "layer_kind": "template_docking",
        "constraint_kind": "safe_residue_pair_coupling_constraint_v1",
        "source_benchmark_file": "",
        "source_constraint_kind": TEMPLATE_DOCKING_KIND,
        "coupling_source_kind": TEMPLATE_DOCKING_SOURCE_KIND,
        "coordinate_truth_used_to_build_constraints": False,
        "native_truth_used_before_coupling_selection": False,
        "external_evolutionary_couplings_used": False,
        "raw_sequence_exposed": False,
        "constraints": constraints,
        "structure_model_used_before_coupling_selection": False,
    }


def _map_template_worker(
    args: tuple[TemplateStructure, str, float, int, float | None, int],
) -> tuple[str, tuple[ContactPair, ...], TemplateInterfaceResult, str | None]:
    template, target_sequence, contact_cutoff_angstrom, minimum_sequence_separation, chemical_threshold, domain_scan_step = args
    mapped_pairs, info = map_template_interface_to_target(
        template,
        target_sequence,
        contact_cutoff_angstrom=contact_cutoff_angstrom,
        minimum_sequence_separation=minimum_sequence_separation,
        chemical_threshold=chemical_threshold,
        domain_scan_step=domain_scan_step,
    )
    return template.template_id, mapped_pairs, info, template.chain_id


def run_template_docking_v0(
    *,
    target_row: RealCoordinateVisualRow,
    templates: Sequence[TemplateStructure],
    long_range_threshold: int = 24,
    minimum_sequence_separation: int = DEFAULT_MIN_SEQUENCE_SEPARATION,
    contact_cutoff_angstrom: float = DEFAULT_CONTACT_CUTOFF_ANGSTROM,
    chemical_threshold: float | None = CHEMICAL_THRESHOLD_AUTO,
    min_template_support: int = 1,
    domain_scan_step: int = 1,
    source_mode: str = "benchmark_row_templates",
    jobs: int = 1,
    use_internal_gap: bool = True,
) -> tuple[TemplateDockingReport, dict[str, object], dict[str, object]]:
    target_sequence = target_row.sequence

    template_rows: list[dict[str, object]] = []
    template_votes: dict[str, tuple[ContactPair, ...]] = {}
    template_weights: dict[str, float] = {}
    vote_to_templates: dict[ContactPair, list[str]] = {}

    process_count = max(1, jobs if jobs > 0 else (cpu_count() or 1))
    templates_to_use = list(templates)
    worker_args = [
        (
            template,
            target_sequence,
            contact_cutoff_angstrom,
            minimum_sequence_separation,
            chemical_threshold,
            domain_scan_step,
        )
        for template in templates_to_use
    ]
    if process_count == 1 or len(templates_to_use) <= 1:
        mapped_results = [_map_template_worker(args) for args in worker_args]
    else:
        with Pool(processes=min(process_count, len(templates_to_use))) as pool:
            mapped_results = pool.map(_map_template_worker, worker_args)

    for template_id, mapped_pairs, info, template_chain in mapped_results:
        template_rows.append(
            {
                "template_id": info.template_id,
                "template_source_accession": info.source_accession,
                "template_source_id": info.source_id,
                "template_kind": info.source_kind,
                "source_path": info.source_path,
                "domain_split": info.domain_split,
                "interface_contact_count": info.interface_contact_count,
                "mapped_contact_count": len(info.mapped_pairs),
                "chem_kept_contact_count": info.chem_kept_contact_count,
                "mapping_mode": info.mapping_mode,
                "alignment_score": info.alignment_score,
                "alignment_identity": info.alignment_identity,
                "template_chain": template_chain,
            }
        )
        if info.chem_kept_contact_count <= 0:
            continue
        template_votes[info.template_id] = tuple(info.voted_pairs)
        template_weights[info.template_id] = 1.0
        for pair in info.voted_pairs:
            vote_to_templates.setdefault(pair, []).append(info.template_id)

    ranked = rank_template_pairs(
        template_votes,
        template_weights,
        long_range_threshold=long_range_threshold,
    )

    winner_id = ""
    winner_pairs: tuple[ContactPair, ...] = tuple()
    winner_source = ""
    winner_score = 0
    if template_votes:
        winner_scores: dict[str, float] = {}
        for template in template_rows:
            template_id = template["template_id"]
            pairs = template_votes.get(str(template_id), tuple())
            if not pairs:
                continue
            pair_score = sum(1 for _ in pairs)
            candidate = (
                pair_score * 1000.0
                + float(template.get("alignment_identity", 0.0)) * 500.0
                + float(template.get("alignment_score", 0.0)) * 0.01
            )
            winner_scores[template_id] = candidate
        if winner_scores:
            winner_id, winner_score = max(winner_scores.items(), key=lambda item: item[1])
            winner_info = next(
                row for row in template_rows if row["template_id"] == winner_id
            )
            winner_source = str(winner_info["template_source_accession"])
            winner_pairs = template_votes.get(winner_id, tuple())

    pair_scores = tuple((pair, score, support) for pair, score, support in ranked)

    min_support = max(int(min_template_support), 1)
    voted_pairs: list[ContactPair] = []
    voted_meta: list[dict[str, object]] = []
    for pair, score, support in pair_scores:
        if support >= min_support:
            voted_pairs.append(pair)
            voted_meta.append(
                {
                    "i": pair[0],
                    "j": pair[1],
                    "support_count": support,
                    "pair_score": score,
                    "template_count": support,
                    "template_ids": tuple(vote_to_templates.get(pair, tuple())),
                }
            )

    # If the strict support filter filtered too much, back off to winner-template only.
    if not voted_pairs and winner_pairs:
        voted_pairs = list(winner_pairs)
        for pair in winner_pairs:
            voted_meta.append(
                {
                    "i": pair[0],
                    "j": pair[1],
                    "support_count": 1,
                    "pair_score": 0.0,
                    "template_count": 1,
                    "template_ids": (winner_id,),
                }
            )

    voted_pairs = sorted(set(voted_pairs))
    winner_pairs = sorted(set(winner_pairs))

    predicted_pairs = tuple(voted_pairs)

    metric = evaluate_contact_prediction(
        native_pairs=target_row.native_contact_pairs(),
        predicted_pairs=predicted_pairs,
        long_range_threshold=long_range_threshold,
    )

    pair_support_count: dict[ContactPair, int] = {
        pair: len(vote_to_templates.get(pair, tuple())) for pair in predicted_pairs
    }
    pair_template_map: dict[ContactPair, tuple[str, ...]] = {
        pair: tuple(vote_to_templates.get(pair, tuple())) for pair in predicted_pairs
    }
    constraint_count = len(template_rows)
    constraint_payload = _build_constraint_payload(
        target_row,
        pairs=predicted_pairs,
        pair_support=pair_support_count,
        pair_templates=pair_template_map,
        template_count=constraint_count,
    )

    voted_pairs_payload = [
        {
            "i": pair[0],
            "j": pair[1],
            "support_count": pair_support_count[pair],
            "vote_count": pair_support_count[pair],
            "template_ids": list(pair_template_map[pair]),
        }
        for pair in sorted(predicted_pairs)
    ]

    winner_payload = [
        {"i": pair[0], "j": pair[1], "rank": index + 1}
        for index, pair in enumerate(winner_pairs)
    ]

    metric_dict = metric.to_dict()
    report = TemplateDockingReport(
        kind=TEMPLATE_DOCKING_KIND,
        source_mode=source_mode,
        target_source_accession=target_row.source_accession,
        target_row_id=target_row.row_id,
        target_sequence_length=target_row.sequence_length,
        template_count=len(template_rows),
        templates_used=sum(1 for rows in template_rows if rows["chem_kept_contact_count"] > 0),
        templates_with_mapped_pairs=sum(1 for rows in template_rows if rows["mapped_contact_count"] > 0),
        winner_template_id=winner_id,
        winner_template_source=winner_source,
        winner_vote_count=winner_score,
        min_template_support=min_support,
        predicted_pair_count=len(predicted_pairs),
        voted_pair_count=len(voted_pairs),
        winner_pair_count=len(winner_pairs),
        long_range_recall=float(metric.long_range_contact_recall),
        native_contact_precision=float(metric.native_contact_precision),
        native_contact_recall=float(metric.native_contact_recall),
        contact_map_f1=float(metric.contact_map_f1),
        false_contact_rate=float(metric.false_contact_rate),
        target_accession=target_row.source_accession,
        certificate_path=None,
        template_rows=tuple(template_rows),
        voted_pairs=tuple(voted_pairs_payload),
        winner_pairs=tuple(winner_payload),
        predicted_pairs=tuple(
            {
                "i": pair[0],
                "j": pair[1],
                "templates": list(pair_template_map[pair]),
            }
            for pair in sorted(predicted_pairs)
        ),
        metric_after_native_audit=metric_dict,
        long_range_threshold=long_range_threshold,
    )

    return report, constraint_payload, constraint_payload
