from __future__ import annotations

"""True local MSA-conditioned co-evolution expansion around DCA anchors.

This module is deliberately different from the earlier proxy experiment:

* If a raw MSA is supplied, local coupling is computed from the filtered
  sub-MSA around each anchor.
* If no raw MSA is supplied, the runner must report ACTION_REQUIRED instead of
  silently falling back to a proxy.
* Native contacts/coordinates are never used before the selected contact map is
  frozen. Native is attached only for audit.

The local score is a normalized mutual-information signal computed only from
sequences where the anchor residues match the query/reference residues. This is
not plmDCA, but it is the real testable local-coevolution channel that can run
without external Python packages.
"""

from dataclasses import asdict, dataclass
from math import log, sqrt
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from pharmacotopology.folding_coarse_grain_md_geometry import CoarseGrainMDContactDecision, _matched_controls_for_report, _pair_hash, _rounded
from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_native_contact_eval import ContactMetricPacket, ContactPair, evaluate_contact_prediction, normalized_contact_pairs
from pharmacotopology.folding_real_coordinate_visual_benchmark import RealCoordinateVisualRow

TRUE_LOCAL_MSA_KIND = "true_local_msa_coevolution_expansion_v0"
TRUE_LOCAL_MSA_MODE = "raw_msa_anchor_conditioned_local_mi"
TRUE_LOCAL_MSA_DECISION_KIND = "true_local_msa_coevolution_contact_decision_v0"
TRUE_LOCAL_MSA_RULE = (
    "requires_raw_msa;top_safe_dca_anchors;filter_sub_msa_by_query_anchor_residues;"
    "compute_normalized_mi_in_anchor_local_window;degree_cap;native_audit_after_selection"
)
TRUE_LOCAL_MSA_CLAIM_RULE = "claim_requires_precision_recall_long_range_ge_0_70_and_control_margins"

LONG_RANGE_THRESHOLD = 24
LOCAL_LONG_RANGE_FLOOR = 30
GAP_CHARS = set("-.")
AA_ALPHABET = tuple("ACDEFGHIKLMNPQRSTVWYXBZUO")
AutoThreshold = Union[float, str]
AUTO_THRESHOLD_VALUES = frozenset({"auto", "internal-gap", "internal_gap", "internalgap"})


@dataclass(frozen=True)
class ParsedMSA:
    path: str
    format: str
    records: Tuple[Tuple[str, str], ...]
    aligned_length: int
    query_sequence: str
    query_to_alignment: Mapping[int, int]

    @property
    def sequence_count(self) -> int:
        return len(self.records)


@dataclass(frozen=True)
class TrueLocalMSAContact:
    row_id: str
    source_accession: str
    i: int
    j: int
    score: float
    channel: str
    parent_anchor_i: int
    parent_anchor_j: int
    offset_i: int
    offset_j: int
    filtered_sequence_count: int
    raw_mi: float
    normalized_mi: float
    anchor_i_residue: str
    anchor_j_residue: str
    query_i_residue: str
    query_j_residue: str
    true_msa_mi_used: bool = True
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class TrueLocalMSARowReport:
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    msa_path: str
    raw_msa_available_for_true_local_mi: bool
    msa_sequence_count: int
    msa_aligned_length: int
    top_safe_dca_anchor_count_requested: int
    safe_anchor_count: int
    local_window: int
    threshold: float
    min_filtered_sequences: int
    candidate_local_pair_count: int
    accepted_local_pair_count: int
    selected_contact_count: int
    selected_long_range_contact_count: int
    selected_contact_map_hash: str
    metric_after_native_audit: ContactMetricPacket
    matched_control_count: int
    best_control_f1_after_audit: float
    best_control_long_range_recall_after_audit: float
    f1_margin_vs_best_control: float
    long_range_recall_margin_vs_best_control: float
    row_true_local_msa_claim_allowed: bool
    row_universal_physical_law_claim_allowed: bool
    row_claim_rejection_reason: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    target_native_truth_attached_after_selection_for_evaluation: bool = True
    learned_geometry_prior_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = True

    def to_dict(self) -> dict:
        payload = asdict(self)
        payload["metric_after_native_audit"] = self.metric_after_native_audit.to_dict()
        return payload


@dataclass(frozen=True)
class TrueLocalMSAPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    msa_path: str
    raw_msa_available_for_true_local_mi: bool
    top_safe_dca_anchor_count_requested: int
    local_window: int
    threshold: float
    min_filtered_sequences: int
    mean_safe_anchor_count: float
    mean_accepted_local_pair_count: float
    mean_native_contact_precision_after_audit: float
    mean_native_contact_recall_after_audit: float
    mean_long_range_contact_recall_after_audit: float
    mean_contact_map_f1_after_audit: float
    mean_f1_margin_vs_best_control: float
    mean_long_range_recall_margin_vs_best_control: float
    true_local_msa_claim_allowed: bool
    universal_physical_law_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    rows: Tuple[TrueLocalMSARowReport, ...]
    contacts: Tuple[TrueLocalMSAContact, ...]
    decisions: Tuple[CoarseGrainMDContactDecision, ...]

    def to_dict(self) -> dict:
        return {
            "kind": self.kind,
            "source_mode": self.source_mode,
            "row_count": self.row_count,
            "decision_rule": self.decision_rule,
            "claim_rule": self.claim_rule,
            "msa_path": self.msa_path,
            "raw_msa_available_for_true_local_mi": self.raw_msa_available_for_true_local_mi,
            "top_safe_dca_anchor_count_requested": self.top_safe_dca_anchor_count_requested,
            "local_window": self.local_window,
            "threshold": self.threshold,
            "min_filtered_sequences": self.min_filtered_sequences,
            "mean_safe_anchor_count": self.mean_safe_anchor_count,
            "mean_accepted_local_pair_count": self.mean_accepted_local_pair_count,
            "mean_native_contact_precision_after_audit": self.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": self.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": self.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": self.mean_contact_map_f1_after_audit,
            "mean_f1_margin_vs_best_control": self.mean_f1_margin_vs_best_control,
            "mean_long_range_recall_margin_vs_best_control": self.mean_long_range_recall_margin_vs_best_control,
            "true_local_msa_claim_allowed": self.true_local_msa_claim_allowed,
            "universal_physical_law_claim_allowed": self.universal_physical_law_claim_allowed,
            "folding_problem_solved": self.folding_problem_solved,
            "claim_rejection_reason": self.claim_rejection_reason,
            "rows": [r.to_dict() for r in self.rows],
            "contacts": [c.to_dict() for c in self.contacts],
            "decisions": [d.to_dict() for d in self.decisions],
        }


def _bounded01(value: float) -> float:
    return round(max(0.0, min(1.0, float(value))), 6)


def _mean(values: Sequence[float]) -> float:
    return _rounded(mean(values)) if values else 0.0


def _ungap(sequence: str) -> str:
    return "".join(ch for ch in sequence.strip().upper() if ch not in GAP_CHARS and not ch.isspace())


def _read_fasta_like(path: Path) -> Tuple[Tuple[str, str], ...]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    chunks: List[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper()))
            name = line[1:].strip() or "sequence_%d" % (len(records) + 1)
            chunks = []
        elif not line.startswith("#"):
            chunks.append(line.replace(" ", ""))
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return tuple(records)


def _read_stockholm(path: Path) -> Tuple[Tuple[str, str], ...]:
    seqs: Dict[str, List[str]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.rstrip()
        if not line or line.startswith("#") or line == "//":
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        name, chunk = parts[0], parts[1]
        if name.startswith("#"):
            continue
        seqs.setdefault(name, []).append(chunk)
    return tuple((name, "".join(chunks).upper()) for name, chunks in seqs.items())


def parse_msa(path: Path, query_sequence: str) -> ParsedMSA:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError("raw MSA file not found: %s" % path)
    suffix = path.suffix.lower()
    if suffix in (".sto", ".stockholm"):
        records = _read_stockholm(path)
        fmt = "stockholm"
    else:
        records = _read_fasta_like(path)
        fmt = "fasta_like"
    if not records:
        raise ValueError("MSA file contains no sequence records: %s" % path)
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError("MSA records must be aligned to equal length; got lengths=%s" % sorted(lengths)[:10])
    query = query_sequence.upper().replace(" ", "")
    query_ungapped = _ungap(query)
    query_record = records[0][1]
    for _name, seq in records:
        if _ungap(seq).startswith(query_ungapped) or query_ungapped.startswith(_ungap(seq)):
            query_record = seq
            break
    mapping: Dict[int, int] = {}
    residue_index = 0
    for aln_index, residue in enumerate(query_record):
        if residue not in GAP_CHARS:
            residue_index += 1
            mapping[residue_index] = aln_index
    missing = [idx for idx in range(1, len(query_ungapped) + 1) if idx not in mapping]
    if missing:
        raise ValueError("MSA query mapping is incomplete; first missing residue=%s" % missing[0])
    return ParsedMSA(
        path=str(path),
        format=fmt,
        records=records,
        aligned_length=next(iter(lengths)),
        query_sequence=query_ungapped,
        query_to_alignment=mapping,
    )


def _safe_anchors(row: RealCoordinateVisualRow, constraints: Sequence[CouplingConstraint], top_n: int) -> Tuple[CouplingConstraint, ...]:
    safe = []
    for c in constraints:
        if (c.row_id != row.row_id and c.source_accession != row.source_accession) or c.i >= c.j:
            continue
        if c.coordinate_truth_used_to_build_constraint or c.native_truth_used_before_coupling_selection or c.structure_model_used:
            continue
        if c.sequence_separation < LONG_RANGE_THRESHOLD:
            continue
        safe.append(c)
    safe.sort(key=lambda c: (-(c.apc_corrected_score or c.confidence), c.rank if c.rank else 999999, c.i, c.j))
    return tuple(safe[: max(1, top_n)])


def _column_values(msa: ParsedMSA, residue_index: int, allowed_record_indexes: Sequence[int]) -> List[str]:
    aln_index = msa.query_to_alignment[residue_index]
    values = []
    for record_index in allowed_record_indexes:
        residue = msa.records[record_index][1][aln_index].upper()
        if residue in GAP_CHARS:
            continue
        values.append(residue)
    return values


def _normalized_mi(xs: Sequence[str], ys: Sequence[str], pseudocount: float = 0.5) -> Tuple[float, float]:
    pairs = [(a, b) for a, b in zip(xs, ys) if a not in GAP_CHARS and b not in GAP_CHARS]
    if len(pairs) < 4:
        return 0.0, 0.0
    # Use the observed local alphabet, not the full amino-acid alphabet.
    # Full-alphabet pseudocounts over-dampen small anchor-conditioned
    # sub-MSAs and erase the very signal this experiment is meant to test.
    alphabet_x = sorted(set(a for a, _b in pairs))
    alphabet_y = sorted(set(b for _a, b in pairs))
    total = float(len(pairs)) + pseudocount * len(alphabet_x) * len(alphabet_y)
    px: Dict[str, float] = dict((a, pseudocount * len(alphabet_y)) for a in alphabet_x)
    py: Dict[str, float] = dict((b, pseudocount * len(alphabet_x)) for b in alphabet_y)
    pxy: Dict[Tuple[str, str], float] = {}
    for a in alphabet_x:
        for b in alphabet_y:
            pxy[(a, b)] = pseudocount
    for a, b in pairs:
        px[a] = px.get(a, 0.0) + 1.0
        py[b] = py.get(b, 0.0) + 1.0
        pxy[(a, b)] = pxy.get((a, b), 0.0) + 1.0
    mi = 0.0
    hx = 0.0
    hy = 0.0
    for a in alphabet_x:
        pa = px[a] / total
        if pa > 0:
            hx -= pa * log(pa)
    for b in alphabet_y:
        pb = py[b] / total
        if pb > 0:
            hy -= pb * log(pb)
    for a in alphabet_x:
        for b in alphabet_y:
            pab = pxy[(a, b)] / total
            if pab <= 0:
                continue
            pa = px[a] / total
            pb = py[b] / total
            mi += pab * log(pab / (pa * pb))
    nmi = mi / max(1e-9, sqrt(hx * hy))
    return _bounded01(nmi), round(mi, 6)


def _resolve_internal_gap_threshold(scores: Sequence[float]) -> float:
    if not scores:
        return 0.0
    if len(scores) == 1:
        return scores[0]
    ordered = sorted(scores, reverse=True)
    gaps = [ordered[index] - ordered[index + 1] for index in range(len(ordered) - 1)]
    if not gaps:
        return ordered[-1]
    max_gap_index = max(range(len(gaps)), key=lambda index: gaps[index])
    return ordered[max_gap_index + 1]


def _resolve_threshold(
    *,
    threshold: AutoThreshold,
    candidate_scores: Sequence[float],
) -> float:
    if isinstance(threshold, str) and threshold.lower() in AUTO_THRESHOLD_VALUES:
        return _bounded01(_resolve_internal_gap_threshold(candidate_scores))
    return _bounded01(float(threshold))


def _records_matching_anchor(msa: ParsedMSA, row: RealCoordinateVisualRow, anchor: CouplingConstraint) -> Tuple[List[int], str, str]:
    ai = msa.query_to_alignment[anchor.i]
    aj = msa.query_to_alignment[anchor.j]
    ref_i = row.sequence[anchor.i - 1].upper()
    ref_j = row.sequence[anchor.j - 1].upper()
    indexes = []
    for idx, (_name, seq) in enumerate(msa.records):
        if ai >= len(seq) or aj >= len(seq):
            continue
        if seq[ai].upper() == ref_i and seq[aj].upper() == ref_j:
            indexes.append(idx)
    return indexes, ref_i, ref_j


def build_true_local_msa_contacts(
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    msa: ParsedMSA,
    top_anchor_count: int = 50,
    window: int = 5,
    threshold: AutoThreshold = 0.5,
    degree_cap: int = 5,
    min_filtered_sequences: int = 16,
) -> Tuple[
    Tuple[TrueLocalMSAContact, ...],
    Tuple[ContactPair, ...],
    Dict[ContactPair, float],
    float,
    int,
]:
    anchors = _safe_anchors(row, constraints, top_n=top_anchor_count)
    accepted: Dict[ContactPair, TrueLocalMSAContact] = {}
    scores: Dict[ContactPair, float] = {}
    local_candidates: list[tuple[float, TrueLocalMSAContact]] = []

    def add(record: TrueLocalMSAContact) -> None:
        pair = record.pair()
        old = accepted.get(pair)
        if old is None or record.score > old.score:
            accepted[pair] = record
            scores[pair] = record.score

    for anchor in anchors:
        filtered, anchor_i_res, anchor_j_res = _records_matching_anchor(msa, row, anchor)
        if len(filtered) < min_filtered_sequences:
            continue
        direct_score = _bounded01(max(float(anchor.confidence), float(anchor.apc_corrected_score or 0.0)))
        add(TrueLocalMSAContact(
            row_id=row.row_id,
            source_accession=row.source_accession,
            i=anchor.i,
            j=anchor.j,
            score=direct_score,
            channel="direct_safe_dca_anchor_with_msa_filter",
            parent_anchor_i=anchor.i,
            parent_anchor_j=anchor.j,
            offset_i=0,
            offset_j=0,
            filtered_sequence_count=len(filtered),
            raw_mi=direct_score,
            normalized_mi=direct_score,
            anchor_i_residue=anchor_i_res,
            anchor_j_residue=anchor_j_res,
            query_i_residue=anchor_i_res,
            query_j_residue=anchor_j_res,
        ))
        for di in range(-window, window + 1):
            for dj in range(-window, window + 1):
                if di == 0 and dj == 0:
                    continue
                i2 = anchor.i + di
                j2 = anchor.j + dj
                if i2 < 1 or j2 > row.sequence_length or i2 >= j2 or j2 - i2 < LOCAL_LONG_RANGE_FLOOR:
                    continue
                xs = _column_values(msa, i2, filtered)
                ys = _column_values(msa, j2, filtered)
                nmi, raw_mi = _normalized_mi(xs, ys)
                # Anchor-conditioned local MI: the original anchor must define
                # the sub-MSA, but selection is decided by the candidate pair's
                # MI in that filtered sub-MSA.
                local_candidates.append((nmi, TrueLocalMSAContact(
                    row_id=row.row_id,
                    source_accession=row.source_accession,
                    i=i2,
                    j=j2,
                    score=nmi,
                    channel="true_local_msa_mi_window",
                    parent_anchor_i=anchor.i,
                    parent_anchor_j=anchor.j,
                    offset_i=di,
                    offset_j=dj,
                    filtered_sequence_count=len(filtered),
                    raw_mi=raw_mi,
                    normalized_mi=nmi,
                    anchor_i_residue=anchor_i_res,
                    anchor_j_residue=anchor_j_res,
                    query_i_residue=row.sequence[i2 - 1].upper(),
                    query_j_residue=row.sequence[j2 - 1].upper(),
                )))
    resolved_threshold = _resolve_threshold(
        threshold=threshold,
        candidate_scores=tuple(score for score, _ in local_candidates),
    )
    for score, record in local_candidates:
        if score >= resolved_threshold:
            add(record)
    ranked = sorted(accepted.values(), key=lambda r: (-r.score, r.channel != "direct_safe_dca_anchor_with_msa_filter", r.i, r.j))
    degrees = [0] * row.sequence_length
    selected: List[ContactPair] = []
    selected_scores: Dict[ContactPair, float] = {}
    for record in ranked:
        i, j = record.pair()
        if degrees[i - 1] >= degree_cap or degrees[j - 1] >= degree_cap:
            continue
        selected.append((i, j))
        selected_scores[(i, j)] = record.score
        degrees[i - 1] += 1
        degrees[j - 1] += 1
    return tuple(ranked), normalized_contact_pairs(selected), selected_scores, resolved_threshold, len(local_candidates)


def _claim(metric: ContactMetricPacket, f1_margin: float, lr_margin: float) -> Tuple[bool, bool, str]:
    claim = (
        metric.native_contact_precision >= 0.70
        and metric.native_contact_recall >= 0.70
        and metric.long_range_contact_recall >= 0.70
        and f1_margin > 0.03
        and lr_margin > 0.03
    )
    reason = "true_local_msa_coevolution_survived_gate" if claim else "true_local_msa_claim_rejected_precision_recall_or_long_range_below_0_70"
    # Uses MSA/evolutionary data, therefore not a sequence-only physical law.
    return claim, False, reason


def _decisions(row: RealCoordinateVisualRow, selected: Sequence[ContactPair], scores: Mapping[ContactPair, float]) -> Tuple[CoarseGrainMDContactDecision, ...]:
    out = []
    for i, j in normalized_contact_pairs(selected):
        out.append(CoarseGrainMDContactDecision(
            kind=TRUE_LOCAL_MSA_DECISION_KIND,
            row_id=row.row_id,
            source_accession=row.source_accession,
            source_mode=TRUE_LOCAL_MSA_MODE,
            i=i,
            j=j,
            sequence_separation=j - i,
            final_distance_angstrom=0.0,
            geometry_contact_score=_bounded01(scores.get((i, j), 0.0)),
            selected=True,
            selected_from_final_structure=False,
            msa_dca_used_before_selection=True,
        ))
    return tuple(out)


def run_true_local_msa_row(
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    msa: ParsedMSA,
    top_anchor_count: int = 50,
    window: int = 5,
    threshold: AutoThreshold = 0.5,
    degree_cap: int = 5,
    min_filtered_sequences: int = 16,
) -> Tuple[TrueLocalMSARowReport, Tuple[TrueLocalMSAContact, ...], Tuple[CoarseGrainMDContactDecision, ...]]:
    (
        contacts,
        selected,
        scores,
        resolved_threshold,
        candidate_local_pair_count,
    ) = build_true_local_msa_contacts(
        row,
        constraints,
        msa,
        top_anchor_count=top_anchor_count,
        window=window,
        threshold=threshold,
        degree_cap=degree_cap,
        min_filtered_sequences=min_filtered_sequences,
    )
    metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=selected)
    best_control_f1, best_control_lr, control_count = _matched_controls_for_report(row=row, selected_pairs=selected)
    f1_margin = _rounded(metric.contact_map_f1 - best_control_f1)
    lr_margin = _rounded(metric.long_range_contact_recall - best_control_lr)
    claim, universal, reason = _claim(metric, f1_margin, lr_margin)
    anchors = _safe_anchors(row, constraints, top_n=top_anchor_count)
    report = TrueLocalMSARowReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        msa_path=msa.path,
        raw_msa_available_for_true_local_mi=True,
        msa_sequence_count=msa.sequence_count,
        msa_aligned_length=msa.aligned_length,
        top_safe_dca_anchor_count_requested=top_anchor_count,
        safe_anchor_count=len(anchors),
        local_window=window,
        threshold=_rounded(resolved_threshold),
        min_filtered_sequences=min_filtered_sequences,
        candidate_local_pair_count=candidate_local_pair_count,
        accepted_local_pair_count=sum(1 for c in contacts if c.channel == "true_local_msa_mi_window"),
        selected_contact_count=len(selected),
        selected_long_range_contact_count=sum(1 for i, j in selected if j - i >= LONG_RANGE_THRESHOLD),
        selected_contact_map_hash=_pair_hash(selected),
        metric_after_native_audit=metric,
        matched_control_count=control_count,
        best_control_f1_after_audit=best_control_f1,
        best_control_long_range_recall_after_audit=best_control_lr,
        f1_margin_vs_best_control=f1_margin,
        long_range_recall_margin_vs_best_control=lr_margin,
        row_true_local_msa_claim_allowed=claim,
        row_universal_physical_law_claim_allowed=universal,
        row_claim_rejection_reason=reason,
    )
    return report, contacts, _decisions(row, selected, scores)


def run_true_local_msa_packet(
    rows: Sequence[RealCoordinateVisualRow],
    constraints: Sequence[CouplingConstraint],
    msa_path: Path,
    evaluation_source_accessions: Optional[Sequence[str]] = None,
    top_anchor_count: int = 50,
    window: int = 5,
    threshold: AutoThreshold = 0.5,
    degree_cap: int = 5,
    min_filtered_sequences: int = 16,
) -> TrueLocalMSAPacket:
    evaluation_set = set(evaluation_source_accessions or ())
    eval_rows = tuple(row for row in rows if not evaluation_set or row.source_accession in evaluation_set)
    if len(eval_rows) != 1:
        raise ValueError("true local MSA runner currently expects one source-accession per raw MSA; got %d rows" % len(eval_rows))
    row = eval_rows[0]
    msa = parse_msa(Path(msa_path), row.sequence)
    reports = []
    all_contacts = []
    all_decisions = []
    report, contacts, decisions = run_true_local_msa_row(
        row,
        constraints,
        msa,
        top_anchor_count=top_anchor_count,
        window=window,
        threshold=threshold,
        degree_cap=degree_cap,
        min_filtered_sequences=min_filtered_sequences,
    )
    reports.append(report)
    all_contacts.extend(contacts[:2000])
    all_decisions.extend(decisions[:1000])
    precisions = [r.metric_after_native_audit.native_contact_precision for r in reports]
    recalls = [r.metric_after_native_audit.native_contact_recall for r in reports]
    long_ranges = [r.metric_after_native_audit.long_range_contact_recall for r in reports]
    f1s = [r.metric_after_native_audit.contact_map_f1 for r in reports]
    f1_margins = [r.f1_margin_vs_best_control for r in reports]
    lr_margins = [r.long_range_recall_margin_vs_best_control for r in reports]
    claim = bool(reports) and all(r.row_true_local_msa_claim_allowed for r in reports)
    reason = "true_local_msa_coevolution_survived_gate" if claim else "true_local_msa_claim_rejected_for_rows:" + ",".join(r.source_accession for r in reports if not r.row_true_local_msa_claim_allowed)
    return TrueLocalMSAPacket(
        kind=TRUE_LOCAL_MSA_KIND,
        source_mode=TRUE_LOCAL_MSA_MODE,
        row_count=len(reports),
        decision_rule=TRUE_LOCAL_MSA_RULE,
        claim_rule=TRUE_LOCAL_MSA_CLAIM_RULE,
        msa_path=str(msa_path),
        raw_msa_available_for_true_local_mi=True,
        top_safe_dca_anchor_count_requested=top_anchor_count,
        local_window=window,
        threshold=_rounded(reports[0].threshold if reports else 0.0),
        min_filtered_sequences=min_filtered_sequences,
        mean_safe_anchor_count=_mean([r.safe_anchor_count for r in reports]),
        mean_accepted_local_pair_count=_mean([r.accepted_local_pair_count for r in reports]),
        mean_native_contact_precision_after_audit=_mean(precisions),
        mean_native_contact_recall_after_audit=_mean(recalls),
        mean_long_range_contact_recall_after_audit=_mean(long_ranges),
        mean_contact_map_f1_after_audit=_mean(f1s),
        mean_f1_margin_vs_best_control=_mean(f1_margins),
        mean_long_range_recall_margin_vs_best_control=_mean(lr_margins),
        true_local_msa_claim_allowed=claim,
        universal_physical_law_claim_allowed=False,
        folding_problem_solved=claim,
        claim_rejection_reason=reason,
        rows=tuple(reports),
        contacts=tuple(all_contacts),
        decisions=tuple(all_decisions),
    )
