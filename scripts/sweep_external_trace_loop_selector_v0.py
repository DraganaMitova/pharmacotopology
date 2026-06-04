from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.artifact_io import write_csv_rows  # noqa: E402
from pharmacotopology.folding_coupling_negative_controls import (  # noqa: E402
    generate_external_coupling_negative_controls,
)
from pharmacotopology.folding_coupling_nucleus_selector import (  # noqa: E402
    SELECTED_EVENTS_PER_ROW,
    TRACE_LOOP_MARGIN_GATE_BLOCKED_FUTURE_MAX,
    TRACE_LOOP_MARGIN_GATE_MIN,
    CouplingNucleusContext,
    build_coupling_nucleus_context,
    coupling_nucleus_score,
    select_coupling_events,
    selector_metrics,
)
from pharmacotopology.folding_evolutionary_constraints import (  # noqa: E402
    CouplingConstraint,
    CouplingDataset,
    compatible_future_event,
    load_coupling_dataset,
)
from pharmacotopology.folding_external_coupling_importer import (  # noqa: E402
    import_external_coupling_dataset,
)
from pharmacotopology.folding_external_coupling_sources import (  # noqa: E402
    SERIOUS_EXTERNAL_COUPLING_POLICY,
)
from pharmacotopology.folding_external_coupling_trace_loop import (  # noqa: E402
    MATCHED_CONTROL_NAMES,
)
from pharmacotopology.folding_nucleus_closure_search import (  # noqa: E402
    NucleusClosureEvent,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_EXTERNAL_COUPLING_FILE = Path(
    "data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
)
DEFAULT_ORACLE_COUPLING_FILE = Path(
    "data/folding_real_coordinate_visual_8_couplings.locked.json"
)
DEFAULT_OUTPUT = Path(
    "first_contact_clean_pharmacotopology_layer_run/"
    "external_trace_loop_selector_sweep_v0.csv"
)


@dataclass(frozen=True)
class SelectorSweepConfig:
    selector_config: str
    score_mode: str
    top_fraction: float
    min_new_top_pairs: int
    min_candidate_top_pairs: int
    min_mean_new_confidence: float
    margin_gated: bool
    max_blocked_future_pressure: float | None
    confidence_mass_weight: float
    rank_mass_weight: float
    top_pair_weight: float
    mean_confidence_weight: float
    nucleus_weight: float


@dataclass(frozen=True)
class RowPairFeatures:
    confidence_by_pair: Mapping[tuple[int, int], float]
    rank_signal_by_pair: Mapping[tuple[int, int], float]
    top_pairs: frozenset[tuple[int, int]]


def _rounded(value: float) -> float:
    return round(value, 6)


def _events_by_row(
    events: Sequence[NucleusClosureEvent],
) -> dict[str, tuple[NucleusClosureEvent, ...]]:
    grouped: dict[str, list[NucleusClosureEvent]] = {}
    for event in events:
        grouped.setdefault(event.row_id, []).append(event)
    return {row_id: tuple(row_events) for row_id, row_events in grouped.items()}


def _row_pair_features(
    constraints: Sequence[CouplingConstraint],
    *,
    top_fraction: float,
) -> RowPairFeatures:
    ranked = tuple(
        sorted(
            constraints,
            key=lambda constraint: (
                -constraint.confidence,
                constraint.i,
                constraint.j,
                constraint.constraint_id,
            ),
        )
    )
    confidence_by_pair = {constraint.pair(): constraint.confidence for constraint in ranked}
    denominator = max(1, len(ranked) - 1)
    rank_signal_by_pair = {
        constraint.pair(): 1.0 - ((rank - 1) / denominator)
        for rank, constraint in enumerate(ranked, start=1)
    }
    top_count = max(1, int(len(ranked) * top_fraction))
    top_pairs = frozenset(constraint.pair() for constraint in ranked[:top_count])
    return RowPairFeatures(
        confidence_by_pair=confidence_by_pair,
        rank_signal_by_pair=rank_signal_by_pair,
        top_pairs=top_pairs,
    )


def _trace_score(
    *,
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
    config: SelectorSweepConfig,
    newly_covered: set[tuple[int, int]],
    features: RowPairFeatures,
    top_new_count: int,
) -> float:
    assessment = context.assessment_by_event_id[event.event_id]
    covered_confidence = sum(features.confidence_by_pair[pair] for pair in newly_covered)
    rank_mass = sum(features.rank_signal_by_pair[pair] for pair in newly_covered)
    mean_confidence = covered_confidence / max(1, len(newly_covered))
    top_pair_signal = top_new_count / max(1, config.min_new_top_pairs or 2)
    return _rounded(
        config.confidence_mass_weight * min(1.0, covered_confidence / 3.0)
        + config.rank_mass_weight * min(1.0, rank_mass / 2.0)
        + config.top_pair_weight * min(1.0, top_pair_signal)
        + config.mean_confidence_weight * mean_confidence
        + config.nucleus_weight * coupling_nucleus_score(event, context)
        + 0.16 * assessment.future_preservation_score
        - 0.10 * assessment.blocked_future_pressure
    )


def select_sweep_trace_loop_events(
    context: CouplingNucleusContext,
    config: SelectorSweepConfig,
) -> tuple[NucleusClosureEvent, ...]:
    constraints_by_row = context.coupling_dataset.constraints_by_row_id()
    competitive_by_row = _events_by_row(context.competitive_events)
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_constraints = tuple(constraints_by_row.get(row.row_id, ()))
        if not row_constraints:
            continue
        features = _row_pair_features(
            row_constraints,
            top_fraction=config.top_fraction,
        )
        uncovered = set(features.confidence_by_pair)
        row_candidates = tuple(competitive_by_row.get(row.row_id, ()))
        row_selected: list[NucleusClosureEvent] = []
        selected_ids: set[str] = set()
        while uncovered and len(row_selected) < SELECTED_EVENTS_PER_ROW:
            scored_candidates: list[tuple[float, NucleusClosureEvent]] = []
            for event in row_candidates:
                if event.event_id in selected_ids:
                    continue
                if any(
                    not compatible_future_event(selected_event, event)
                    for selected_event in row_selected
                ):
                    continue
                event_pairs = set(event.candidate_region_pairs())
                newly_covered = event_pairs & uncovered
                if not newly_covered:
                    continue
                assessment = context.assessment_by_event_id[event.event_id]
                if (
                    config.max_blocked_future_pressure is not None
                    and assessment.blocked_future_pressure
                    > config.max_blocked_future_pressure
                ):
                    continue
                if (
                    config.margin_gated
                    and context.coupling_decoy_margin_by_event_id[event.event_id]
                    < TRACE_LOOP_MARGIN_GATE_MIN
                ):
                    continue
                top_new_count = len(newly_covered & features.top_pairs)
                if top_new_count < config.min_new_top_pairs:
                    continue
                top_candidate_count = len(event_pairs & features.top_pairs)
                if top_candidate_count < config.min_candidate_top_pairs:
                    continue
                mean_confidence = sum(
                    features.confidence_by_pair[pair] for pair in newly_covered
                ) / max(1, len(newly_covered))
                if mean_confidence < config.min_mean_new_confidence:
                    continue
                scored_candidates.append(
                    (
                        _trace_score(
                            event=event,
                            context=context,
                            config=config,
                            newly_covered=newly_covered,
                            features=features,
                            top_new_count=top_new_count,
                        ),
                        event,
                    )
                )
            if not scored_candidates:
                break
            _, chosen = max(
                scored_candidates,
                key=lambda item: (
                    item[0],
                    context.assessment_by_event_id[item[1].event_id].direct_support_score,
                    -item[1].segment_a_start,
                    -item[1].segment_b_start,
                    item[1].event_id,
                ),
            )
            row_selected.append(chosen)
            selected_ids.add(chosen.event_id)
            uncovered -= set(chosen.candidate_region_pairs())
        selected.extend(row_selected)
    return tuple(selected)


def _score_modes() -> Mapping[str, tuple[float, float, float, float, float]]:
    return {
        "current_like": (0.48, 0.0, 0.0, 0.0, 0.28),
        "rank_concentrated": (0.26, 0.24, 0.16, 0.08, 0.20),
        "mean_confidence": (0.20, 0.16, 0.12, 0.22, 0.22),
        "rank_heavy": (0.18, 0.34, 0.18, 0.08, 0.18),
    }


def build_sweep_configs() -> tuple[SelectorSweepConfig, ...]:
    configs: list[SelectorSweepConfig] = []
    for score_mode, weights in _score_modes().items():
        (
            confidence_mass_weight,
            rank_mass_weight,
            top_pair_weight,
            mean_confidence_weight,
            nucleus_weight,
        ) = weights
        for margin_gated in (False, True):
            for top_fraction in (0.10, 0.20, 0.30, 0.40):
                for min_new_top_pairs in (0, 1, 2):
                    for min_candidate_top_pairs in (0, 1):
                        for min_mean_new_confidence in (0.0, 0.20, 0.35):
                            name = (
                                f"{score_mode}"
                                f"_top{top_fraction:.2f}"
                                f"_new{min_new_top_pairs}"
                                f"_cand{min_candidate_top_pairs}"
                                f"_mean{min_mean_new_confidence:.2f}"
                                f"_{'margin' if margin_gated else 'open'}"
                            )
                            configs.append(
                                SelectorSweepConfig(
                                    selector_config=name,
                                    score_mode=score_mode,
                                    top_fraction=top_fraction,
                                    min_new_top_pairs=min_new_top_pairs,
                                    min_candidate_top_pairs=min_candidate_top_pairs,
                                    min_mean_new_confidence=min_mean_new_confidence,
                                    margin_gated=margin_gated,
                                    max_blocked_future_pressure=(
                                        TRACE_LOOP_MARGIN_GATE_BLOCKED_FUTURE_MAX
                                        if margin_gated
                                        else None
                                    ),
                                    confidence_mass_weight=confidence_mass_weight,
                                    rank_mass_weight=rank_mass_weight,
                                    top_pair_weight=top_pair_weight,
                                    mean_confidence_weight=mean_confidence_weight,
                                    nucleus_weight=nucleus_weight,
                                )
                            )
    return tuple(configs)


def _context_by_name(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    external_dataset: CouplingDataset,
) -> dict[str, CouplingNucleusContext]:
    controls = generate_external_coupling_negative_controls(
        rows=rows,
        dataset=external_dataset,
    )
    datasets = {"external_real": external_dataset}
    datasets.update({name: control.dataset for name, control in controls.items()})
    return {
        name: build_coupling_nucleus_context(rows=rows, coupling_dataset=dataset)
        for name, dataset in datasets.items()
    }


def _control_field(prefix: str, control_name: str) -> str:
    return f"{prefix}_{control_name.replace('external_', '')}"


def run_selector_sweep(
    *,
    benchmark_file: Path,
    external_coupling_file: Path,
    oracle_coupling_file: Path,
    output: Path,
    top_results_output: Path | None,
    max_configs: int | None,
    progress_every: int,
) -> tuple[Path, Path | None]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    import_result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=external_coupling_file,
        policy=SERIOUS_EXTERNAL_COUPLING_POLICY,
    )
    contexts = _context_by_name(rows=rows, external_dataset=import_result.dataset)
    external_context = contexts["external_real"]
    physical_selected = select_coupling_events(
        external_context,
        selector_name="physical_rerank",
    )
    physical_metric = selector_metrics(
        external_context,
        selector_name="physical_rerank",
        selected_events=physical_selected,
    )
    oracle_context = build_coupling_nucleus_context(
        rows=rows,
        coupling_dataset=load_coupling_dataset(oracle_coupling_file),
    )
    oracle_selected = select_coupling_events(
        oracle_context,
        selector_name="coupling_trace_loop",
    )
    oracle_metric = selector_metrics(
        oracle_context,
        selector_name="oracle_coordinate_positive_control",
        selected_events=oracle_selected,
    )
    oracle_recall_floor = _rounded(0.50 * oracle_metric.long_range_contact_recall)

    configs = build_sweep_configs()
    if max_configs is not None:
        configs = configs[:max_configs]

    rows_out: list[dict[str, object]] = []
    for index, config in enumerate(configs, start=1):
        if progress_every and (
            index == 1 or index == len(configs) or index % progress_every == 0
        ):
            print(
                f"sweeping selector config {index}/{len(configs)}: "
                f"{config.selector_config}",
                file=sys.stderr,
            )
        external_selected = select_sweep_trace_loop_events(external_context, config)
        external_metric = selector_metrics(
            external_context,
            selector_name=config.selector_config,
            selected_events=external_selected,
        )
        control_metrics = {}
        for control_name in MATCHED_CONTROL_NAMES:
            control_context = contexts[control_name]
            control_selected = select_sweep_trace_loop_events(control_context, config)
            control_metrics[control_name] = selector_metrics(
                control_context,
                selector_name=f"{config.selector_config}_{control_name}",
                selected_events=control_selected,
            )
        control_false_rates = [
            metric.false_nucleus_rate for metric in control_metrics.values()
        ]
        control_precisions = [
            metric.contact_cluster_precision for metric in control_metrics.values()
        ]
        control_enrichments = [
            metric.real_vs_decoy_coupling_enrichment_ratio
            for metric in control_metrics.values()
        ]
        max_control_enrichment = max(control_enrichments) if control_enrichments else 0.0
        enrichment_ratio = (
            _rounded(
                external_metric.real_vs_decoy_coupling_enrichment_ratio
                / max_control_enrichment
            )
            if external_metric.selected_event_count > 0 and max_control_enrichment
            else 0.0
        )
        beats_physical = (
            external_metric.selected_event_count > 0
            and external_metric.false_nucleus_rate < physical_metric.false_nucleus_rate
            and external_metric.contact_cluster_precision
            > physical_metric.contact_cluster_precision
        )
        beats_matched_controls = (
            external_metric.selected_event_count > 0
            and all(
                external_metric.false_nucleus_rate < value
                for value in control_false_rates
            )
            and all(
                external_metric.contact_cluster_precision > value
                for value in control_precisions
            )
        )
        meets_oracle_recall_floor = (
            external_metric.selected_event_count > 0
            and external_metric.long_range_contact_recall >= oracle_recall_floor
        )
        passes_enrichment_gate = enrichment_ratio > 1.25
        row = {
            **asdict(config),
            "external_selected_event_count": external_metric.selected_event_count,
            "external_false_nucleus_rate": external_metric.false_nucleus_rate,
            "external_cluster_precision": external_metric.contact_cluster_precision,
            "external_long_range_recall": external_metric.long_range_contact_recall,
            "external_constraint_recall": external_metric.coupling_constraint_recall,
            "external_real_vs_decoy_enrichment_ratio": (
                external_metric.real_vs_decoy_coupling_enrichment_ratio
            ),
            "external_vs_max_control_enrichment_ratio": enrichment_ratio,
            "physical_false_nucleus_rate": physical_metric.false_nucleus_rate,
            "physical_cluster_precision": physical_metric.contact_cluster_precision,
            "oracle_recall_floor": oracle_recall_floor,
            "beats_physical": beats_physical,
            "beats_matched_controls": beats_matched_controls,
            "meets_oracle_recall_floor": meets_oracle_recall_floor,
            "passes_enrichment_gate": passes_enrichment_gate,
            "external_probe_passed": (
                beats_physical
                and beats_matched_controls
                and meets_oracle_recall_floor
                and passes_enrichment_gate
            ),
        }
        for control_name, metric in control_metrics.items():
            row[_control_field("control_selected", control_name)] = (
                metric.selected_event_count
            )
            row[_control_field("control_false", control_name)] = (
                metric.false_nucleus_rate
            )
            row[_control_field("control_precision", control_name)] = (
                metric.contact_cluster_precision
            )
            row[_control_field("control_recall", control_name)] = (
                metric.long_range_contact_recall
            )
            row[_control_field("control_enrichment", control_name)] = (
                metric.real_vs_decoy_coupling_enrichment_ratio
            )
        rows_out.append(row)

    rows_out.sort(
        key=lambda row: (
            not bool(row["external_probe_passed"]),
            not bool(row["beats_matched_controls"]),
            not bool(row["beats_physical"]),
            not bool(row["meets_oracle_recall_floor"]),
            not bool(row["passes_enrichment_gate"]),
            float(row["external_false_nucleus_rate"]),
            -float(row["external_cluster_precision"]),
            -float(row["external_vs_max_control_enrichment_ratio"]),
        )
    )
    write_csv_rows(rows_out, output)
    if top_results_output is not None:
        top_results_output.parent.mkdir(parents=True, exist_ok=True)
        top_results_output.write_text(
            json.dumps(rows_out[:25], ensure_ascii=False, indent=2, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
    print(f"wrote {len(rows_out)} selector configs to {output}")
    for row in rows_out[:10]:
        print(
            row["selector_config"],
            "pass=",
            row["external_probe_passed"],
            "beats_controls=",
            row["beats_matched_controls"],
            "false=",
            row["external_false_nucleus_rate"],
            "precision=",
            row["external_cluster_precision"],
            "recall=",
            row["external_long_range_recall"],
            "enrichment_vs_control=",
            row["external_vs_max_control_enrichment_ratio"],
        )
    return output, top_results_output


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep external-only trace-loop selector gates against matched "
            "negative controls. This is exploratory and does not unlock claims."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--external-coupling-file",
        default=str(DEFAULT_EXTERNAL_COUPLING_FILE),
    )
    parser.add_argument(
        "--oracle-coupling-file",
        default=str(DEFAULT_ORACLE_COUPLING_FILE),
    )
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--top-results-output")
    parser.add_argument(
        "--max-configs",
        type=int,
        help="Limit the deterministic config grid for quick smoke runs.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
        help="Print progress to stderr every N selector configs; use 0 to silence.",
    )
    args = parser.parse_args()
    run_selector_sweep(
        benchmark_file=Path(args.benchmark_file),
        external_coupling_file=Path(args.external_coupling_file),
        oracle_coupling_file=Path(args.oracle_coupling_file),
        output=Path(args.output),
        top_results_output=(
            Path(args.top_results_output) if args.top_results_output else None
        ),
        max_configs=args.max_configs,
        progress_every=args.progress_every,
    )


if __name__ == "__main__":
    main()
