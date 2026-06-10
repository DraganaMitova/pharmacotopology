from pathlib import Path

from pharmacotopology.folding_contact_law_features import contact_law_feature_rows_for_row
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_iterative_distance_geometry_diffusion import (
    ITERATIVE_DISTANCE_GEOMETRY_DIFFUSION_KIND,
    run_iterative_distance_geometry_diffusion,
)
from pharmacotopology.folding_nucleus_closure_search import nucleus_closure_events_for_row
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_iterative_distance_geometry_diffusion_is_native_free_and_bounded() -> None:
    rows = load_real_coordinate_visual_rows(REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json")
    row = next(item for item in rows if item.source_accession == "1PGA:A")
    coupling_dataset = load_coupling_dataset(
        REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
    )
    constraints = tuple(item for item in coupling_dataset.constraints if item.source_accession == row.source_accession)
    features = contact_law_feature_rows_for_row(row)
    events = nucleus_closure_events_for_row(row, features)

    packet = run_iterative_distance_geometry_diffusion(
        row=row,
        constraints=constraints,
        events=events,
        iteration_count=1,
        attraction_budget_fraction=0.9,
        final_budget_fraction=1.1,
        max_degree=5,
    )

    assert packet.kind == ITERATIVE_DISTANCE_GEOMETRY_DIFFUSION_KIND
    assert packet.final_pair_count > 0
    assert packet.iteration_count == 1
    assert not packet.coordinate_truth_used_before_selection
    assert not packet.native_truth_used_before_selection
    assert not packet.alphafold_used_before_selection
    assert not packet.structure_template_used_before_selection
    assert not packet.raw_sequence_exposed
    assert all(not row_score.native_truth_used_before_selection for row_score in packet.scores)
