import gzip
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.external_dca_runner import (  # noqa: E402
    FAMILY_WIDE_PFAM_APC_METHOD,
    PfamDomainMapping,
    read_stockholm_sample,
    run_pfam_apc_covariation_for_row,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


BENCHMARK_8 = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"


def _write_stockholm(path: Path) -> Path:
    rows = ["# STOCKHOLM 1.0"]
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    for index in range(40):
        left = alphabet[index % len(alphabet)]
        right = alphabet[(index * 3) % len(alphabet)]
        sequence = f"{left}CDEFGHIK{right}MNPQRSTVWY"
        rows.append(f"seq{index:03d} {sequence}")
    rows.append("//")
    with gzip.open(path, "wt", encoding="utf-8") as file:
        file.write("\n".join(rows) + "\n")
    return path


def test_pfam_apc_runner_scores_gzipped_stockholm_without_truth(tmp_path) -> None:
    row = load_real_coordinate_visual_rows(BENCHMARK_8)[0]
    alignment = _write_stockholm(tmp_path / "PFTEST.full.sto.gz")

    sequences, total_seen = read_stockholm_sample(alignment, max_records=25)
    result = run_pfam_apc_covariation_for_row(
        row=row,
        mappings=(
            PfamDomainMapping(
                pfam_id="PFTEST",
                name="Synthetic",
                description="Synthetic external alignment",
                start=1,
                end=20,
                coverage=1.0,
            ),
        ),
        alignment_dir=tmp_path,
        max_records=25,
        minimum_sequence_separation=5,
    )

    assert len(sequences) == 25
    assert total_seen == 40
    assert result.dca_status == "dca_available"
    assert result.covariation_method == FAMILY_WIDE_PFAM_APC_METHOD
    assert result.sample_depth == 25
    assert result.accepted_pair_count > 0
    assert result.pairs[0].i < result.pairs[0].j
