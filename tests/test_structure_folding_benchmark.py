import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.folding_structure_benchmark import (  # noqa: E402
    STRUCTURE_BENCHMARK_KIND,
    build_falsification_report,
    build_structure_benchmark_report,
    load_real10_with_structure_evidence,
    parse_pdb_contact_graph,
    sequence_order_control_rows,
    structure_signature_from_contact_graph,
)


BENCHMARK_FILE = ROOT / "data" / "folding_benchmarks_real_10.locked.json"
STRUCTURE_EVIDENCE_FILE = (
    ROOT / "data" / "folding_benchmarks_real_10_structure_evidence.json"
)


PDB_FIXTURE = """\
HELIX    1   A ALA A    1  LEU A    4  1                                   4
SHEET    1   A 2 VAL A   7  THR A   9  0
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C
ATOM      2  CA  LEU A   2       1.500   0.000   0.000  1.00 10.00           C
ATOM      3  CA  GLU A   3       3.000   0.000   0.000  1.00 10.00           C
ATOM      4  CA  LEU A   4       4.500   0.000   0.000  1.00 10.00           C
ATOM      5  CA  GLY A   5       6.000   0.000   0.000  1.00 10.00           C
ATOM      6  CA  LYS A   6       7.500   0.000   0.000  1.00 10.00           C
ATOM      7  CA  VAL A   7       0.000   4.000   0.000  1.00 10.00           C
ATOM      8  CA  ILE A   8       1.500   4.000   0.000  1.00 10.00           C
ATOM      9  CA  THR A   9       3.000   4.000   0.000  1.00 10.00           C
ATOM     10  CA  PHE A  10       4.500   4.000   0.000  1.00 10.00           C
END
"""


def test_pdb_contact_graph_extractor_reads_structure_signals() -> None:
    graph = parse_pdb_contact_graph(PDB_FIXTURE, chain_id="A")
    signature, features = structure_signature_from_contact_graph(
        graph,
        expected_sequence_length=10,
        reference_sequence="ALELGKVITF",
    )

    assert len(graph.residues) == 10
    assert graph.contact_pairs
    assert features["contact_count"] > 0
    assert features["helix_fraction"] == 0.4
    assert features["sheet_fraction"] == 0.3
    assert 0.0 <= signature.contact_map_closure <= 1.0
    assert 0.0 <= signature.long_range_contact_order <= 1.0


def test_tracked_structure_evidence_keeps_truth_channels_separate() -> None:
    payload = json.loads(STRUCTURE_EVIDENCE_FILE.read_text(encoding="utf-8"))
    rows = payload["references"]

    assert payload["benchmark_kind"] == STRUCTURE_BENCHMARK_KIND
    assert payload["coordinate_rows"] == 8
    assert payload["disorder_reference_rows"] == 2
    assert payload["folding_problem_solved"] is False
    assert payload["folding_solution_claim_created"] is False
    assert len(rows) == 10
    assert sum(row["evidence_kind"] == "coordinate_contact_graph" for row in rows) == 8
    assert sum(row["evidence_kind"] == "disorder_reference" for row in rows) == 2
    assert all(
        row["structure_topology_signature_kind"]
        != "locked_broad_class_prototype_from_external_label"
        for row in rows
    )


def test_structure_report_exposes_revision_and_order_sensitivity_gap() -> None:
    references, evidence = load_real10_with_structure_evidence(
        BENCHMARK_FILE,
        STRUCTURE_EVIDENCE_FILE,
    )
    report = build_structure_benchmark_report(
        references,
        evidence,
        source_benchmark_file=BENCHMARK_FILE,
        structure_evidence_file=STRUCTURE_EVIDENCE_FILE,
    )
    control_rows = sequence_order_control_rows(references)
    falsification = build_falsification_report(report, control_rows)

    assert report["benchmark_kind"] == STRUCTURE_BENCHMARK_KIND
    assert report["source_label_benchmark_kind"] == "real_external_label_benchmark_v0"
    assert report["prediction_vs_structure_accuracy"] == 0.3
    assert report["prediction_vs_label_accuracy"] == 0.4
    assert report["structure_vs_label_agreement_rate"] == 0.9
    assert report["sequence_order_sensitivity_score"] == 0.0
    assert report["composition_only_warning"] is True
    assert report["revision_required"] is True
    assert report["claim_allowed"] is False
    assert report["folding_problem_solved"] is False

    assert len(control_rows) == 50
    assert all(row["same_composition"] for row in control_rows)
    assert all(row["control_sequence_written"] is False for row in control_rows)
    assert all("control_sequence" not in row for row in control_rows)
    assert falsification["tau_transfer_falsification"] is False
    assert falsification["rho_revision_needed"] is True
    assert falsification["omega_terminal_claim_allowed"] is False
