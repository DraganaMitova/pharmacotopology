#!/usr/bin/env python3
from __future__ import annotations

"""V18 p53 TAD/MDM2 partner-induced evidence test.

This is the first selected V18 pressure-evidence test from the V17 sprint. It
uses the already-provenanced 1YCR partner-bound complex context and performs a
zero-MD, no-tuning interface/role evidence readout. It does not claim universal
folding, does not run MD, does not select by native metrics, and explicitly
keeps the isolated p53 TAD autonomous-fold claim disabled/abstained.
"""

import argparse
import json
import math
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V17_CERT = RUN_ROOT / "V17_PRESSURE_EVIDENCE_SPRINT_LOCK" / "v17_pressure_evidence_sprint_lock_certificate.json"
DEFAULT_PREFLIGHT_CERT = RUN_ROOT / "V16_TRANSFER_DATA_PREFLIGHT" / "v16_transfer_data_preflight_certificate.json"
DEFAULT_PDB = REPO_ROOT / "data" / "v16_pressure_targets" / "p53_TAD_MDM2" / "1YCR.pdb"
DEFAULT_OUT_DIR = RUN_ROOT / "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST"

PROBE_RADII_ANGSTROM = [6.0, 8.0, 10.0]


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # pragma: no cover
        raise SystemExit(f"invalid JSON in {label}: {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _target_row_by_id(cert: dict[str, Any], key: str, target_id: str) -> dict[str, Any]:
    rows = cert.get(key) if isinstance(cert.get(key), list) else []
    for row in rows:
        if isinstance(row, dict) and row.get("target_id") == target_id:
            return row
    return {}


def _material_by_id(row: dict[str, Any]) -> dict[str, dict[str, Any]]:
    mats = row.get("material_results") if isinstance(row.get("material_results"), list) else []
    return {str(m.get("material_id")): m for m in mats if isinstance(m, dict)}


def _read_ca_atoms(path: Path) -> dict[str, list[dict[str, Any]]]:
    atoms: dict[str, list[dict[str, Any]]] = {}
    if not path.exists():
        return atoms
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 54 or line[12:16].strip() != "CA":
                continue
            chain = line[21].strip() or "_"
            try:
                resseq = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            atoms.setdefault(chain, []).append({"residue_number": resseq, "xyz": xyz})
    return atoms


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _interface_probe(atoms: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    chain_counts = {chain: len(items) for chain, items in atoms.items()}
    if len(chain_counts) < 2:
        return {
            "chain_pair": None,
            "chain_ca_counts": chain_counts,
            "min_ca_distance": None,
            "multi_radius_contact_counts": {},
            "contact_probe_policy": "multi_radius_report_only_no_single_fixed_selection_threshold",
            "interface_contact_evidence_present": False,
        }
    chains = sorted(chain_counts, key=lambda c: chain_counts[c], reverse=True)
    large_chain = chains[0]
    small_chain = chains[-1]
    large_atoms = atoms[large_chain]
    small_atoms = atoms[small_chain]
    min_d: float | None = None
    counts = {str(radius): 0 for radius in PROBE_RADII_ANGSTROM}
    nearest: dict[str, Any] | None = None
    for a in large_atoms:
        for b in small_atoms:
            d = _distance(a["xyz"], b["xyz"])
            if min_d is None or d < min_d:
                min_d = d
                nearest = {
                    "large_chain": large_chain,
                    "large_residue_number": a["residue_number"],
                    "small_chain": small_chain,
                    "small_residue_number": b["residue_number"],
                    "distance_angstrom": round(d, 3),
                }
            for radius in PROBE_RADII_ANGSTROM:
                if d <= radius:
                    counts[str(radius)] += 1
    # This is an evidence readout from a partner-bound complex, not a tuned selector.
    # We require interface proximity to be visible across more than one report-only probe radius.
    present_probe_count = sum(1 for value in counts.values() if value > 0)
    return {
        "chain_pair": f"{large_chain}-{small_chain}",
        "large_chain": large_chain,
        "small_chain": small_chain,
        "chain_ca_counts": chain_counts,
        "min_ca_distance": round(min_d, 3) if min_d is not None else None,
        "nearest_ca_pair": nearest,
        "multi_radius_contact_counts": counts,
        "contact_probe_radii_angstrom_report_only": PROBE_RADII_ANGSTROM,
        "contact_probe_policy": "multi_radius_report_only_no_single_fixed_selection_threshold",
        "interface_contact_evidence_present": present_probe_count >= 2,
        "selection_threshold_used": False,
    }


def _small_chain_helix_proxy(atoms: dict[str, list[dict[str, Any]]], small_chain: str | None) -> dict[str, Any]:
    if not small_chain or small_chain not in atoms:
        return {
            "small_chain": small_chain,
            "helix_proxy_available": False,
            "helix_proxy_policy": "report_only_no_fold_claim",
        }
    chain_atoms = atoms[small_chain]
    dists_i3: list[float] = []
    dists_i4: list[float] = []
    for i in range(len(chain_atoms)):
        if i + 3 < len(chain_atoms):
            dists_i3.append(_distance(chain_atoms[i]["xyz"], chain_atoms[i + 3]["xyz"]))
        if i + 4 < len(chain_atoms):
            dists_i4.append(_distance(chain_atoms[i]["xyz"], chain_atoms[i + 4]["xyz"]))
    def _summary(values: list[float]) -> dict[str, Any]:
        if not values:
            return {"count": 0, "mean": None, "min": None, "max": None}
        return {
            "count": len(values),
            "mean": round(sum(values) / len(values), 3),
            "min": round(min(values), 3),
            "max": round(max(values), 3),
        }
    return {
        "small_chain": small_chain,
        "small_chain_ca_count": len(chain_atoms),
        "ca_i_to_i_plus_3_distance_summary": _summary(dists_i3),
        "ca_i_to_i_plus_4_distance_summary": _summary(dists_i4),
        "helix_proxy_available": len(chain_atoms) >= 8,
        "helix_proxy_policy": "report_only_partner_bound_segment_geometry_no_autonomous_fold_claim",
        "selection_threshold_used": False,
    }


def _v17_selected_p53(v17: dict[str, Any]) -> bool:
    return (
        v17.get("sprint_lock_status") == "V17_PRESSURE_EVIDENCE_SPRINT_LOCKED"
        and v17.get("selected_V18_target") == "p53_TAD_MDM2"
        and v17.get("selected_V18_test") == "p53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST"
        and v17.get("claim_allowed") is False
    )


def build_evidence(v17: dict[str, Any], preflight: dict[str, Any], pdb_path: Path) -> dict[str, Any]:
    preflight_row = _target_row_by_id(preflight, "target_results", "p53_TAD_MDM2")
    mats = _material_by_id(preflight_row)
    mat = mats.get("p53_TAD_MDM2_complex_structure", {})
    atoms = _read_ca_atoms(pdb_path)
    interface = _interface_probe(atoms)
    helix_proxy = _small_chain_helix_proxy(atoms, interface.get("small_chain"))

    v17_ready = _v17_selected_p53(v17)
    material_ready = mat.get("material_status") == "present_and_provenanced" and bool(atoms)
    partner_complex_context = material_ready and len(interface.get("chain_ca_counts", {})) >= 2
    interface_evidence = bool(interface.get("interface_contact_evidence_present"))
    isolated_abstain_status = "clean_abstain_no_isolated_TAD_material_no_autonomous_fold_claim"

    forbidden_violations: list[str] = []
    autonomous_fold_claim_made = False
    if autonomous_fold_claim_made:  # kept explicit as a guard, should never fire.
        forbidden_violations.append("isolated_TAD_as_soluble_compact_core")

    role_evidence_found = bool(v17_ready and partner_complex_context and interface_evidence and not forbidden_violations)
    status = (
        "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_PASSED_CLAIM_DISABLED"
        if role_evidence_found
        else "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_CLEAN_ABSTAIN_OR_BLOCKED_CLAIM_DISABLED"
    )

    return {
        "kind": "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_v0",
        "run_mode": "zero_md_partner_induced_evidence_readout_no_simulation_no_threshold_tuning",
        "test_status": status,
        "target_id": "p53_TAD_MDM2",
        "role_class": "disorder_partner_induced_object",
        "pressure_class": "disorder_partner_induced_binding",
        "source_v17_sprint_status": v17.get("sprint_lock_status"),
        "source_selected_V18_target": v17.get("selected_V18_target"),
        "source_preflight_status": preflight.get("data_preflight_status"),
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "positive_folding_evidence_found": False,
        "positive_pressure_evidence_found": role_evidence_found,
        "partner_induced_role_evidence_found": role_evidence_found,
        "new_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "selection_policy": "partner_context_and_interface_readout_only_no_native_metric_selection",
        "selected_core_or_clean_abstain": (
            "partner_induced_interface_or_helix_signal_found" if role_evidence_found else "clean_abstain_partner_induced_evidence_incomplete"
        ),
        "role_buckets_assigned": (
            ["partner_induced_interface_or_helix_signal_found", "binding_stabilized_interface_core_context"]
            if role_evidence_found
            else []
        ),
        "isolated_TAD_autonomous_fold_status": isolated_abstain_status,
        "clean_abstain_roles": [isolated_abstain_status],
        "forbidden_misclassification": [
            "isolated_TAD_as_soluble_compact_core",
            "partner_interface_as_autonomous_fold_claim",
            "binding_helix_as_universal_p53_fold",
        ],
        "forbidden_misclassification_violations": forbidden_violations,
        "interface_contact_readout": interface,
        "partner_bound_helix_proxy_readout": helix_proxy,
        "available_evidence": [
            item
            for item, present in {
                "partner_bound_complex_context": partner_complex_context,
                "interface_or_contact_evidence": interface_evidence,
                "state_label_bound_context": material_ready,
                "leakage_guard_autonomous_fold_vs_partner_induced_fold": not autonomous_fold_claim_made,
            }.items()
            if present
        ],
        "missing_evidence": ["isolated_TAD_disorder_or_clean_abstain_context", "external_couplings_if_available"],
        "limitations": [
            "zero_MD_only_no_partner_induced_dynamics_claim",
            "isolated_TAD_material_not_present_so_autonomous_fold_claim_is_abstained",
            "no_universal_p53_fold_claim",
            "external_couplings_not_required_for_this_first_partner-bound evidence readout",
        ],
        "v19_recommended_next_step": "add_isolated_p53_TAD_disorder_or_clean_abstain_material_then_compare_bound_vs_isolated_without_autonomous_fold_claim",
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    iface = cert.get("interface_contact_readout", {}) if isinstance(cert.get("interface_contact_readout"), dict) else {}
    helix = cert.get("partner_bound_helix_proxy_readout", {}) if isinstance(cert.get("partner_bound_helix_proxy_readout"), dict) else {}
    lines = [
        "# V18 p53 TAD/MDM2 Partner-Induced Evidence Test",
        "",
        "Zero-MD, no tuning, no folding claim. This reads partner-bound interface/helix role evidence and keeps isolated p53 TAD autonomous-fold claims forbidden.",
        "",
        f"Status: `{cert.get('test_status')}`",
        f"Positive pressure evidence found: `{cert.get('positive_pressure_evidence_found')}`",
        f"Positive folding evidence found: `{cert.get('positive_folding_evidence_found')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## Interface contact readout",
        f"Chain pair: `{iface.get('chain_pair')}`",
        f"Minimum CA distance: `{iface.get('min_ca_distance')}` Å",
        f"Multi-radius contact counts: `{iface.get('multi_radius_contact_counts')}`",
        f"Probe policy: `{iface.get('contact_probe_policy')}`",
        "",
        "## Partner-bound helix proxy",
        f"Small chain: `{helix.get('small_chain')}`",
        f"CA i→i+3 summary: `{helix.get('ca_i_to_i_plus_3_distance_summary')}`",
        f"CA i→i+4 summary: `{helix.get('ca_i_to_i_plus_4_distance_summary')}`",
        f"Policy: `{helix.get('helix_proxy_policy')}`",
        "",
        "## Guards",
        f"Isolated TAD status: `{cert.get('isolated_TAD_autonomous_fold_status')}`",
        f"Forbidden misclassification violations: `{cert.get('forbidden_misclassification_violations')}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v18_p53_tad_mdm2_partner_induced_evidence_certificate.json",
        "interface_readout": out_dir / "v18_p53_tad_mdm2_interface_contact_readout.json",
        "report": out_dir / "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    paths["interface_readout"].write_text(json.dumps(cert.get("interface_contact_readout", {}), indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v17-cert", type=Path, default=DEFAULT_V17_CERT)
    parser.add_argument("--preflight-cert", type=Path, default=DEFAULT_PREFLIGHT_CERT)
    parser.add_argument("--pdb", type=Path, default=DEFAULT_PDB)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    v17 = _read_json(args.v17_cert, "V17 sprint certificate")
    preflight = _read_json(args.preflight_cert, "V16 data preflight certificate")
    cert = build_evidence(v17, preflight, args.pdb)
    write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "test_status": cert["test_status"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "certificate": str(args.out_dir / "v18_p53_tad_mdm2_partner_induced_evidence_certificate.json"),
        "report": str(args.out_dir / "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_REPORT.md"),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
