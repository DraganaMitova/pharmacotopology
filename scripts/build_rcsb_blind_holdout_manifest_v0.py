from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from difflib import SequenceMatcher
from pathlib import Path
from typing import Mapping, Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_contact_topology import sha256_sequence  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_topology import normalize_sequence  # noqa: E402


RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_POLYMER_ENTITY_URL = (
    "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
)
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
PDBE_PFAM_MAPPING_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb_id}"
DEFAULT_OUTPUT = Path("data/blind_holdout_manifest_v0.locked.json")
DEFAULT_EXCLUDE_BENCHMARKS = (
    Path("data/folding_real_coordinate_visual_8.locked.json"),
    Path("data/folding_real_coordinate_holdout_1ubq.locked.json"),
)
STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _sha256_text(text: str) -> str:
    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def _selection_hash(seed: str, identifier: str) -> str:
    return hashlib.sha256(f"{seed}:{identifier}".encode("utf-8")).hexdigest()


def _download_json(url: str, *, timeout: int = 60) -> Mapping[str, object]:
    with urllib.request.urlopen(url, timeout=timeout) as response:
        parsed = json.loads(response.read().decode("utf-8"))
    if not isinstance(parsed, Mapping):
        raise ValueError(f"{url} did not return a JSON object")
    return parsed


def _post_json(url: str, payload: Mapping[str, object], *, timeout: int = 60) -> object:
    body = json.dumps(payload, sort_keys=True).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=body,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def _search_query(*, rows: int) -> dict[str, object]:
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_entity_polymer_type",
                        "operator": "exact_match",
                        "value": "Protein",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": (
                            "rcsb_entry_info.polymer_entity_count_protein"
                        ),
                        "operator": "equals",
                        "value": 1,
                    },
                },
            ],
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "scoring_strategy": "combined",
        },
    }


def _existing_exclusions(
    benchmark_files: Sequence[Path],
) -> tuple[set[str], set[str], tuple[str, ...]]:
    source_accessions: set[str] = set()
    sequence_hashes: set[str] = set()
    sequences: list[str] = []
    for path in benchmark_files:
        if not path.exists():
            continue
        try:
            rows = load_real_coordinate_visual_rows(path)
        except ValueError:
            continue
        for row in rows:
            source_accessions.add(row.source_accession.upper())
            sequence_hashes.add(row.sequence_sha256)
            sequences.append(row.sequence)
    return source_accessions, sequence_hashes, tuple(sequences)


def _too_similar_to_excluded(
    sequence: str,
    excluded_sequences: Sequence[str],
    *,
    max_similarity: float,
) -> bool:
    for excluded in excluded_sequences:
        ratio = SequenceMatcher(
            None,
            sequence,
            excluded,
            autojunk=False,
        ).ratio()
        if ratio >= max_similarity:
            return True
    return False


def _identifier_parts(identifier: str) -> tuple[str, str]:
    if "_" not in identifier:
        raise ValueError(f"invalid RCSB polymer_entity identifier: {identifier}")
    pdb_id, entity_id = identifier.split("_", 1)
    return pdb_id.upper(), entity_id


def _resolution(entry: Mapping[str, object]) -> Optional[float]:
    info = entry.get("rcsb_entry_info", {})
    if not isinstance(info, Mapping):
        return None
    values = info.get("resolution_combined", [])
    if not isinstance(values, list) or not values:
        return None
    try:
        return float(min(float(value) for value in values))
    except (TypeError, ValueError):
        return None


def _method(entry: Mapping[str, object]) -> str:
    methods = entry.get("exptl", [])
    if not isinstance(methods, list) or not methods:
        return ""
    first = methods[0]
    if not isinstance(first, Mapping):
        return ""
    return str(first.get("method", "")).upper()


def _pfam_ids_for_chain(pdb_id: str, chain_id: str) -> tuple[str, ...]:
    parsed = _download_json(PDBE_PFAM_MAPPING_URL.format(pdb_id=pdb_id.lower()))
    entry = parsed.get(pdb_id.lower(), {})
    if not isinstance(entry, Mapping):
        return ()
    pfams = entry.get("Pfam", {})
    if not isinstance(pfams, Mapping):
        return ()
    hits: list[str] = []
    for pfam_id, info in pfams.items():
        if not isinstance(info, Mapping):
            continue
        mappings = info.get("mappings", [])
        if not isinstance(mappings, list):
            continue
        for mapping in mappings:
            if not isinstance(mapping, Mapping):
                continue
            if (
                mapping.get("chain_id") == chain_id
                or mapping.get("struct_asym_id") == chain_id
            ):
                hits.append(str(pfam_id))
                break
    return tuple(sorted(set(hits)))


def _candidate_target(
    *,
    identifier: str,
    selection_rank: int,
    selection_hash: str,
    min_length: int,
    max_length: int,
    max_resolution: float,
    allowed_methods: set[str],
    excluded_source_accessions: set[str],
    excluded_sequence_hashes: set[str],
    excluded_sequences: Sequence[str],
    max_excluded_sequence_similarity: float,
    require_pfam_mapping: bool,
    benchmark_dir: Path,
    coupling_dir: Path,
) -> dict[str, object] | None:
    pdb_id, entity_id = _identifier_parts(identifier)
    entity = _download_json(
        RCSB_POLYMER_ENTITY_URL.format(pdb_id=pdb_id, entity_id=entity_id)
    )
    entry = _download_json(RCSB_ENTRY_URL.format(pdb_id=pdb_id))
    identifiers = entity.get("rcsb_polymer_entity_container_identifiers", {})
    if not isinstance(identifiers, Mapping):
        return None
    chain_ids = tuple(str(chain) for chain in identifiers.get("auth_asym_ids", ()))
    if len(chain_ids) != 1:
        return None
    chain_id = chain_ids[0]
    source_accession = f"{pdb_id}:{chain_id}".upper()
    if source_accession in excluded_source_accessions:
        return None

    entity_poly = entity.get("entity_poly", {})
    if not isinstance(entity_poly, Mapping):
        return None
    if str(entity_poly.get("rcsb_entity_polymer_type", "")) != "Protein":
        return None
    sequence = normalize_sequence(
        str(entity_poly.get("pdbx_seq_one_letter_code_can", ""))
    )
    if not set(sequence) <= STANDARD_AA:
        return None
    if not min_length <= len(sequence) <= max_length:
        return None
    sequence_hash = sha256_sequence(sequence)
    if sequence_hash in excluded_sequence_hashes:
        return None
    if _too_similar_to_excluded(
        sequence,
        excluded_sequences,
        max_similarity=max_excluded_sequence_similarity,
    ):
        return None

    method = _method(entry)
    if method not in allowed_methods:
        return None
    resolution = _resolution(entry)
    if resolution is None or resolution > max_resolution:
        return None

    pfam_ids: tuple[str, ...] = ()
    if require_pfam_mapping:
        try:
            pfam_ids = _pfam_ids_for_chain(pdb_id, chain_id)
        except (urllib.error.URLError, TimeoutError, ValueError):
            return None
        if not pfam_ids:
            return None

    slug = f"blind_{selection_rank:03d}_{pdb_id.lower()}_{chain_id.lower()}"
    accession_info = entry.get("rcsb_accession_info", {})
    release_date = ""
    if isinstance(accession_info, Mapping):
        release_date = str(accession_info.get("initial_release_date", ""))
    uniprot_ids = identifiers.get("uniprot_ids", ())
    if not isinstance(uniprot_ids, list):
        uniprot_ids = []
    return {
        "target_id": slug,
        "selection_rank": selection_rank,
        "selection_hash": selection_hash,
        "pdb_id": pdb_id,
        "entity_id": entity_id,
        "chain_id": chain_id,
        "source_accession": source_accession,
        "sequence_length": len(sequence),
        "sequence_sha256": sequence_hash,
        "source_database": "RCSB_PDB",
        "experimental_method": method,
        "resolution_angstrom": round(resolution, 3),
        "initial_release_date": release_date,
        "uniprot_ids": tuple(str(item) for item in uniprot_ids),
        "pfam_ids": pfam_ids,
        "pfam_mapping_required": require_pfam_mapping,
        "coordinate_source_url": f"https://files.rcsb.org/download/{pdb_id}.pdb",
        "rcsb_polymer_entity_url": RCSB_POLYMER_ENTITY_URL.format(
            pdb_id=pdb_id,
            entity_id=entity_id,
        ),
        "rcsb_entry_url": RCSB_ENTRY_URL.format(pdb_id=pdb_id),
        "benchmark_file": str(benchmark_dir / f"{slug}.locked.json"),
        "external_coupling_file": str(
            coupling_dir / f"{slug}_query_centered_apc.locked.json"
        ),
        "raw_sequence_exposed": False,
        "coordinate_truth_used_to_build_constraints": False,
        "native_truth_used_before_coupling_selection": False,
        "oracle_constraint_control": False,
    }


def build_rcsb_blind_holdout_manifest_v0(
    *,
    output: Path,
    max_targets: int,
    search_rows: int,
    seed: str,
    min_length: int,
    max_length: int,
    max_resolution: float,
    allowed_methods: set[str],
    exclude_benchmark_files: Sequence[Path],
    max_excluded_sequence_similarity: float,
    require_pfam_mapping: bool,
    benchmark_dir: Path,
    coupling_dir: Path,
) -> Path:
    query = _search_query(rows=search_rows)
    response = _post_json(RCSB_SEARCH_URL, query)
    if not isinstance(response, Mapping):
        raise ValueError("RCSB search response must be a JSON object")
    result_set = response.get("result_set", [])
    if not isinstance(result_set, list):
        raise ValueError("RCSB search response result_set must be a list")
    identifiers = [
        str(item.get("identifier", ""))
        for item in result_set
        if isinstance(item, Mapping) and item.get("identifier")
    ]
    ranked_identifiers = sorted(
        identifiers,
        key=lambda identifier: _selection_hash(seed, identifier),
    )
    excluded_accessions, excluded_hashes, excluded_sequences = _existing_exclusions(
        exclude_benchmark_files
    )
    targets: list[dict[str, object]] = []
    rejection_counts: dict[str, int] = {}
    for identifier in ranked_identifiers:
        if len(targets) >= max_targets:
            break
        selection_rank = len(targets) + 1
        try:
            target = _candidate_target(
                identifier=identifier,
                selection_rank=selection_rank,
                selection_hash=_selection_hash(seed, identifier),
                min_length=min_length,
                max_length=max_length,
                max_resolution=max_resolution,
                allowed_methods=allowed_methods,
                excluded_source_accessions=excluded_accessions,
                excluded_sequence_hashes=excluded_hashes,
                excluded_sequences=excluded_sequences,
                max_excluded_sequence_similarity=max_excluded_sequence_similarity,
                require_pfam_mapping=require_pfam_mapping,
                benchmark_dir=benchmark_dir,
                coupling_dir=coupling_dir,
            )
        except (ValueError, urllib.error.URLError, TimeoutError) as exc:
            rejection_counts[exc.__class__.__name__] = (
                rejection_counts.get(exc.__class__.__name__, 0) + 1
            )
            continue
        if target is None:
            rejection_counts["filtered"] = rejection_counts.get("filtered", 0) + 1
            continue
        targets.append(target)

    payload = {
        "manifest_kind": "rcsb_blind_external_holdout_manifest_v0",
        "manifest_frozen": True,
        "selector_frozen_after_manifest": True,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "seed": seed,
        "target_count": len(targets),
        "requested_target_count": max_targets,
        "search_rows": search_rows,
        "source_api": {
            "rcsb_search_url": RCSB_SEARCH_URL,
            "rcsb_polymer_entity_url": RCSB_POLYMER_ENTITY_URL,
            "rcsb_entry_url": RCSB_ENTRY_URL,
            "pdbe_pfam_mapping_url": PDBE_PFAM_MAPPING_URL,
        },
        "selection_policy": {
            "return_type": "polymer_entity",
            "single_chain_only": True,
            "protein_entity_count_per_entry": 1,
            "min_length": min_length,
            "max_length": max_length,
            "max_resolution_angstrom": max_resolution,
            "allowed_methods": tuple(sorted(allowed_methods)),
            "require_pfam_mapping": require_pfam_mapping,
            "exclude_benchmark_files": tuple(str(path) for path in exclude_benchmark_files),
            "excluded_source_accession_count": len(excluded_accessions),
            "excluded_sequence_sha256_count": len(excluded_hashes),
            "excluded_sequence_similarity_count": len(excluded_sequences),
            "max_excluded_sequence_similarity": (
                max_excluded_sequence_similarity
            ),
            "deterministic_hash_selection": True,
            "manual_target_selection_after_scoring": False,
        },
        "rcsb_search_query_sha256": _sha256_text(
            json.dumps(query, sort_keys=True)
        ),
        "rejection_counts": rejection_counts,
        "raw_sequence_exposed": False,
        "coordinate_truth_used_to_build_constraints": False,
        "native_truth_used_before_coupling_selection": False,
        "oracle_constraint_control": False,
        "targets": targets,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a frozen RCSB-derived blind holdout manifest for external "
            "folding-nucleus batteries. This writes target metadata only; it "
            "does not build coordinate locks or coupling files."
        )
    )
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--max-targets", type=int, default=50)
    parser.add_argument("--search-rows", type=int, default=5000)
    parser.add_argument("--seed", default="blind_holdout_manifest_v0")
    parser.add_argument("--min-length", type=int, default=50)
    parser.add_argument("--max-length", type=int, default=250)
    parser.add_argument("--max-resolution", type=float, default=2.5)
    parser.add_argument("--max-excluded-sequence-similarity", type=float, default=0.35)
    parser.add_argument(
        "--allowed-method",
        action="append",
        default=["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"],
    )
    parser.add_argument(
        "--exclude-benchmark-file",
        action="append",
        default=[str(path) for path in DEFAULT_EXCLUDE_BENCHMARKS],
    )
    parser.add_argument(
        "--allow-without-pfam-mapping",
        action="store_true",
    )
    parser.add_argument(
        "--benchmark-dir",
        default="data/blind_coordinate_holdouts_v0",
    )
    parser.add_argument(
        "--coupling-dir",
        default="data/blind_external_couplings_v0",
    )
    args = parser.parse_args()
    output = build_rcsb_blind_holdout_manifest_v0(
        output=Path(args.output),
        max_targets=args.max_targets,
        search_rows=args.search_rows,
        seed=args.seed,
        min_length=args.min_length,
        max_length=args.max_length,
        max_resolution=args.max_resolution,
        allowed_methods={method.upper() for method in args.allowed_method},
        exclude_benchmark_files=tuple(
            Path(path) for path in args.exclude_benchmark_file
        ),
        max_excluded_sequence_similarity=args.max_excluded_sequence_similarity,
        require_pfam_mapping=not args.allow_without_pfam_mapping,
        benchmark_dir=Path(args.benchmark_dir),
        coupling_dir=Path(args.coupling_dir),
    )
    print(output)


if __name__ == "__main__":
    main()
