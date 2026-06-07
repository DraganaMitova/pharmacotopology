from __future__ import annotations

"""Bounded internet-facing MSA-free structure predictor adapters.

The project cannot assume that ESMFold/OmegaFold/SPIRED are installed locally.
This module provides narrow, auditable adapters that can try public or
key-backed single-sequence services when the runtime has internet access.  The
adapters deliberately do not persist the submitted sequence, do not call
AlphaFold endpoints, and never loop/retry indefinitely.
"""

import json
import socket
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping


ESM_ATLAS_FOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
BIOLM_ESMFOLD_URL = "https://biolm.ai/api/v3/esmfold/predict/"


@dataclass(frozen=True)
class InternetStructurePredictionAttempt:
    predictor_id: str
    url: str
    output_pdb_path: str
    attempted: bool
    success: bool
    status: str
    timeout_seconds: float
    http_status: int | None = None
    pdb_atom_line_count: int = 0
    pdb_ca_line_count: int = 0
    response_bytes: int = 0
    error_tail: str = ""
    sequence_sent_over_network: bool = True
    raw_sequence_persisted: bool = False
    alphafold_endpoint_used: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _tail(text: str, *, limit: int = 1200) -> str:
    return text if len(text) <= limit else text[-limit:]


def _count_pdb_lines(text: str) -> tuple[int, int]:
    atoms = 0
    cas = 0
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atoms += 1
        if line[12:16].strip() == "CA":
            cas += 1
    return atoms, cas


def _looks_like_pdb(text: str, *, minimum_ca_count: int) -> bool:
    atoms, cas = _count_pdb_lines(text)
    return atoms > 0 and cas >= max(1, minimum_ca_count)


def _write_pdb_if_valid(
    *,
    pdb_text: str,
    output_pdb_path: Path,
    minimum_ca_count: int,
) -> tuple[bool, int, int]:
    atom_count, ca_count = _count_pdb_lines(pdb_text)
    if atom_count <= 0 or ca_count < max(1, minimum_ca_count):
        return False, atom_count, ca_count
    output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    output_pdb_path.write_text(pdb_text, encoding="utf-8")
    return True, atom_count, ca_count


def _extract_pdb_strings(obj: Any) -> Iterable[str]:
    """Yield plausible PDB strings from an arbitrary JSON-ish response."""

    if isinstance(obj, str):
        if "ATOM" in obj or "HETATM" in obj:
            yield obj
        return
    if isinstance(obj, Mapping):
        for value in obj.values():
            yield from _extract_pdb_strings(value)
        return
    if isinstance(obj, list):
        for value in obj:
            yield from _extract_pdb_strings(value)


def post_esm_atlas_fold_sequence(
    *,
    sequence: str,
    output_pdb_path: Path,
    timeout_seconds: float = 120.0,
    url: str = ESM_ATLAS_FOLD_URL,
    minimum_ca_fraction: float = 0.80,
) -> InternetStructurePredictionAttempt:
    """POST a single sequence to the ESM Atlas ESMFold endpoint.

    The endpoint returns PDB text on success.  This function is intentionally a
    single bounded HTTP call with no retry loop.
    """

    timeout = max(1.0, float(timeout_seconds))
    minimum_ca_count = int(max(1, len(sequence) * minimum_ca_fraction))
    output_pdb_path = output_pdb_path.resolve()
    request = urllib.request.Request(
        url,
        data=sequence.encode("ascii"),
        headers={"Content-Type": "text/plain"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310 - explicit user-requested endpoint
            body = response.read()
            text = body.decode("utf-8", errors="replace")
            success, atom_count, ca_count = _write_pdb_if_valid(
                pdb_text=text,
                output_pdb_path=output_pdb_path,
                minimum_ca_count=minimum_ca_count,
            )
            return InternetStructurePredictionAttempt(
                predictor_id="esm_atlas_esmfold_single_sequence",
                url=url,
                output_pdb_path=str(output_pdb_path),
                attempted=True,
                success=success,
                status="prediction_available" if success else "response_not_valid_pdb",
                timeout_seconds=timeout,
                http_status=getattr(response, "status", None),
                pdb_atom_line_count=atom_count,
                pdb_ca_line_count=ca_count,
                response_bytes=len(body),
            )
    except (TimeoutError, socket.timeout) as exc:
        return InternetStructurePredictionAttempt(
            predictor_id="esm_atlas_esmfold_single_sequence",
            url=url,
            output_pdb_path=str(output_pdb_path),
            attempted=True,
            success=False,
            status="request_timed_out",
            timeout_seconds=timeout,
            error_tail=_tail(f"{type(exc).__name__}: {exc}"),
        )
    except urllib.error.HTTPError as exc:
        try:
            error_text = exc.read().decode("utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            error_text = str(exc)
        return InternetStructurePredictionAttempt(
            predictor_id="esm_atlas_esmfold_single_sequence",
            url=url,
            output_pdb_path=str(output_pdb_path),
            attempted=True,
            success=False,
            status="http_error",
            timeout_seconds=timeout,
            http_status=exc.code,
            error_tail=_tail(error_text),
        )
    except Exception as exc:  # noqa: BLE001 - propagate through audit artifact, not hang/crash
        return InternetStructurePredictionAttempt(
            predictor_id="esm_atlas_esmfold_single_sequence",
            url=url,
            output_pdb_path=str(output_pdb_path),
            attempted=True,
            success=False,
            status="request_failed",
            timeout_seconds=timeout,
            error_tail=_tail(f"{type(exc).__name__}: {exc}"),
        )


def post_biolm_esmfold(
    *,
    sequence: str,
    output_pdb_path: Path,
    api_key: str | None,
    timeout_seconds: float = 180.0,
    url: str = BIOLM_ESMFOLD_URL,
    minimum_ca_fraction: float = 0.80,
) -> InternetStructurePredictionAttempt:
    """POST to BioLM's ESMFold endpoint when a user supplies an API key."""

    timeout = max(1.0, float(timeout_seconds))
    output_pdb_path = output_pdb_path.resolve()
    if not api_key:
        return InternetStructurePredictionAttempt(
            predictor_id="biolm_esmfold_single_sequence",
            url=url,
            output_pdb_path=str(output_pdb_path),
            attempted=False,
            success=False,
            status="skipped_missing_BIOLM_API_KEY",
            timeout_seconds=timeout,
            sequence_sent_over_network=False,
        )
    payload = json.dumps({"items": [{"sequence": sequence}]}).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=payload,
        headers={"Authorization": f"Token {api_key}", "Content-Type": "application/json"},
        method="POST",
    )
    minimum_ca_count = int(max(1, len(sequence) * minimum_ca_fraction))
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310
            body = response.read()
            text = body.decode("utf-8", errors="replace")
            candidates: list[str] = []
            try:
                candidates.extend(_extract_pdb_strings(json.loads(text)))
            except json.JSONDecodeError:
                candidates.append(text)
            for candidate in candidates:
                success, atom_count, ca_count = _write_pdb_if_valid(
                    pdb_text=candidate,
                    output_pdb_path=output_pdb_path,
                    minimum_ca_count=minimum_ca_count,
                )
                if success:
                    return InternetStructurePredictionAttempt(
                        predictor_id="biolm_esmfold_single_sequence",
                        url=url,
                        output_pdb_path=str(output_pdb_path),
                        attempted=True,
                        success=True,
                        status="prediction_available",
                        timeout_seconds=timeout,
                        http_status=getattr(response, "status", None),
                        pdb_atom_line_count=atom_count,
                        pdb_ca_line_count=ca_count,
                        response_bytes=len(body),
                    )
            atom_count, ca_count = _count_pdb_lines(text)
            return InternetStructurePredictionAttempt(
                predictor_id="biolm_esmfold_single_sequence",
                url=url,
                output_pdb_path=str(output_pdb_path),
                attempted=True,
                success=False,
                status="response_not_valid_pdb",
                timeout_seconds=timeout,
                http_status=getattr(response, "status", None),
                pdb_atom_line_count=atom_count,
                pdb_ca_line_count=ca_count,
                response_bytes=len(body),
            )
    except (TimeoutError, socket.timeout) as exc:
        return InternetStructurePredictionAttempt(
            predictor_id="biolm_esmfold_single_sequence",
            url=url,
            output_pdb_path=str(output_pdb_path),
            attempted=True,
            success=False,
            status="request_timed_out",
            timeout_seconds=timeout,
            error_tail=_tail(f"{type(exc).__name__}: {exc}"),
        )
    except urllib.error.HTTPError as exc:
        try:
            error_text = exc.read().decode("utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            error_text = str(exc)
        return InternetStructurePredictionAttempt(
            predictor_id="biolm_esmfold_single_sequence",
            url=url,
            output_pdb_path=str(output_pdb_path),
            attempted=True,
            success=False,
            status="http_error",
            timeout_seconds=timeout,
            http_status=exc.code,
            error_tail=_tail(error_text),
        )
    except Exception as exc:  # noqa: BLE001
        return InternetStructurePredictionAttempt(
            predictor_id="biolm_esmfold_single_sequence",
            url=url,
            output_pdb_path=str(output_pdb_path),
            attempted=True,
            success=False,
            status="request_failed",
            timeout_seconds=timeout,
            error_tail=_tail(f"{type(exc).__name__}: {exc}"),
        )
