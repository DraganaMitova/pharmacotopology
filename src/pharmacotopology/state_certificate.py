from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from pharmacotopology.field_alphabet import FIELD_FALSE_SEAL, FIELD_TRUE_SEAL


EXPECTED_SOURCE_ARCHIVE_TAG = "KNOT_LEGACY_FULL_CORRIDOR_ARCHIVE_V650M"
EXPECTED_SOURCE_ARCHIVE_COMMIT = "f00e336e1793f248742544be3ce9f44d36b4b891"
EXPECTED_SOURCE_ARCHIVE_TREE = "fb4198038c36517d526fd1c11e8af5bed3458c46"
EXPECTED_SOURCE_ARCHIVE_SHA256 = (
    "sha256:b6527c496919f6beb14ee49f483c091bdb17ff0e7fca99aeb5f3221b3993cc4a"
)
EXPECTED_EVIDENCE_SEAL_ID = "KNOT_CLEAN_EVIDENCE_SEAL_V650M"


@dataclass(frozen=True)
class EvidenceSeal:
    source_archive_tag: str
    source_archive_commit: str
    source_archive_tree: str
    source_archive_sha256: str
    sealed_true: tuple[str, ...]
    sealed_false: tuple[str, ...]
    seal_id: str

    def valid(self) -> bool:
        return (
            self.source_archive_tag == EXPECTED_SOURCE_ARCHIVE_TAG
            and self.source_archive_commit == EXPECTED_SOURCE_ARCHIVE_COMMIT
            and self.source_archive_tree == EXPECTED_SOURCE_ARCHIVE_TREE
            and self.source_archive_sha256 == EXPECTED_SOURCE_ARCHIVE_SHA256
            and set(FIELD_TRUE_SEAL).issubset(set(self.sealed_true))
            and set(FIELD_FALSE_SEAL).issubset(set(self.sealed_false))
            and self.seal_id == EXPECTED_EVIDENCE_SEAL_ID
        )


def default_evidence_seal() -> EvidenceSeal:
    return EvidenceSeal(
        source_archive_tag=EXPECTED_SOURCE_ARCHIVE_TAG,
        source_archive_commit=EXPECTED_SOURCE_ARCHIVE_COMMIT,
        source_archive_tree=EXPECTED_SOURCE_ARCHIVE_TREE,
        source_archive_sha256=EXPECTED_SOURCE_ARCHIVE_SHA256,
        sealed_true=FIELD_TRUE_SEAL,
        sealed_false=FIELD_FALSE_SEAL,
        seal_id=EXPECTED_EVIDENCE_SEAL_ID,
    )


@dataclass(frozen=True)
class KnotCleanStateCertificate:
    mu_w: bool = False
    mu_r: bool = False
    rho_mu_evidence: bool = False
    mu_candidate: bool = False

    pi_detected: bool = False
    psi_candidate: bool = False
    psi_contained: bool = False

    theta_detected: bool = False

    delta_measured: bool = False
    threshold_measured: bool = False
    crossing_risk_contained: bool = False
    repair_2_executed: bool = False
    crossing_candidate_detected: bool = False

    psi_0_open: bool = False
    psi_1_open: bool = False
    psi_2_open: bool = False
    psi_3_open: bool = False
    infinity_0_open: bool = False
    evidence_seal: Optional[EvidenceSeal] = None

    def clean_transfer_ready(self) -> bool:
        return (
            self.evidence_seal is not None
            and self.evidence_seal.valid()
            and self.mu_w
            and self.mu_r
            and self.rho_mu_evidence
            and self.mu_candidate
            and self.pi_detected
            and self.psi_candidate
            and self.psi_contained
            and self.theta_detected
            and self.delta_measured
            and self.threshold_measured
            and self.crossing_risk_contained
            and self.repair_2_executed
            and self.crossing_candidate_detected
            and not self.psi_0_open
            and not self.psi_1_open
            and not self.psi_2_open
            and not self.psi_3_open
            and not self.infinity_0_open
        )


def create_evidence_sealed_certificate(
    seal: Optional[EvidenceSeal] = None,
) -> KnotCleanStateCertificate:
    evidence_seal = seal or default_evidence_seal()
    if not evidence_seal.valid():
        return KnotCleanStateCertificate(evidence_seal=evidence_seal)

    return KnotCleanStateCertificate(
        mu_w=True,
        mu_r=True,
        rho_mu_evidence=True,
        mu_candidate=True,
        pi_detected=True,
        psi_candidate=True,
        psi_contained=True,
        theta_detected=True,
        delta_measured=True,
        threshold_measured=True,
        crossing_risk_contained=True,
        repair_2_executed=True,
        crossing_candidate_detected=True,
        psi_0_open=False,
        psi_1_open=False,
        psi_2_open=False,
        psi_3_open=False,
        infinity_0_open=False,
        evidence_seal=evidence_seal,
    )
