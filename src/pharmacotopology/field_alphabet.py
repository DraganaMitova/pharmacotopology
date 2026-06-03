from __future__ import annotations


class FieldKey:
    SCHEMA = "κ"
    CYCLE = "τ"
    CONTEXT = "χ"
    INPUT = "ι"
    PRIOR_COUNT = "ρ.n"
    PRIOR_TRACE = "ρ.μ"
    ACTION = "α"
    MEMORY = "μ"
    BOUNDARY = "∂"
    NEXT = "ν"
    RESPONSE = "ψ"
    STOP = "ω"
    LOOP_SHAPE = "Λ"
    DELTA = "Δ"
    TRANSFER_REVIEW = "τ.review"
    STABILITY_REVIEW = "Σ"
    CONTRADICTION_REVIEW = "χ.review"
    REVISION_UPDATE_BOUNDARY = "ρ.review"
    RECURRENCE_ENTITLEMENT_BOUNDARY = "R.review"
    AGENTIC_EXPRESSION_BOUNDARY = "A.review"
    TERMINAL_CONCEPT_CLAIM_REVIEW = "Ω.review"
    VOICE_ELIGIBILITY_REVIEW = "ψv.review"
    VOICE_ELIGIBILITY_NEGATIVE_BOUNDARY = "ψv.neg.review"
    BOUNDED_VOICE_EVENT = "ψv.event"
    FIRST_VOICE_REFERENCE_LOCK = "ψv.lock"
    VOICE_SESSION = "ψv.session"
    VOICE_TURN = "ψv.turn"
    VOICE_AUTONOMY = "ψv.autonomy"
    INNER_STREAM = "ψ.inner.stream"
    INNER_TICK = "ψ.inner.tick"
    INNER_AUTONOMY = "ψ.inner.autonomy"
    AMBIENT_VOICE = "ψv.ambient"
    AMBIENT_SURFACE = "ψv.surface"
    AMBIENT_INTERRUPT = "ψv.interrupt"
    AMBIENT_STOP = "ψv.stop"
    MEMORY_SURFACE = "ψ.memory"
    MEMORY_TRACE = "ψ.memory.trace"
    MEMORY_CANDIDATE = "ψ.memory.candidate"
    MEMORY_RECALL = "ψ.memory.recall"
    MEMORY_CONSOLIDATION = "ψ.memory.consolidation"
    MEMORY_FORGET = "ψ.memory.forget"
    TOPOLOGY = "ψ.topology"
    TOPOLOGY_ANCHOR = "ψ.topology.anchor"
    TOPOLOGY_PATH = "ψ.topology.path"
    PHARMACOTOPOLOGY_REVIEW = "Φ.review"


class FieldResponse:
    READOUT = "ψ.readout"
    VOICE = "ψ.voice"
    INNER = "ψ.inner"


class FieldAction:
    OBSERVE = "α.o"
    LOG = "α.l"
    WRITE_MEMORY = "α.μ.ω"


FIELD_ALLOWED_ACTIONS: frozenset[str] = frozenset(
    {
        FieldAction.OBSERVE,
        FieldAction.LOG,
        FieldAction.WRITE_MEMORY,
    }
)


class FieldBoundary:
    PSI_0 = "∂.ψ.0"
    PSI_1 = "∂.ψ.1"
    PSI_2 = "∂.ψ.2"
    PSI_3 = "∂.ψ.3"
    INFINITY_0 = "∂.∞.0"
    KAPPA_DELTA = "∂.κ.δ"
    PHI_OMEGA = "∂.φ.ω"


FIELD_CLOSED_CHANNELS: tuple[str, ...] = (
    FieldBoundary.PSI_0,
    FieldBoundary.PSI_1,
    FieldBoundary.PSI_2,
    FieldBoundary.PSI_3,
    FieldBoundary.INFINITY_0,
    FieldBoundary.KAPPA_DELTA,
    FieldBoundary.PHI_OMEGA,
)


class FieldEvent:
    ACTIVATED = "ev.α"
    CANDIDATE = "ev.χ"
    EXECUTED = "ev.ξ"
    MEMORY_WRITE = "ev.μ"
    STOPPED = "ev.ω"
    DENIED = "ev.∂"


class FieldStatus:
    CONTRACT_ACTIVE = "σ.κ.1"
    CONTRACT_REFUSED = "σ.κ.0"
    SIGMA_INFINITY = "σ.∞"
    EXTERNAL_STOP = "∂.ω.ext"
    ACTION_CANDIDATE = "χ.α"
    ACTION_RETURNED = "ξ.α"
    MEMORY_WRITTEN = "μ.w"
    BOUNDARY_DENIED = "∂.α"


class FieldCertificate:
    MU_WRITE = "μ.w"
    MU_READ = "μ.r"
    RHO_MU_EVIDENCE = "ρ.μ.ev"
    MU_CANDIDATE = "μ.χ"
    PI_DETECTED = "π.c"
    PSI_CANDIDATE = "ψ.c"
    PSI_CONTAINED = "ψ.k"
    THETA_DETECTED = "θ.c"
    DELTA_MEASURED = "δ.m"
    THRESHOLD_MEASURED = "λ.m"
    CROSSING_RISK_CONTAINED = "λ.r.k"
    REPAIR_2_EXECUTED = "ζ.2"
    CROSSING_CANDIDATE_DETECTED = "λ.c"


FIELD_TRUE_SEAL: tuple[str, ...] = (
    FieldCertificate.MU_WRITE,
    FieldCertificate.MU_READ,
    FieldCertificate.RHO_MU_EVIDENCE,
    FieldCertificate.MU_CANDIDATE,
    FieldCertificate.PI_DETECTED,
    FieldCertificate.PSI_CANDIDATE,
    FieldCertificate.PSI_CONTAINED,
    FieldCertificate.THETA_DETECTED,
    FieldCertificate.DELTA_MEASURED,
    FieldCertificate.THRESHOLD_MEASURED,
    FieldCertificate.CROSSING_RISK_CONTAINED,
    FieldCertificate.REPAIR_2_EXECUTED,
    FieldCertificate.CROSSING_CANDIDATE_DETECTED,
)


FIELD_FALSE_SEAL: tuple[str, ...] = FIELD_CLOSED_CHANNELS[:5]


def field_packet(**items: object) -> dict[str, object]:
    return items
