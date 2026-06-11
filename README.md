# Pharmacotopology

Mechanism-language experiments for hard protein-folding regimes.

This repository tests a modest but useful idea: before claiming structure-level success, can a system produce sealed, falsifiable mechanism packets from non-coordinate evidence, then survive post-seal holdout validation?

It is not a universal protein-folding solver. It is not an AlphaFold replacement. It is a small, audit-heavy research workspace for asking whether operator-level biological explanations can be made reproducible, source-separated, and hard to overclaim.

## Current Frontier

Recent live-target packets now cover three different hard classes:

| Version | Target | Class | Status |
| --- | --- | --- | --- |
| V44 | FUS-LC residues 1-214 | IDP / phase separation | sealed packet, review-required |
| V45 | TDP-43 LCD residues 274-414 | IDP / phase separation replication | sealed packet, review-required |
| V46 | CFTR F508del | membrane multidomain folding / proteostasis | committed and synced |
| V47 | RfaH-CTD residues 101-162 | metamorphic alpha/beta fold switching | local review candidate |

Every packet keeps:

```text
folding_problem_solved = false
```

The strongest current claim is narrower:

```text
sealed mechanism-language packets can generate and validate operator-level predictions
across several hard folding regimes without coordinate leakage before prediction
```

## How The Runs Work

Each live-target attack follows the same discipline:

1. Acquire allowed non-coordinate sources.
2. Write a sealed prediction packet and hash it.
3. Keep coordinates, native contacts, AlphaFold-style models, internal runtime artifacts, and holdouts out of prediction.
4. Open independent holdout evidence only after sealing.
5. Score operator buckets, perturbation predictions, leakage controls, and forbidden-claim guards.
6. Emit a certificate and a short report.

The useful output is not a picture of a fold. It is a claim-controlled packet: what mechanism is proposed, what would falsify it, which shortcuts were blocked, and which post-seal evidence supported or contradicted it.

## Reproduce The Latest Packets

From the repository root:

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 scripts/build_v46_cftr_f508del_membrane_multidomain_sources_v0.py
python3 scripts/run_v46_cftr_f508del_membrane_multidomain_attack_v0.py
python3 -m pytest -q tests/test_v46_cftr_f508del_membrane_multidomain_attack.py
```

V47 is intentionally review-gated before commit. Once approved, it should be added with its protocol, focused tests, generated certificate, and report.

## Repository Map

- `scripts/` builds source manifests, sealed packets, validation artifacts, and printed reports.
- `tests/` contains focused regression tests for each packet's scientific and leakage gates.
- `docs/` records the protocol and claim boundary for each version.
- `data/live_unsolved_targets/` stores source manifests, sealed predictions, post-seal holdouts, and validation JSON.
- `first_contact_clean_pharmacotopology_layer_run/` stores run certificates and human-readable reports.

## Scientific Boundary

Useful claims are welcome here. Grand claims are not.

Allowed:

```text
This run produced a sealed, source-separated mechanism packet that passed its
defined holdout and leakage controls.
```

Not allowed:

```text
Universal protein folding is solved.
Coordinates were predicted de novo.
External review is unnecessary.
```

The work is most interesting when it stays humble: each packet should make the next experiment clearer, not make the conclusion louder.
