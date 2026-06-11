# V54 Simulated Mutation And Condition Response

V54 adds perturbation simulation. The simulator predicts direction rather than exact atomic coordinates.

## Required Comparisons

- Wild type or reference state.
- Known damaging mutation.
- Known rescue mutation or condition.
- Neutral or control perturbation.
- Wrong-region perturbation.

## Directional Examples

- FUS/TDP aromatic or charge perturbations shift phase tendency.
- CFTR F508del weakens NBD1/interface/proteostasis routing.
- CFTR corrector-like context improves the routing metric.
- RfaH partner/release perturbations shift alpha/beta basin occupancy.
- ORF6 C-terminal perturbation weakens host-interface hijacking.

Pass condition: correct perturbations move the relevant state variable in the expected direction, while wrong perturbations do not look equally correct.
