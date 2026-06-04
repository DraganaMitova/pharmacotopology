# KNOT Particle-Like Closure Toy Universe Boundary

This branch adds a bounded simulation question:

Can simple mark/field/topology rules produce stable particle-like forms?

It does not claim that real particles have been discovered, modeled, or explained.
The output is a toy-universe emergence trace only.

## Primitive

Each lattice node carries:

- `state`: activation or mark intensity.
- `phase`: cyclic orientation used for resonance and winding.
- `charge`: signed polarity, limited to `-1`, `0`, or `1`.
- `memory`: local attractor residue from prior coherent closure.
- `links`: implicit lattice links from the configured neighbor offsets.

## Update Rule

The simulation applies a small deterministic rule set:

- Difference spreads through neighbor activation and opposite-charge pressure.
- Similar phases resonate weakly.
- Unstable activation decays.
- Closed eight-cell loops stabilize when occupancy, charge coherence, phase
  winding, and memory are all high enough.
- Attention is modeled as selection pressure that projects a coherent closure
  into a neighboring loop template.
- Optional decay emission lets a sufficiently stable closure shed paired
  daughter templates. These daughters are still only counted if the ordinary
  closure detector tracks them over time.
- Stress emission replaces the artificial emission clock with a tracked
  internal pressure accumulator. A closure emits when its own persistence,
  mass pressure, memory, winding, and crowding cross a threshold.
- Optional resource budgeting makes branching non-free. When too many loose
  closures compete for field support, reinforcement weakens and over-branching
  can collapse. Opposite-charge closures inside a bound radius receive a small
  resource discount as a first toy analogue of composite stability.
- Optional confinement marks daughter closures as fragment-like. Fragment-like
  tracks are penalized when isolated after a grace period, but keep normal field
  support when they remain near opposite-charge closure partners.

Speed is not treated as the source of emergence. The rule set is the hypothesis;
faster runs only reveal its consequences sooner.

## First Milestone

The first goal is:

Can a seed become a stable moving loop?

The runner reports this through:

- `stable_loop_count`
- `moving_loop_count`
- `max_track_persistence`
- `first_goal_met`

## One-Seed Decay Question

A single seed is meaningful, but only for a narrower question:

Can one coherent closure persist, decay, shed residue, and create daughter
closures under the same field rules?

Without a decay or emission channel, one seed can only persist, move, or die.
It cannot produce an ecology by itself. With decay emission enabled, the useful
measurement is whether daughters:

- persist as stable closures
- form opposite-charge or opposite-winding families
- annihilate, decay, or remain confined near the parent
- conserve charge-like balance or drift into runaway growth

Runaway branching is a failed or over-fertile rule regime, not evidence of real
particle creation.

Periodic emission is a comparison control, not the deep target. The deeper
question uses `decay_emission_mode="stress"`:

Can the field generate birth times from internal closure pressure instead of a
fixed external clock?

In stress mode, prime-looking birth times are not meaningful by themselves. The
question is whether intervals become structured by internal dynamics, crowding,
resource collapse, or bound-state survival.

## Quark-Like Toy Criterion

The simulation may describe a result as a quark-like toy analogue only when all
of the following are true inside the toy field:

- fragment-like daughter tracks are born from an initial closure
- stable unbound fragment-like tracks are absent or strongly suppressed
- bound fragment-like tracks persist near opposite-charge partners
- the run reports `confined_fragment_track_count > 0`
- the run still keeps `known_particle_claim_created = false`

This does not identify real quarks. It means the toy field produced confined
substructure: fragments whose survival depends on local composite binding.

## Bounded Claims

Allowed:

- particle-like stable localized loops
- one-seed daughter-loop emergence inside the toy field
- confined fragment-like substructure inside the toy field
- charge-like polarity inside the toy field
- spin-like phase winding inside the toy field
- collision and annihilation counters inside the toy field
- repeatable loop families inside the toy field

Not allowed:

- known-particle discovery claims
- electron, proton, neutron, photon, or Standard Model identification
- physical mass, physical charge, or physical spin claims
- claims that acceleration alone creates particles
