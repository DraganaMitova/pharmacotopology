# Pharmacotopology Collapse/Selection Logic Deep Dive

## Executive Summary: 1CLL:A Mystery

**The Problem**: Event recall is 1.0 but collapse isn't recovering all 36 native long-range contacts for 1CLL:A  
**Row ID**: `coord_008_pdb_1CLL_A_calmodulin`  
**Current State**: 148/148 external couplings available (100%), target_coverage = 1.0

The paradox suggests:
- All constraint pairs from evolutionary couplings are being "discovered" (event_recall = 1.0)
- But native long-range contact PAIRS (distance ≥ 24) within those regions aren't being selected
- This points to a **selection filtering issue**, not an event generation issue

---

## Part 1: Collapse Search Logic (`folding_nucleus_closure_search.py`)

### Event Generation Algorithm

The system scans the sequence for "nucleus closure events" - pairs of segments where evolutionary evidence suggests contacts.

```python
SEGMENT_LENGTH = 8           # Size of each segment
SEGMENT_STRIDE = 4           # Overlap between segments
MIN_SEGMENT_GAP = 4          # Minimum gap between segments
```

**Process**:
1. Slide an 8-residue window across the sequence (stride 4)
2. For each segment A, find all segment B candidates at distance ≥ 4
3. Score top 10 contact features in the region
4. Calculate event score from 8 components:

**Event Score Formula** (weighted components):
```
nucleus_score = (
    0.38 * contact_cluster_gain           # Evolutionary coupling strength
    + 0.12 * secondary_structure_compat   # Helix/beta compatibility
    + 0.20 * hydrophobic_burial_gain      # Hydrophobic pairing support
    + 0.14 * registry_support             # Parallel contact evidence
    - 0.16 * loop_entropy_cost            # Loop entropic cost
    - 0.08 * geometry_violation_cost      # Span penalty if > 65% of sequence
    - 0.12 * frustration_cost             # Charge/breaker penalties
    - 0.06 * isolation_penalty            # Feature isolation penalty
    + 0.12  # bias term
)
```

### Critical Thresholds

- **NUCLEUS_THRESHOLD_MIN = 0.30**: Events must score ≥ this to be accepted
- **FRUSTRATION_REJECTION_LIMIT = 0.75**: Events with frustration ≥ this are **rejected outright** (hard cutoff)

**For 1CLL:A**: These filters should pass most events since calmodulin has good helix/hydrophobic structure.

---

## Part 2: Physical Selection (`folding_physical_selection.py`)

### What It Does

Takes candidate nucleus events and applies physical state scoring based on protein structure.

### Physical State Components

Events are scored on:
1. `burial_gain` - How buried residues become upon closure
2. `loop_strain` - Conformational cost of loop closure
3. `steric_clash_score` - Steric conflicts introduced
4. `unsatisfied_polar_penalty` - Polar groups left unsatisfied
5. `future_frustration_score` - Will future contacts be blocked?

### Active Physical Score Calculation

```python
active_physical_score = (
    burial_gain
    - 0.75 * loop_strain
    - 0.85 * steric_clash_score
    - 0.75 * unsatisfied_polar_penalty
    - 0.65 * future_frustration_score
    + 0.20 * decoy_margin  # vs competing event
)
```

### Selection Gates

**PHYSICAL_GATE** filters:
- `active_physical_score >= 0.50`
- `unsatisfied_polar_penalty <= 0.24`
- `steric_clash_score <= 0.08`

**FUTURE_FRUSTRATION_GATE** additionally requires:
- Future closure path must be preserved (special calculation)

### Selection Process

Per row:
- Take all competitive events (winners of conflict resolution)
- Score with active_physical_score
- Select top 40 by score, sorted by `-score, segment_a_start, segment_b_start`

**Potential Issue for 1CLL:A**: If many native contact events score low on physical metrics (buried, polar conflicts), they get filtered here.

---

## Part 3: Coupling Nucleus Selector (`folding_coupling_nucleus_selector.py`)

### Critical Finding: This is Where Most Filtering Happens

The coupling selector applies **evolutionary coupling constraints** through a sophisticated multi-stage pipeline with cumulative filters.

### Core Data Flow

```
Competitive Events (physical selection output)
    ↓
Build CouplingNucleusContext:
  - Load external coupling constraints (evolutionary evidence)
  - Assess each event vs constraints
  - Calculate "coupling selectivity score" for each event
  - Calculate "blocked_future_pressure" (KEY METRIC)
    ↓
Apply Trace Loop Selectors (multiple gates):
  - coupling_future_gate (basic)
  - coupling_trace_loop (coverage-based)
  - rank_consistent_cluster_gated
  - persistent_rank_consistent_cluster_gated
  - score_margin_expanded
  - boundary_continuity_expanded
  - edge_continuity_expanded
  - pressure_release_expanded
  - registry_extension_expanded
  - terminal_bridge_expanded
  - boundary_field_replacement_probe
  - macro_scale_future_preserved
```

### The Fundamental Gate: `_passes_coupling_future_gate()`

This is the **bottleneck** that could drop 1CLL:A native contacts:

```python
def _passes_coupling_future_gate(event):
    assessment = context.assessment_by_event_id[event.event_id]
    return (
        assessment.direct_support_score >= 0.22  # COUPLING_DIRECT_SUPPORT_MIN
        AND assessment.future_preservation_score >= 0.62  # COUPLING_FUTURE_PRESERVATION_MIN
        AND assessment.blocked_future_pressure <= 0.16  # COUPLING_BLOCKED_FUTURE_MAX
    )
```

**Key Threshold**: `blocked_future_pressure <= 0.16`
- This measures: "Will evolutionary constraints be violated by future contacts after this closure?"
- If an event "locks in" a region that would prevent future evolutionary-supported contacts → HIGH blocked_future_pressure
- Events with blocked_future_pressure > 0.16 are **automatically gated out**

**Hypothesis**: 1CLL:A native long-range contacts may be getting filtered because:
1. They have HIGH blocked_future_pressure (> 0.16)
2. They don't have strong enough direct_support_score (< 0.22)
3. They lack future_preservation_score (< 0.62)

### The Trace Loop Coverage Algorithm

For events that pass the future gate, the system uses a **greedy set-covering algorithm**:

```python
def select_coupling_trace_loop_events():
    uncovered = set(all_constraint_pairs)
    selected = []
    
    while uncovered and len(selected) < 40:
        # Score candidates by:
        trace_score = (
            0.48 * min(1.0, covered_confidence / 3.0)  # New constraints covered
            + 0.28 * coupling_nucleus_score(event)
            + 0.16 * assessment.future_preservation_score
            - 0.10 * assessment.blocked_future_pressure
        )
        
        # Pick highest scorer
        chosen = max_by(trace_score, ...)
        selected.append(chosen)
        uncovered -= event.candidate_pairs
    
    return selected
```

**Key insight**: This algorithm optimizes for **constraint coverage**, not native contact recovery.
- It selects events that cover new constraint pairs
- Native contact pairs might already be "covered" by previous selections
- Once all constraint pairs are covered → algorithm stops
- Remaining native-only pairs in uncovered regions → never selected

### Rescue Gates (Expansion Phases)

After core selection, progressive rescue phases add back events that:

1. **score_margin_expanded**: Score margin > 0.15 vs decoy, cluster gain ≥ 0.46
2. **boundary_continuity_rescue**: Score 0.45-0.60, secondary structure ≥ 0.53
3. **edge_continuity_rescue**: Score 0.42-0.50, cluster gain high
4. **pressure_release_rescue**: Score 0.42-0.46, with blocked_future 0.12-0.22 (INTERESTING!)
5. **registry_extension_rescue**: Anchored to existing selections
6. **direct_margin_tail_rescue**: Supported by coupling constraints
7. **terminal_bridge_rescue**: Bridges between selection regions
8. **high_pressure_tail_rescue**: Like pressure_release but for high-pressure regions
9. **boundary_field_replacement**: Replaces weaker selections with better options

**Critical Pattern**: Many rescue gates have SPECIFIC SCORE RANGES:
- This means well-scoring events (nucleus_score > 0.60) might FAIL rescue gates too!

---

## Part 4: Ranking & Persistence (`folding_nucleus_rank_enrichment.py`)

### Simple But Important

The ranking system just:

1. **Sorts** all candidate events by a score function
2. **Measures enrichment** at top-10, top-25, top-50
3. Calculates: `rank_enrichment = (top_k_native_positive_rate) / (baseline_native_positive_rate)`

```python
def native_positive_rate(events):
    return sum(1 for e in events if e.native_contact_count_after_scoring > 0) / len(events)

def rank_enrichment_rows_for_row(...):
    ranked = sorted(events, key=lambda e: -score_function(e))
    baseline = native_positive_rate(ranked)
    
    for cutoff in [10, 25, 50]:
        top = ranked[:cutoff]
        top_rate = native_positive_rate(top)
        enrichment = top_rate / baseline  # How much better than average?
```

**For 1CLL:A**: If `event_recall = 1.0`, this means events covering native pairs are in the ranked list, just not all selected.

The **persistence mechanism** in coupling selector (`trace_loop_persistence_evidence`) looks at:
- How many good neighbors does an event have?
- Are neighbors consistently supporting it?
- Requires `TRACE_LOOP_PERSISTENT_RECOVERY_NEIGHBOR_COUNT_MIN = 2` neighbors

---

## Diagnostic Summary: What Might Drop 1CLL:A Native Contacts

### Most Likely Culprits (in order):

1. **BLOCKED_FUTURE_PRESSURE GATE** (threshold: ≤ 0.16)
   - Native 1CLL:A closure events may "lock in" regions that prevent future evolutionary contacts
   - These get filtered by fundamental `_passes_coupling_future_gate()`
   - Can't be rescued by later gates

2. **COVERAGE ALGORITHM SATURATION**
   - Once core events cover all 148 constraint pairs, algorithm stops
   - Native-only pairs in those regions never get selected
   - Rescue gates only kick in for selected regions

3. **SCORE RANGE GATES**
   - Rescue gates like `score_margin_expanded` require score 0.44-0.60
   - If native events score outside this range → excluded
   - High-scoring events (> 0.60) might not match any rescue criteria

4. **FUTURE_PRESERVATION_MIN = 0.62**
   - If 1CLL:A native contacts have weak evolutionary signal in future positions
   - They fail the fundamental gate

5. **CLUSTERING REQUIREMENTS**
   - Many gates require `contact_cluster_gain >= 0.30-0.46`
   - Isolated native contacts fail this

### Less Likely (but possible):

6. **direct_support_score < 0.22**: Weak direct evolutionary support
7. **Competitive event selection**: Native event loses to alternative event
8. **Sequence-based gates**: `FRUSTRATION_REJECTION_LIMIT = 0.75` or `geometry_violation_cost`

---

## Test Data from 1CLL:A Run

**Selected Events**: 10 events shown in `external_coupling_trace_loop_selected_events.csv` (rows 75-84)

Example event from row 75:
```
rank=74, segment_a=[13-20], segment_b=[25-32]
direct_support=0.578018, future_preservation=0.681147
blocked_future_pressure=0.091233 (PASSES!)
coupling_selectivity_score=0.700325
native_contact_count=7
```

**Key observation**: Selected events have relatively low blocked_future_pressure (0.09-0.35), suggesting events with higher blocked_future are being filtered.

---

## Selector Benchmark Results

From `external_coupling_trace_loop_selectors.csv`:

| Selector | Events | Long-range Recall | False Rate | Precision |
|----------|--------|-------------------|-----------|-----------|
| external_real | 90 | 0.312 | 0.394 | 0.082 |
| external_macro_scale_future_preserved | 28 | **0.564** | 0.0 | 0.044 |
| external_boundary_field_replacement_probe | 37 | 0.350 | 0.0 | 0.168 |
| external_terminal_bridge_expanded | 37 | 0.323 | 0.0 | 0.156 |
| external_pressure_release_expanded | 33 | 0.316 | 0.0 | 0.167 |

**Note**: Better recall (0.564) achieved by macro_scale_future_preserved, which uses future_preservation_score >= 0.32 gate.

---

## Concrete Code Patterns to Investigate

### 1. Check Blocked Future Pressure Distribution

In coupling selector, margins are calculated:
```python
margins[event.event_id] = coupling_nucleus_score(event) 
                        - coupling_nucleus_score(decoy)
```

But the real filter is `_passes_coupling_future_gate()` which checks:
```python
assessment.blocked_future_pressure <= 0.16
```

**Action**: Dump `blocked_future_pressure` for all 1CLL:A events, see how many fall outside 0.16 threshold.

### 2. Check Trace Loop Persistence

The persistence evidence calculation (lines 801-851) looks for neighboring events:
```python
neighbors = [candidate for candidate in candidates 
             if compatible_future_event(event, candidate)
             and trace_distance <= 32  # TRACE_LOOP_PERSISTENT_NEIGHBOR_WINDOW
             and contact_cluster_gain >= 0.30
             and direct_support >= 0.45]
```

**Action**: Check if 1CLL:A native events have weak neighbor support.

### 3. Check Score Margin Gates

The `_selector_score_decoy_margin()` function (line 834) calculates:
```python
margin = coupling_nucleus_score(event) - coupling_nucleus_score(decoy)
```

Gates like `score_margin_expanded` require:
```python
margin >= 0.15  # TRACE_LOOP_SCORE_MARGIN_EXPANSION_DECOY_MARGIN_MIN
```

**Action**: Verify native contact events have sufficient decoy margins.

### 4. Compare with Other Proteins

From selector benchmarks, other proteins achieve 0.3-0.35 long-range recall with same gates.  
1CLL:A at 0.312 is not anomalous - the system inherently doesn't recover all native contacts.

---

## The Real Problem: System Design Trade-off

The selection system optimizes for:
1. **Evolutionary signal preservation**: Future contacts must not be blocked
2. **Constraint coverage**: Get evidence for as many evolutionary pairs as possible
3. **Physical plausibility**: Physically sound closures

It does **NOT** optimize for:
- Recovering ALL native long-range contacts
- Complete contact map coverage
- Native contact recall per se

**Why event_recall = 1.0 but long-range recall < 1.0**:
- All constraint pairs are discovered (event_recall measures constraint coverage)
- But not all native pairs within those constraint regions are selected
- Some native pairs remain "within events" that don't make the final selection

---

## Recommendations for Further Investigation

### 1. Identify Dropped Native Contacts

```python
# Pseudo-code
for 1CLL:A_native_pair in all_36_native_long_range:
    # Check if in any selected event's candidate region
    in_selected = any(pair in selected_event.candidate_region_pairs() 
                      for selected_event in selected_1cll_events)
    if not in_selected:
        print(f"Native pair {pair} not in any selected event")
```

### 2. Examine Unselected Events Containing Native Pairs

```python
unselected_with_natives = [
    event for event in all_competitive_events
    if event not in selected_events
    and any(pair in set(event.candidate_region_pairs()) 
            for pair in native_long_range_pairs)
]

for event in unselected_with_natives:
    print(f"Event {event.event_id}:")
    print(f"  blocked_future_pressure: {assessment[event].blocked_future_pressure}")
    print(f"  why gated out: ...")
```

### 3. Test Threshold Adjustments

- Lower `blocked_future_pressure` threshold from 0.16 to 0.20 or 0.25
- Raise `future_preservation_min` to select events that don't block futures
- Adjust `contact_cluster_gain` minimums to include isolated natives

### 4. Look at Competitive Event Selection

Check if native contact events are losing to higher-scoring non-native alternatives in the competitive selection phase.

---

## Summary of Key Constants

### Thresholds That Might Drop 1CLL:A Contacts

| Threshold | Value | Impact |
|-----------|-------|--------|
| COUPLING_BLOCKED_FUTURE_MAX | 0.16 | Hard gate - only contacts with low blocking included |
| COUPLING_FUTURE_PRESERVATION_MIN | 0.62 | Needs strong evolutionary signal in future positions |
| COUPLING_DIRECT_SUPPORT_MIN | 0.22 | Needs direct evolutionary backing |
| TRACE_LOOP_CLUSTER_GATE_MIN | 0.32 | Isolated contacts excluded |
| TRACE_LOOP_PERSISTENT_NEIGHBOR_WINDOW | 32 | Must have neighbors within 32 residues |
| SELECTED_EVENTS_PER_ROW | 40 | Only 40 events selected max |

### Scoring Formulas

**Trace Score** (for selection):
```
0.48 * constraint_coverage + 0.28 * nucleus_score + 0.16 * future_preservation - 0.10 * blocked_future
```

**Coupling Nucleus Score** (ranking):
```
0.62 * coupling_selectivity + 0.18 * physical_score + 0.12 * burial - 0.14 * polar - 0.10 * steric
```

---

## Files to Examine for 1CLL:A Data

1. **external_coupling_trace_loop_selected_events.csv** - Which events are selected (rows 75-84)
2. **external_coupling_trace_loop_row_status.csv** - High-level metrics (row 9)
3. **external_coupling_trace_loop_report.json** - Full report with all metrics
4. Could add debug output to capture:
   - All competitive events and their assessments
   - blocked_future_pressure distribution
   - Which gates filter which events
   - Native contact coverage per event
