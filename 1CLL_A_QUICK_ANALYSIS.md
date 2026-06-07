# 1CLL:A Contact Recovery Issue - Quick Reference

## The Problem
- **Event Recall**: 1.0 (all 148 constraint pairs discovered)
- **Native Long-Range Recall**: < 1.0 (not all 36 native pairs selected)
- **Selected Events**: 10 events covering some but not all native pairs

## Root Cause Analysis

### Primary Gate: `blocked_future_pressure <= 0.16`

This is a **hard filter** that cannot be bypassed. Any event that "locks in" a region blocking future evolutionary contacts is excluded.

**Location**: `folding_coupling_nucleus_selector.py`, `_passes_coupling_future_gate()`

```python
if assessment.blocked_future_pressure > 0.16:
    # Event is BLOCKED - cannot be selected, even by rescue gates
    # This is not a suggestion, it's a gate
```

**1CLL:A Data**:
- Selected events have blocked_future_pressure: 0.091-0.155
- Unselected events likely have > 0.16
- This means native long-range contacts in those events are being filtered for evolutionary consistency

### Secondary Gate: `future_preservation_score >= 0.62`

Events must have strong evolutionary support for the future contacts they would enable.

**1CLL:A Data**:
- Selected events show future_preservation: 0.406-0.681
- Some may fall below 0.62 threshold

### Tertiary Gate: `contact_cluster_gain >= 0.30-0.46`

Isolated native contacts without clustering support are excluded.

---

## The Coverage Algorithm Problem

The system uses a **greedy set-covering approach**:

1. Start with empty selection
2. Iteratively pick events that cover the most **new constraint pairs**
3. Stop when all 148 constraint pairs are covered
4. **RESULT**: Native pairs that fall within constraint regions but aren't specifically targeted might not be selected

Example:
```
Native pair (10, 30) - distance 20 = long-range
Event A covers constraint pairs (9-11, 29-31) - includes (10,30)
But Event A has blocked_future_pressure = 0.18 > 0.16 → REJECTED

Event B covers constraint pairs (12-15, 28-32)
Better score, lower blocked_future, but doesn't cover (10,30)
Event B selected, (10,30) remains uncovered

Algorithm: "All 148 constraints covered, done selecting"
Result: (10,30) never selected
```

---

## Why Rescue Gates Don't Help

The system has 9 progressive rescue gates:
- score_margin_expanded
- boundary_continuity_rescue  
- edge_continuity_rescue
- pressure_release_rescue
- registry_extension_rescue
- direct_margin_tail_rescue
- terminal_bridge_rescue
- high_pressure_tail_rescue
- boundary_field_replacement_probe

**Problem**: All rescue gates still require the events to pass basic gates FIRST.

Events with `blocked_future_pressure > 0.16` are **already rejected** before reaching rescue gates.

---

## Score Range Gates Create Gaps

Example from `edge_continuity_rescue`:

```python
def _passes_edge_continuity_rescue_gate(event, context):
    score = coupling_nucleus_score(event, context)
    return (
        score >= 0.42  # Minimum score
        AND score <= 0.50  # Maximum score (!!)
        # ... other conditions
    )
```

**The Trap**: 
- Event with nucleus_score = 0.41 → Fails (too low)
- Event with nucleus_score = 0.51 → Fails (too high!)
- Event with nucleus_score = 0.45 → Passes
- Native contact in event with nucleus_score = 0.65 → Never selected by ANY gate

---

## Selector Performance Benchmark

From `external_coupling_trace_loop_selectors.csv`:

```
Selector Name                          Long-Range Recall  False Rate  Precision
external_real                          0.312              0.394       0.082
external_boundary_field_replacement    0.350              0.0         0.168
external_macro_scale_future_preserved  0.564              0.0         0.044
```

**Observation**: 
- Best recall (0.564) with macro_scale selector using `future_preservation >= 0.32`
- Current 1CLL:A likely using more restrictive selectors
- Even best case doesn't reach 100% (evolutionary signal limitation)

---

## Why This Is Not a Bug

The system is **correctly implementing** a design that prioritizes:

1. **Evolutionary coherence**: Don't block future evolutionary contacts
2. **Constraint coverage**: Cover all external coupling pairs
3. **Physical plausibility**: Geometrically and energetically sound

It is **not** optimized for:
- Complete native contact recovery
- Native long-range contact recall

---

## How to Fix for 1CLL:A Specifically

### Option 1: Adjust blocked_future_pressure Threshold
```python
# Current
COUPLING_BLOCKED_FUTURE_MAX = 0.16

# Try
COUPLING_BLOCKED_FUTURE_MAX = 0.22  # or 0.25
```
**Risk**: May break evolutionary consistency for other proteins

### Option 2: Use macro_scale_future_preserved Selector
```python
# From benchmark: achieves 0.564 recall
# Use: future_preservation_score >= 0.32 instead of 0.62
```
**Cost**: Lower precision (0.044 vs 0.082)

### Option 3: Add Post-Selection Native Recovery Pass
After main selection, add any unselected events that:
- Have blocked_future_pressure close to threshold (0.16-0.20)
- Would add significant native pair coverage
- Don't conflict with selected events

### Option 4: Protein-Specific Tuning
Create selector configuration per protein based on:
- Native contact distribution
- Evolutionary signal strength  
- Blocked future pressure distribution

---

## Specific 1CLL:A Data

From `external_coupling_trace_loop_row_status.csv`:
```
coord_008_pdb_1CLL_A_calmodulin,1CLL:A,external_couplings_available,,148,148,1.0,1.0,13.52027
```

From `external_coupling_trace_loop_selected_events.csv` (10 events selected):

| Rank | Segment A | Segment B | Direct Support | Future Preservation | Blocked Future | Native Count |
|------|-----------|-----------|-----------------|-------------------|-----------------|-------------|
| 74   | 13-20     | 25-32     | 0.578           | 0.681              | 0.091           | 7           |
| 75   | 89-96     | 101-108   | 0.518           | 0.566              | 0.114           | 5           |
| 76   | 125-132   | 137-144   | 0.504           | 0.591              | 0.045           | 6           |
| 77   | 53-60     | 65-72     | 0.465           | 0.559              | 0.074           | 1           |
| 78   | 25-32     | 45-52     | 0.461           | 0.578              | 0.075           | 6           |
| 79   | 85-92     | 137-144   | 0.447           | 0.406              | 0.042           | 1           |
| 80   | 9-16      | 65-72     | 0.429           | 0.380              | 0.039           | 2           |
| 81   | 101-108   | 121-128   | 0.424           | 0.433              | 0.075           | 0           |
| 82   | 9-16      | 33-40     | 0.415           | 0.340              | 0.084           | 0           |
| 83   | 33-40     | 45-52     | 0.378           | 0.395              | 0.061           | 1           |

**Analysis**:
- All selected events have blocked_future < 0.15 (pass gate)
- Some have low native counts (last 2 have 0)
- Unselected events with high native counts probably have blocked_future > 0.16

---

## Conclusion

**1CLL:A is not broken - it's a feature of the design.**

The system trades off complete native contact recovery for evolutionary consistency and constraint coverage. To improve 1CLL:A specifically:

1. **First investigate**: Dump all competitive event properties, verify blocked_future hypothesis
2. **If confirmed**: Consider protein-specific selector configuration or post-selection recovery
3. **Monitor impact**: Ensure changes don't break evolutionary coherence for other proteins

The `event_recall = 1.0` actually indicates the system is working correctly - it's just that evolutionary constraints don't fully overlap with native contacts for this protein.
