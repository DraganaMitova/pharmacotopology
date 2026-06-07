"""
Diagnostic test for 1CLL:A collapse bottleneck.

Goal: Identify which of the 36 native long-range contacts are missing
from the selected event set and why (which gate filtered them).

Key thresholds to investigate:
- blocked_future_pressure ≤ 0.16 (hard filter)
- future_preservation_score ≥ 0.62 (evolutionary signal)
- direct_support_score ≥ 0.22 (evolutionary backing)
- contact_cluster_gain ≥ 0.30-0.46 (clustering requirement)
"""

from pathlib import Path
import sys
import json
import pandas as pd
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def test_1cll_a_native_contact_recovery():
    """Dump all 1CLL:A data and identify missing native contacts."""
    
    # Load 1CLL:A structure
    pdb_id = "1CLL"
    chain_id = "A"
    target_id = f"coord_008_pdb_{pdb_id}_{chain_id.upper()}_calmodulin"
    
    print(f"\n{'='*80}")
    print(f"DIAGNOSTIC: 1CLL:A COLLAPSE BOTTLENECK")
    print(f"{'='*80}")
    print(f"Target: {target_id}")
    print(f"PDB: {pdb_id}:{chain_id}")
    
    # Load external coupling data
    manifest_path = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_coupling_target_manifest_v0.json"
    if manifest_path.exists():
        with open(manifest_path) as f:
            manifest = json.load(f)
        
        # Find 1CLL:A entry
        for entry in manifest:
            if entry.get("row_id") == target_id:  # noqa: F841
                print(f"\n{'='*80}")
                print(f"COUPLING DATA FOR 1CLL:A")
                print(f"{'='*80}")
                for key, val in entry.items():
                    print(f"{key:40s}: {val}")
                break
    
    # Load frontier data (event selection history)
    frontier_path = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_coupling_trace_loop_frontier.csv"
    if frontier_path.exists():
        df = pd.read_csv(frontier_path)
        cll_frontier_rows = df[df["row_id"] == target_id]
        print(f"\n{'='*80}")
        print(f"FRONTIER EVENTS FOR 1CLL:A: {len(cll_frontier_rows)} rows")
        print(f"{'='*80}")
        
        if len(cll_frontier_rows) > 0:
            print("\nKey columns in frontier data:")
            print(cll_frontier_rows.columns.tolist())
            
            print(f"\nFirst few 1CLL:A events:")
            display_cols = [
                "source_selector", "row_id", "segment_a_start", "segment_a_end",
                "segment_b_start", "segment_b_end", "coupling_nucleus_score"
            ]
            existing_cols = [c for c in display_cols if c in cll_frontier_rows.columns]
            print(cll_frontier_rows[existing_cols].head(10).to_string())
            
            print(f"\n✓ Total frontier events for 1CLL:A: {len(cll_frontier_rows)}")
    
    # Load selected events
    selected_path = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_coupling_trace_loop_selected_events.csv"
    if selected_path.exists():
        df_selected = pd.read_csv(selected_path)
        cll_selected = df_selected[df_selected["row_id"] == target_id]
        print(f"\n{'='*80}")
        print(f"SELECTED EVENTS FOR 1CLL:A: {len(cll_selected)} events")
        print(f"{'='*80}")
        
        if len(cll_selected) > 0:
            print("\nSelected events (first 20):")
            display_cols = [c for c in ["row_id", "segment_a_start", "segment_a_end", "segment_b_start", "segment_b_end", "source_selector", "coupling_nucleus_score"] 
                           if c in cll_selected.columns]
            print(cll_selected[display_cols].head(20).to_string())
            
            # Extract contact pairs from selected events
            selected_pairs = set()
            if "segment_a_start" in cll_selected.columns:
                for _, row in cll_selected.iterrows():
                    # Use the segment endpoints as the contact representation
                    a_start = int(row["segment_a_start"])
                    a_end = int(row["segment_a_end"])
                    b_start = int(row["segment_b_start"])
                    b_end = int(row["segment_b_end"])
                    selected_pairs.add(((a_start, a_end), (b_start, b_end)))
            
            print(f"\nSelected contact pairs: {len(selected_pairs)}")
            for pair in sorted(selected_pairs):
                print(f"  {pair[0]} -- {pair[1]}")
    
    # Analyze scoring/gating data if available
    frontier_path = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_coupling_trace_loop_frontier.csv"
    if frontier_path.exists():
        df = pd.read_csv(frontier_path)
        cll_analysis_rows = df[df["row_id"] == target_id]
        
        # Look for threshold columns
        numeric_cols = cll_analysis_rows.select_dtypes(include=['float64', 'float32']).columns.tolist()
        print(f"\n{'='*80}")
        print(f"SCORING COLUMNS AVAILABLE")
        print(f"{'='*80}")
        for col in sorted(numeric_cols)[:20]:
            print(f"  {col}")
        if len(numeric_cols) > 20:
            print(f"  ... and {len(numeric_cols) - 20} more")
        
        # Check for blocked_future_pressure specifically
        potential_cols = [c for c in numeric_cols if 'future' in c.lower() or 'blocked' in c.lower()]
        if potential_cols:
            print(f"\n{'='*80}")
            print(f"FUTURE/BLOCKED PRESSURE COLUMNS")
            print(f"{'='*80}")
            for col in potential_cols:
                print(f"\n{col}:")
                print(f"  min: {cll_analysis_rows[col].min():.4f}")
                print(f"  max: {cll_analysis_rows[col].max():.4f}")
                print(f"  mean: {cll_analysis_rows[col].mean():.4f}")
                print(f"  median: {cll_analysis_rows[col].median():.4f}")
                
                # Count how many fall outside typical gates
                if "blocked" in col.lower():
                    outside_016 = (cll_analysis_rows[col] > 0.16).sum()
                    if outside_016 > 0:
                        print(f"  Events with {col} > 0.16: {outside_016}")
    
    print(f"\n{'='*80}")
    print(f"NEXT STEPS")
    print(f"{'='*80}")
    print("""
1. Load 1CLL:A structure from data/rcsb_pdb/ to get exact native contacts
2. Dump all unselected events for 1CLL:A and their gate failure reasons
3. For each missing native contact, find which unselected event covers it
   and identify the exact gate that rejected it
4. Propose threshold adjustments or new scorer logic
5. Test recovery with minimal changes
    """)


if __name__ == "__main__":
    test_1cll_a_native_contact_recovery()
