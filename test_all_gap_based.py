import sys
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

def test_all_proteins_gap_based():
    df = pd.read_csv('first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_frontier.csv')
    
    results = []
    
    for row_id in df['row_id'].unique():
        prot = df[df['row_id'] == row_id]
        
        scores = sorted([row['future_preservation_score'] for _, row in prot.iterrows()])
        
        if len(scores) < 2:
            gap_threshold = 0.0
        else:
            gaps = [scores[i+1] - scores[i] for i in range(len(scores)-1)]
            max_gap_index = max(range(len(gaps)), key=lambda i: gaps[i])
            gap_threshold = scores[max_gap_index + 1] # boundary is the value just above the gap
            
        pass_count = sum(1 for s in scores if s >= gap_threshold)
        
        print(f"{row_id[:25]:<25} | Gap Threshold: {gap_threshold:.4f} | Pass: {pass_count}/{len(prot)}")

if __name__ == "__main__":
    test_all_proteins_gap_based()
