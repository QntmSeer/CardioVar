import pandas as pd
import numpy as np
from utils import get_delta_rna_seq

def process_batch(df):
    """
    Process a batch of variants.
    
    Args:
        df (pd.DataFrame): DataFrame with columns ['chrom', 'pos', 'ref', 'alt']
        
    Returns:
        pd.DataFrame: Results with max delta and top position.
    """
    results = []
    
    required_cols = ['chrom', 'pos', 'ref', 'alt']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"Input CSV must contain columns: {required_cols}")
        
    for index, row in df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        
        # Run prediction
        rel_coords, delta_rna = get_delta_rna_seq(chrom, pos, ref, alt)
        
        # Extract summary metrics
        max_abs_idx = np.argmax(np.abs(delta_rna))
        max_delta = delta_rna[max_abs_idx]
        max_pos = rel_coords[max_abs_idx]
        
        results.append({
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "max_delta_rna": round(max_delta, 4),
            "max_delta_pos_rel": max_pos
        })
        
    return pd.DataFrame(results)
