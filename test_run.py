import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from dotenv import load_dotenv

# Load environment variables
load_dotenv()
API_KEY = os.getenv("ALPHAGENOME_API_KEY")

print(f"API Key loaded: {bool(API_KEY)}")
if API_KEY:
    print(f"API Key start: {API_KEY[:4]}...")

# Configure plotting
sns.set_theme(style="whitegrid")
plt.rcParams['figure.figsize'] = (12, 6)

def get_delta_rna_seq(chrom, pos, ref, alt, window_size=100):
    """
    Fetch or calculate Delta RNA-seq values for a variant.
    """
    # --- SYNTHETIC DATA GENERATION FOR DEMO ---
    x = np.arange(-window_size, window_size + 1)
    direction = np.random.choice([-1, 1]) 
    signal = 3.5 * direction * np.exp(-0.02 * (x)**2) 
    noise = np.random.normal(0, 0.3, len(x))
    delta_rna = signal + noise
    return x, delta_rna

def plot_variant_effect(chrom, pos, ref, alt):
    """
    Analyze and plot the effect of a variant.
    """
    print(f"Analyzing variant: {chrom}:{pos} {ref}->{alt}...")
    
    rel_coords, delta_rna = get_delta_rna_seq(chrom, pos, ref, alt)
    
    top_idx = np.argmax(np.abs(delta_rna))
    top_pos = rel_coords[top_idx]
    top_val = delta_rna[top_idx]
    
    plt.figure(figsize=(12, 6))
    plt.plot(rel_coords, delta_rna, label='Delta RNA-seq', color='#2c3e50', linewidth=2)
    plt.axhline(0, color='gray', linestyle='--', alpha=0.5)
    plt.axvline(0, color='#e74c3c', linestyle=':', label='Variant Position (0)')
    plt.scatter(top_pos, top_val, color='#e74c3c', s=100, zorder=5, label=f'Max Effect ({top_val:.2f})')
    
    plt.title(f"Variant Effect Prediction: {chrom}:{pos} {ref}->{alt}")
    plt.xlabel("Relative Genomic Coordinate (bp)")
    plt.ylabel("Delta RNA-seq Level")
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    # Save instead of show
    filename = f"variant_effect_{chrom}_{pos}.png"
    plt.savefig(filename)
    print(f"Plot saved to {filename}")

# Test Case
try:
    plot_variant_effect("chr22", 36191400, "A", "C")
    print("Test run completed successfully.")
except Exception as e:
    print(f"Test run failed: {e}")
