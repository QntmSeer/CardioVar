import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from dotenv import load_dotenv

# Load environment variables
load_dotenv()
API_KEY = os.getenv("ALPHAGENOME_API_KEY")

# Configure plotting style
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams['figure.figsize'] = (12, 6)
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

def get_delta_rna_seq(chrom, pos, ref, alt, window_size=100):
    """
    Fetch or calculate Delta RNA-seq values for a variant.
    
    Args:
        chrom (str): Chromosome (e.g., 'chr22')
        pos (int): Position
        ref (str): Reference allele
        alt (str): Alternative allele
        window_size (int): Window size around variant to analyze
        
    Returns:
        tuple: (relative_coordinates, delta_values)
    """
    # TODO: Replace with actual alphagenome API call
    # Example: prediction = alphagenome.predict(chrom, pos, ref, alt)
    
    # --- SYNTHETIC DATA GENERATION ---
    np.random.seed(pos % 1000) # Seed for reproducibility based on position
    x = np.arange(-window_size, window_size + 1)
    
    # Generate a synthetic signal: a peak near the variant plus noise
    direction = 1 if (pos % 2 == 0) else -1
    signal = 3.5 * direction * np.exp(-0.02 * (x)**2) 
    
    # Add some random noise
    noise = np.random.normal(0, 0.3, len(x))
    
    delta_rna = signal + noise
    # ------------------------------------------
    
    return x, delta_rna

def get_gene_info(chrom, pos):
    """
    Mock function to get gene info. 
    In production, use pyensembl.
    """
    # Mock database
    mock_genes = {
        "chr22": [("MYH9", 36190000, 36200000), ("APOL1", 36200000, 36210000)],
        "chr1": [("PCSK9", 12345000, 12350000)]
    }
    
    hits = []
    if chrom in mock_genes:
        for gene, start, end in mock_genes[chrom]:
            dist = min(abs(pos - start), abs(pos - end))
            if start <= pos <= end:
                hits.append({"Gene": gene, "Distance": 0, "Type": "Intragenic"})
            else:
                hits.append({"Gene": gene, "Distance": dist, "Type": "Intergenic"})
    
    if not hits:
        hits.append({"Gene": "Unknown", "Distance": -1, "Type": "N/A"})
        
    return hits

def plot_deltas(rel_coords, delta_rna, chrom, pos, ref, alt, line_color='#2c3e50', highlight_color='#e74c3c'):
    """
    Create a matplotlib figure for the variant effect.
    """
    # Identify Top Hits
    top_idx = np.argmax(np.abs(delta_rna))
    top_pos = rel_coords[top_idx]
    top_val = delta_rna[top_idx]
    
    fig, ax = plt.subplots(figsize=(10, 5))
    
    # Main signal line
    sns.lineplot(x=rel_coords, y=delta_rna, color=line_color, linewidth=2.5, ax=ax, label='$\Delta$ RNA-seq')
    
    # Reference lines
    ax.axhline(0, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(0, color=highlight_color, linestyle=':', alpha=0.8, label='Variant (0)')
    
    # Highlight top hit
    ax.scatter(top_pos, top_val, color=highlight_color, s=150, zorder=5, edgecolor='white', linewidth=1.5)
    
    # Annotation
    ax.annotate(f'Max: {top_val:.2f}',
                 (top_pos, top_val),
                 xytext=(10, 10), textcoords='offset points',
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=highlight_color, alpha=0.9),
                 arrowprops=dict(arrowstyle='->', connectionstyle="arc3,rad=.2", color=highlight_color))
    
    # Styling - Minimal
    ax.set_title(f"Variant Impact: {chrom}:{pos} {ref}→{alt}", fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel("Relative Genomic Coordinate (bp)", fontsize=12)
    ax.set_ylabel("$\Delta$ RNA-seq Level", fontsize=12)
    ax.legend(loc='upper right', frameon=False)
    
    # Remove unnecessary spines for minimal look
    sns.despine(trim=True)
    
    plt.tight_layout()
    return fig

def plot_gene_structure(chrom, pos, window_size=100):
    """
    Plot a mock gene structure (exons/introns) around the variant.
    """
    fig, ax = plt.subplots(figsize=(10, 2))
    
    # Mock Exons (relative coordinates)
    # Let's say there's an exon from -50 to +20
    exons = [(-50, 70)] 
    
    # Draw Intron Line
    ax.plot([-window_size, window_size], [0, 0], color='black', linewidth=1)
    
    # Draw Exons
    for start, length in exons:
        rect = plt.Rectangle((start, -0.2), length, 0.4, facecolor='#3498db', edgecolor='black', alpha=0.7)
        ax.add_patch(rect)
        ax.text(start + length/2, 0, "Exon 1", ha='center', va='center', color='white', fontweight='bold')
        
    # Variant Marker
    ax.axvline(0, color='#e74c3c', linestyle=':', linewidth=2, label='Variant')
    
    # Styling
    ax.set_xlim(-window_size, window_size)
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    ax.set_xlabel("Relative Genomic Coordinate (bp)")
    ax.set_title("Gene Structure Context (Mock)", fontsize=12)
    sns.despine(left=True, bottom=False)
    
    plt.tight_layout()
    return fig

def plot_conservation(chrom, pos, window_size=100):
    """
    Plot mock conservation scores (PhyloP) around the variant.
    """
    x = np.arange(-window_size, window_size + 1)
    # Generate mock conservation scores (noisy positive values indicate conservation)
    np.random.seed(pos)
    scores = np.random.normal(0.5, 1.0, len(x))
    scores[window_size-10:window_size+10] += 2.0 # Conserved region around variant
    
    fig, ax = plt.subplots(figsize=(10, 1.5))
    
    # Area plot
    ax.fill_between(x, scores, 0, where=(scores>0), color='#27ae60', alpha=0.6, label='Conserved')
    ax.fill_between(x, scores, 0, where=(scores<0), color='#95a5a6', alpha=0.3, label='Neutral/Accelerated')
    
    # Variant Marker
    ax.axvline(0, color='#e74c3c', linestyle=':', linewidth=2)
    
    # Styling
    ax.set_xlim(-window_size, window_size)
    ax.set_ylabel("PhyloP")
    ax.set_title("Evolutionary Conservation", fontsize=12)
    sns.despine(bottom=True)
    ax.set_xticks([]) # Hide x-axis ticks as it aligns with above plots
    
    plt.tight_layout()
    return fig
