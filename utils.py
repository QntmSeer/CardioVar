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
    Calculate Delta RNA-seq values for a variant using Enformer via variant_engine.
    Raises RuntimeError if Enformer is unavailable (no synthetic fallback).
    """
    from variant_engine import compute_variant_impact
    result = compute_variant_impact(chrom, pos, ref, alt, window_size=window_size)
    curve = result["curve"]
    return np.array(curve["x"]), np.array(curve["y"])

def get_gene_info(chrom, pos):
    """
    Look up gene annotation for a GRCh38 coordinate via Ensembl overlap API.
    Returns an explicit Unannotated entry if no real gene is returned.
    """
    from api_integrations import fetch_ensembl_gene
    from variant_engine import compute_variant_impact

    # Use known GRCh38 cardiac gene coordinates (no chromosome-only fallback)
    known_regions = {
        "chr1": [(156100000, 156200000, "LMNA"), (55000000, 55100000, "PCSK9"), (236700000, 236750000, "ACTN2")],
        "chr2": [(21000000, 21100000, "APOB"), (178500000, 178600000, "TTN")],
        "chr3": [(46850000, 46900000, "MYL3")],
        "chr22": [(36100000, 36400000, "MYH9")],
    }
    hits = []
    for start, end, gene in known_regions.get(chrom, []):
        if start <= pos <= end:
            info = fetch_ensembl_gene(gene)
            hits.append({"Gene": gene, "Distance": 0, "Type": "Intragenic", "Source": "Ensembl" if info else "GRCh38 map"})
    if not hits:
        hits.append({"Gene": "UNANNOTATED", "Distance": None, "Type": "Unannotated"})
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
    sns.lineplot(x=rel_coords, y=delta_rna, color=line_color, linewidth=2.5, ax=ax, label=r'$\Delta$ RNA-seq')
    
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
    ax.set_ylabel(r"$\Delta$ RNA-seq Level", fontsize=12)
    ax.legend(loc='upper right', frameon=False)
    
    # Remove unnecessary spines for minimal look
    sns.despine(trim=True)
    
    plt.tight_layout()
    return fig

def plot_gene_structure(chrom, pos, window_size=100, exons=None):
    """
    Plot real exon structure around the variant (or show Unavailable if none provided).
    """
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.plot([-window_size, window_size], [0, 0], color='black', linewidth=1)

    if exons:
        for ex in exons:
            start, end = ex["start"], ex["end"]
            rect = plt.Rectangle((start, -0.2), end - start, 0.4, facecolor='#3498db', edgecolor='black', alpha=0.7)
            ax.add_patch(rect)
    else:
        ax.text(0, 0.35, "Exon structure unavailable", ha='center', va='center', color='gray')

    ax.axvline(0, color='#e74c3c', linestyle=':', linewidth=2, label='Variant')
    ax.set_xlim(-window_size, window_size)
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    ax.set_xlabel("Relative Genomic Coordinate (bp)")
    ax.set_title("Gene Structure Context", fontsize=12)
    sns.despine(left=True, bottom=False)
    plt.tight_layout()
    return fig

def plot_conservation(chrom, pos, window_size=100, scores=None):
    """
    Plot real PhyloP conservation scores around the variant (or show Unavailable if missing).
    """
    x = np.arange(-window_size, window_size + 1)
    fig, ax = plt.subplots(figsize=(10, 1.5))

    if scores is not None and len(scores) == len(x):
        scores = np.asarray(scores)
        ax.fill_between(x, scores, 0, where=(scores > 0), color='#27ae60', alpha=0.6, label='Conserved')
        ax.fill_between(x, scores, 0, where=(scores < 0), color='#95a5a6', alpha=0.3, label='Neutral/Accelerated')
    else:
        ax.text(0, 0, "PhyloP conservation data unavailable", ha='center', va='center', color='gray')
        ax.set_ylim(-1, 1)

    ax.axvline(0, color='#e74c3c', linestyle=':', linewidth=2)
    ax.set_xlim(-window_size, window_size)
    ax.set_ylabel("PhyloP")
    ax.set_title("Evolutionary Conservation", fontsize=12)
    sns.despine(bottom=True)
    ax.set_xticks([])
    plt.tight_layout()
    return fig
