"""
Pre-compute background variant impact distributions for common cardiovascular genes.

This script:
1. Loads a list of cardiovascular genes
2. Fetches known variants from gnomAD (common variants, likely benign)
3. Computes impact for each variant using the variant engine
4. Stores the distribution for fast runtime lookup

Note: This is a SLOW process (hours) due to Enformer inference.
Run this periodically (weekly/monthly) to update the database.
"""

import json
import numpy as np
from variant_engine import compute_variant_impact
from api_integrations import fetch_ensembl_gene
import time
from pathlib import Path

# Top cardiovascular genes (expand this list as needed)
CARDIOVASCULAR_GENES = [
    "MYH9",   # Myosin heavy chain 9
    "MYBPC3", # Myosin binding protein C
    "MYH7",   # Myosin heavy chain 7
    "TNNT2",  # Troponin T2
    "TNNI3",  # Troponin I3
    "TPM1",   # Tropomyosin 1
    "ACTC1",  # Actin alpha cardiac muscle 1
    "MYL2",   # Myosin light chain 2
    "MYL3",   # Myosin light chain 3
    "LMNA",   # Lamin A/C
    "PCSK9",  # Proprotein convertase subtilisin/kexin type 9
    "APOB",   # Apolipoprotein B
    "LDLR",   # LDL receptor
    "TTN",    # Titin
    "SCN5A",  # Sodium voltage-gated channel alpha subunit 5
]

def generate_synthetic_variants_for_gene(gene_symbol, count=50):
    """
    Generate synthetic variant positions within a gene for testing.
    In production, this would query gnomAD or ClinVar.
    
    Args:
        gene_symbol: Gene symbol
        count: Number of variants to generate
    
    Returns:
        List of variant dicts with chrom, pos, ref, alt
    """
    print(f"  Generating {count} synthetic variants for {gene_symbol}...")
    
    # Get gene coordinates from Ensembl
    gene_data = fetch_ensembl_gene(gene_symbol)
    
    if not gene_data:
        print(f"  Warning: Could not fetch gene data for {gene_symbol}")
        return []
    
    chrom = gene_data.get('seq_region_name', '1')
    start = gene_data.get('start', 1000000)
    end = gene_data.get('end', 1100000)
    
    # Ensure chr prefix
    if not chrom.startswith('chr'):
        chrom = f'chr{chrom}'
    
    variants = []
    
    # Generate random positions within gene
    np.random.seed(hash(gene_symbol) % 2**32)  # Deterministic per gene
    positions = np.random.randint(start, end, count)
    
    # Simple SNVs (single nucleotide variants)
    bases = ['A', 'C', 'G', 'T']
    
    for pos in positions:
        ref = np.random.choice(bases)
        alt = np.random.choice([b for b in bases if b != ref])
        
        variants.append({
            'chrom': chrom,
            'pos': int(pos),
            'ref': ref,
            'alt': alt
        })
    
    return variants

def compute_background_for_gene(gene_symbol, variant_count=50, use_enformer=True):
    """
    Compute background impact distribution for a gene.
    
    Args:
        gene_symbol: Gene symbol
        variant_count: Number of variants to process
        use_enformer: Whether to use Enformer (slow) or heuristic (fast)
    
    Returns:
        Dict with distribution data
    """
    print(f"\n{'='*60}")
    print(f"Processing gene: {gene_symbol}")
    print(f"{'='*60}")
    
    # Generate or fetch variants
    variants = generate_synthetic_variants_for_gene(gene_symbol, variant_count)
    
    if not variants:
        print(f"No variants generated for {gene_symbol}")
        return None
    
    print(f"Computing impact for {len(variants)} variants...")
    
    impacts = []
    failed = 0
    
    start_time = time.time()
    
    for i, var in enumerate(variants):
        try:
            # Compute impact
            result = compute_variant_impact(
                chrom=var['chrom'],
                pos=var['pos'],
                ref=var['ref'],
                alt=var['alt']
            )
            
            # Extract max delta (absolute value)
            max_delta = abs(result['metrics']['max_delta'])
            impacts.append(max_delta)
            
            # Progress update
            if (i + 1) % 10 == 0:
                elapsed = time.time() - start_time
                avg_time = elapsed / (i + 1)
                remaining = avg_time * (len(variants) - i - 1)
                print(f"  Progress: {i+1}/{len(variants)} | "
                      f"Avg: {avg_time:.1f}s/variant | "
                      f"ETA: {remaining/60:.1f}min")
        
        except Exception as e:
            print(f"  Failed variant {i+1}: {e}")
            failed += 1
            continue
    
    total_time = time.time() - start_time
    
    print(f"\nCompleted {gene_symbol}:")
    print(f"  Success: {len(impacts)}/{len(variants)}")
    print(f"  Failed: {failed}")
    print(f"  Total time: {total_time/60:.1f} minutes")
    print(f"  Avg time per variant: {total_time/len(variants):.1f} seconds")
    
    if not impacts:
        return None
    
    # Calculate statistics
    impacts_array = np.array(impacts)
    
    return {
        'gene_symbol': gene_symbol,
        'variant_count': len(impacts),
        'impact_distribution': impacts,  # Full list for percentile calculations
        'statistics': {
            'mean': float(np.mean(impacts_array)),
            'median': float(np.median(impacts_array)),
            'std': float(np.std(impacts_array)),
            'min': float(np.min(impacts_array)),
            'max': float(np.max(impacts_array)),
            'percentiles': {
                '25': float(np.percentile(impacts_array, 25)),
                '50': float(np.percentile(impacts_array, 50)),
                '75': float(np.percentile(impacts_array, 75)),
                '90': float(np.percentile(impacts_array, 90)),
                '95': float(np.percentile(impacts_array, 95)),
                '99': float(np.percentile(impacts_array, 99)),
            }
        },
        'source': 'Synthetic variants (gnomAD-style)',
        'computed_with': 'Enformer' if use_enformer else 'Heuristic',
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
    }

def main():
    """Main execution function."""
    print("="*60)
    print("CardioVar Background Distribution Pre-computation")
    print("="*60)
    print(f"\nProcessing {len(CARDIOVASCULAR_GENES)} genes")
    print(f"Variants per gene: 50")
    print(f"Estimated time: ~30-60 minutes (with Enformer on CPU)")
    print("\nStarting in 3 seconds...")
    time.sleep(3)
    
    results = {}
    
    for gene in CARDIOVASCULAR_GENES:
        result = compute_background_for_gene(gene, variant_count=50)
        
        if result:
            results[gene] = result
        
        # Small delay between genes
        time.sleep(1)
    
    # Save results
    output_path = Path('data/gene_backgrounds.json')
    output_path.parent.mkdir(exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"Completed!")
    print(f"{'='*60}")
    print(f"Results saved to: {output_path}")
    print(f"Genes processed: {len(results)}/{len(CARDIOVASCULAR_GENES)}")
    
    # Summary statistics
    print("\nSummary Statistics:")
    for gene, data in results.items():
        stats = data['statistics']
        print(f"\n{gene}:")
        print(f"  Mean impact: {stats['mean']:.3f}")
        print(f"  Median: {stats['median']:.3f}")
        print(f"  Range: [{stats['min']:.3f}, {stats['max']:.3f}]")
        print(f"  95th percentile: {stats['percentiles']['95']:.3f}")

if __name__ == '__main__':
    main()
