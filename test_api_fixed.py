# test_api_fixed.py
import json
import logging
from api_integrations import (
    fetch_gnomad_frequency,
    fetch_ucsc_phylop,
    fetch_protein_domains,
    fetch_gtex_expression,
    reset_fallback_flag,
    fallback_used
)
import api_integrations

# Configure logging to show debug output
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')

def pretty_print(label, data):
    print(f"\n=== {label} ===")
    if isinstance(data, (dict, list)):
        print(json.dumps(data, indent=2))
    else:
        print(data)

def main():
    chrom, pos, ref, alt = "chr22", 36191400, "A", "C"
    gene = "MYH9"

    print("\n\n>>> TESTING DEFAULT MODE (force_live=False) - Should use fallback if API fails <<<")
    reset_fallback_flag()
    
    # gnomAD frequency
    freq = fetch_gnomad_frequency(chrom, pos, ref, alt, force_live=False)
    pretty_print("gnomAD frequency (Fallback Allowed)", freq)
    print(f"Fallback used: {api_integrations.fallback_used}")

    # PhyloP scores (±100 bp)
    phylop = fetch_ucsc_phylop(chrom, pos - 100, pos + 100, force_live=False)
    pretty_print("UCSC PhyloP (200 bp window) (Fallback Allowed)", phylop[:10] if phylop else None)

    # Protein domains
    domains = fetch_protein_domains(gene, force_live=False)
    pretty_print("Protein domains (MYH9) (Fallback Allowed)", domains)

    # GTEx expression
    expr = fetch_gtex_expression(gene, force_live=False)
    pretty_print("GTEx expression (MYH9) (Fallback Allowed)", expr)

    print("\n\n>>> TESTING FORCE LIVE MODE (force_live=True) - Should return None if API fails <<<")
    reset_fallback_flag()
    
    # gnomAD frequency
    freq_live = fetch_gnomad_frequency(chrom, pos, ref, alt, force_live=True)
    pretty_print("gnomAD frequency (Force Live)", freq_live)
    print(f"Fallback used: {api_integrations.fallback_used}")

    # PhyloP scores
    phylop_live = fetch_ucsc_phylop(chrom, pos - 100, pos + 100, force_live=True)
    pretty_print("UCSC PhyloP (Force Live)", phylop_live[:10] if phylop_live else None)

if __name__ == "__main__":
    main()
