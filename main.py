import argparse
import matplotlib.pyplot as plt
from utils import get_delta_rna_seq, plot_deltas, get_gene_info

def main():
    parser = argparse.ArgumentParser(description="CardioVar: Variant Effect Prediction CLI")
    parser.add_argument("--chrom", type=str, required=True, help="Chromosome (e.g., chr22)")
    parser.add_argument("--pos", type=int, required=True, help="Position")
    parser.add_argument("--ref", type=str, required=True, help="Reference allele")
    parser.add_argument("--alt", type=str, required=True, help="Alternate allele")
    parser.add_argument("--output", type=str, default="output.png", help="Output filename for plot")
    
    args = parser.parse_args()
    
    print(f"Running analysis for {args.chrom}:{args.pos} {args.ref}->{args.alt}")
    
    # 1. Get Data
    rel_coords, delta_rna = get_delta_rna_seq(args.chrom, args.pos, args.ref, args.alt)
    
    # 2. Get Genes
    genes = get_gene_info(args.chrom, args.pos)
    print("Nearby Genes:")
    for g in genes:
        print(f" - {g['Gene']} ({g['Type']}, Dist: {g['Distance']})")
    
    # 3. Plot
    fig = plot_deltas(rel_coords, delta_rna, args.chrom, args.pos, args.ref, args.alt)
    fig.savefig(args.output)
    print(f"Plot saved to {args.output}")

if __name__ == "__main__":
    main()
