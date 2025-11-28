import time
import logging
from api_integrations import (
    fetch_gnomad_frequency,
    fetch_ensembl_gene,
    fetch_ucsc_phylop,
    fetch_gtex_expression,
    fetch_gene_structure,
    fetch_protein_domains,
    reset_fallback_flag,
)

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

chrom, pos, ref, alt = "chr22", 36191400, "A", "C"
gene = "MYH9"

functions = [
    ("gnomAD", lambda: fetch_gnomad_frequency(chrom, pos, ref, alt, force_live=False)),
    ("Ensembl", lambda: fetch_ensembl_gene(gene, force_live=False)),
    ("PhyloP", lambda: fetch_ucsc_phylop(chrom, pos-100, pos+100, force_live=False)),
    ("GTEx", lambda: fetch_gtex_expression(gene, force_live=False)),
    ("GeneStructure", lambda: fetch_gene_structure(chrom, pos, force_live=False)),
    ("ProteinDomains", lambda: fetch_protein_domains(gene, force_live=False)),
]

def time_call(func):
    start = time.time()
    result = func()
    duration = time.time() - start
    return result, duration

print("--- Warm up cache (first live calls) ---")
reset_fallback_flag()
for name, fn in functions:
    _, t = time_call(fn)
    print(f"{name} warm-up: {t:.4f}s")

print("\n--- Cached calls (force_live=False) ---")
reset_fallback_flag()
for name, fn in functions:
    _, t = time_call(fn)
    print(f"{name} cached: {t:.4f}s")

print("\n--- Live calls (force_live=True) ---")
reset_fallback_flag()
for name, fn in [
    ("gnomAD", lambda: fetch_gnomad_frequency(chrom, pos, ref, alt, force_live=True)),
    ("Ensembl", lambda: fetch_ensembl_gene(gene, force_live=True)),
    ("PhyloP", lambda: fetch_ucsc_phylop(chrom, pos-100, pos+100, force_live=True)),
    ("GTEx", lambda: fetch_gtex_expression(gene, force_live=True)),
    ("GeneStructure", lambda: fetch_gene_structure(chrom, pos, force_live=True)),
    ("ProteinDomains", lambda: fetch_protein_domains(gene, force_live=True)),
]:
    _, t = time_call(fn)
    print(f"{name} live: {t:.4f}s")
