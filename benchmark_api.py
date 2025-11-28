import time
from api_integrations import (
    fetch_gnomad_frequency,
    fetch_ensembl_gene,
    fetch_ucsc_phylop,
    fetch_gtex_expression,
    fetch_gene_structure,
    fetch_protein_domains,
    reset_fallback_flag,
)

chrom, pos, ref, alt = "chr22", 36191400, "A", "C"
gene = "MYH9"

def benchmark(func, *args, **kwargs):
    start = time.time()
    result = func(*args, **kwargs)
    duration = time.time() - start
    return result, duration

# Warm up cache
reset_fallback_flag()
fetch_gnomad_frequency(chrom, pos, ref, alt, force_live=False)
fetch_ensembl_gene(gene, force_live=False)
fetch_ucsc_phylop(chrom, pos-100, pos+100, force_live=False)
fetch_gtex_expression(gene, force_live=False)
fetch_gene_structure(chrom, pos, force_live=False)
fetch_protein_domains(gene, force_live=False)

# Benchmark first call (cache miss)
reset_fallback_flag()
_, t1 = benchmark(fetch_gnomad_frequency, chrom, pos, ref, alt, force_live=False)
_, t2 = benchmark(fetch_ensembl_gene, gene, force_live=False)
_, t3 = benchmark(fetch_ucsc_phylop, chrom, pos-100, pos+100, force_live=False)
_, t4 = benchmark(fetch_gtex_expression, gene, force_live=False)
_, t5 = benchmark(fetch_gene_structure, chrom, pos, force_live=False)
_, t6 = benchmark(fetch_protein_domains, gene, force_live=False)

# Benchmark second call (cache hit)
reset_fallback_flag()
_, t1c = benchmark(fetch_gnomad_frequency, chrom, pos, ref, alt, force_live=False)
_, t2c = benchmark(fetch_ensembl_gene, gene, force_live=False)
_, t3c = benchmark(fetch_ucsc_phylop, chrom, pos-100, pos+100, force_live=False)
_, t4c = benchmark(fetch_gtex_expression, gene, force_live=False)
_, t5c = benchmark(fetch_gene_structure, chrom, pos, force_live=False)
_, t6c = benchmark(fetch_protein_domains, gene, force_live=False)

print("--- Timing (seconds) ---")
print(f"gnomAD first: {t1:.4f}, second (cached): {t1c:.4f}")
print(f"Ensembl first: {t2:.4f}, second (cached): {t2c:.4f}")
print(f"PhyloP first: {t3:.4f}, second (cached): {t3c:.4f}")
print(f"GTEx first: {t4:.4f}, second (cached): {t4c:.4f}")
print(f"Gene structure first: {t5:.4f}, second (cached): {t5c:.4f}")
print(f"Protein domains first: {t6:.4f}, second (cached): {t6c:.4f}")
