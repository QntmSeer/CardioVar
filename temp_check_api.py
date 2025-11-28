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

reset_fallback_flag()
print('--- Testing default (cached/fallback) calls ---')
print('gnomAD:', fetch_gnomad_frequency('chr22', 36191400, 'A', 'C', force_live=False))
print('Ensembl gene MYH9:', fetch_ensembl_gene('MYH9', force_live=False))
print('PhyloP:', fetch_ucsc_phylop('chr22', 36191300, 36191500, force_live=False)[:5])
print('GTEx MYH9:', fetch_gtex_expression('MYH9', force_live=False)[:3])
print('Gene structure:', fetch_gene_structure('chr22', 36191400, force_live=False)[:2])
print('Protein domains MYH9:', fetch_protein_domains('MYH9', force_live=False))

reset_fallback_flag()
print('\n--- Testing force_live=True calls ---')
print('gnomAD live:', fetch_gnomad_frequency('chr22', 36191400, 'A', 'C', force_live=True))
print('Ensembl gene live:', fetch_ensembl_gene('MYH9', force_live=True))
print('PhyloP live:', fetch_ucsc_phylop('chr22', 36191300, 36191500, force_live=True)[:5])
print('GTEx live:', fetch_gtex_expression('MYH9', force_live=True)[:3])
print('Gene structure live:', fetch_gene_structure('chr22', 36191400, force_live=True)[:2])
print('Protein domains live:', fetch_protein_domains('MYH9', force_live=True))
