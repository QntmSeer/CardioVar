import time, logging
from api_integrations import fetch_gnomad_frequency, reset_fallback_flag
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

# Known variant present in gnomAD (ACE gene rs699)
chrom = 'chr1'
pos = 155236508
ref = 'A'
alt = 'G'

reset_fallback_flag()
print('--- First call (force_live=True) ---')
start = time.time()
freq = fetch_gnomad_frequency(chrom, pos, ref, alt, force_live=True)
print('Frequency:', freq)
print('Time:', time.time() - start)

reset_fallback_flag()
print('\n--- Second call (force_live=False, should hit cache) ---')
start = time.time()
freq2 = fetch_gnomad_frequency(chrom, pos, ref, alt, force_live=False)
print('Frequency:', freq2)
print('Time:', time.time() - start)
