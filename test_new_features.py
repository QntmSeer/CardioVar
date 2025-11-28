from api_integrations import fetch_clinvar_variants, fetch_dbsnp_variants
import time

print("Testing direct function calls...")
try:
    print("Fetching ClinVar...")
    res = fetch_clinvar_variants("chr22", 36191400, "A", "C", force_live=True)
    print("ClinVar result:", res)
except Exception as e:
    import traceback
    traceback.print_exc()

try:
    print("Fetching dbSNP...")
    res = fetch_dbsnp_variants("chr22", 36191400, "A", "C", force_live=True)
    print("dbSNP result:", res)
except Exception as e:
    import traceback
    traceback.print_exc()
