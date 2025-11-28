from api import batch_impact, BatchRequest, VariantRequest
import json

variants = [
    {"chrom": "chr22", "pos": 36191400, "ref": "A", "alt": "C"},
    {"chrom": "chr1", "pos": 155236508, "ref": "A", "alt": "G"},
    {"chrom": "chrX", "pos": 1000, "ref": "A", "alt": "T"}
]

# Create Pydantic models
req_variants = [VariantRequest(**v) for v in variants]
batch_req = BatchRequest(variants=req_variants)

print("Running batch_impact directly...")
try:
    results = batch_impact(batch_req)
    print("Success!")
    print(json.dumps(results, indent=2))
except Exception as e:
    print(f"Function failed: {e}")
    import traceback
    traceback.print_exc()
