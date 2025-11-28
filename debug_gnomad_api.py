# Fix Windows console encoding
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import requests
import json

# Test gnomAD GraphQL API directly
graphql_url = "https://gnomad.broadinstitute.org/api"

# Test with LMNA pathogenic variant
variant_id = "1-156137666-C-T"  # LMNA c.1621C>T

query = """
query VariantClinVar($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variantId
    clinvar {
      clinicalSignificance
      goldStars
      reviewStatus
      variationId
      releaseDate
    }
  }
}
"""

# Try different datasets
datasets = ["gnomad_r4", "gnomad_r3", "gnomad_r2_1"]

print("Testing gnomAD GraphQL API directly...")
print(f"Variant: {variant_id}")
print("="*60)

for dataset in datasets:
    print(f"\nTrying dataset: {dataset}")
    
    variables = {
        "variantId": variant_id,
        "datasetId": dataset
    }
    
    try:
        resp = requests.post(
            graphql_url,
            json={"query": query, "variables": variables},
            timeout=15,
            headers={"Content-Type": "application/json"}
        )
        
        print(f"Status: {resp.status_code}")
        
        if resp.status_code == 200:
            data = resp.json()
            
            if "errors" in data:
                print(f"GraphQL Errors: {json.dumps(data['errors'], indent=2)}")
            
            if "data" in data:
                variant_data = data.get("data", {}).get("variant")
                if variant_data:
                    print(f"✅ Variant found!")
                    print(f"Variant ID: {variant_data.get('variantId')}")
                    
                    clinvar = variant_data.get("clinvar")
                    if clinvar:
                        print(f"✅ ClinVar data found!")
                        print(json.dumps(clinvar, indent=2))
                    else:
                        print(f"⚠️  No ClinVar data for this variant")
                else:
                    print(f"⚠️  Variant not found in {dataset}")
        else:
            print(f"Error: {resp.text[:200]}")
            
    except Exception as e:
        print(f"Exception: {e}")

# Also test the variant query without ClinVar to see if variant exists
print("\n" + "="*60)
print("Testing if variant exists in gnomAD (without ClinVar)...")
print("="*60)

simple_query = """
query SimpleVariant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variantId
    chrom
    pos
    ref
    alt
    genome {
      ac
      an
      af
    }
  }
}
"""

for dataset in ["gnomad_r4", "gnomad_r3"]:
    print(f"\nDataset: {dataset}")
    variables = {
        "variantId": variant_id,
        "datasetId": dataset
    }
    
    try:
        resp = requests.post(
            graphql_url,
            json={"query": simple_query, "variables": variables},
            timeout=15,
            headers={"Content-Type": "application/json"}
        )
        
        if resp.status_code == 200:
            data = resp.json()
            if "data" in data and data["data"].get("variant"):
                print(f"✅ Variant exists!")
                print(json.dumps(data["data"]["variant"], indent=2))
            else:
                print(f"⚠️  Variant not in {dataset}")
                if "errors" in data:
                    print(f"Errors: {data['errors']}")
    except Exception as e:
        print(f"Exception: {e}")
