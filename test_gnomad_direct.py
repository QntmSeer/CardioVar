# Fix Windows console encoding
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import requests
import json

# Test gnomAD GraphQL directly
graphql_url = "https://gnomad.broadinstitute.org/api"

query = """
query GeneInfo($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
  gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
    gene_id
    symbol
    name
    chrom
    start
    stop
  }
}
"""

variables = {
    "geneSymbol": "SCN5A",
    "referenceGenome": "GRCh38"
}

print("Testing gnomAD GraphQL for SCN5A...")
print("="*60)

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
            print(f"❌ GraphQL Errors:")
            print(json.dumps(data["errors"], indent=2))
        
        if "data" in data:
            gene_data = data.get("data", {}).get("gene")
            if gene_data:
                print(f"\n✅ Gene Data Found!")
                print(json.dumps(gene_data, indent=2))
            else:
                print(f"\n⚠️  No gene data")
    else:
        print(f"Error: {resp.text[:200]}")
        
except Exception as e:
    print(f"Exception: {e}")
