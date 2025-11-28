"""
Test the CardioVar API endpoints on localhost.
"""
import requests
import json

API_URL = "http://localhost:8000"

print("="*60)
print("CardioVar API Test")
print("="*60)

# Test 1: Variant Impact
print("\n1. Testing /variant-impact endpoint...")
print("   Variant: chr22:36191400 A->C")

payload = {
    "assembly": "GRCh38",
    "chrom": "chr22",
    "pos": 36191400,
    "ref": "A",
    "alt": "C",
    "window_size": 100
}

try:
    response = requests.post(f"{API_URL}/variant-impact", json=payload, timeout=120)
    response.raise_for_status()
    data = response.json()
    
    print("   [SUCCESS] Variant impact computed!")
    print(f"   Gene: {data['metrics']['gene_symbol']}")
    print(f"   Max Delta: {data['metrics']['max_delta']:.4f}")
    print(f"   gnomAD Freq: {data['metrics']['gnomad_freq']:.6f}")
    print(f"   Percentile: {data['metrics']['percentile']:.1f}%")
    print(f"   Exons found: {len(data['tracks']['exons'])}")
    
    if data['tracks']['exons']:
        print(f"   First exon: {data['tracks']['exons'][0]['id']}")
    
except Exception as e:
    print(f"   [ERROR] {e}")

# Test 2: Gene Annotations
print("\n2. Testing /gene-annotations endpoint...")
print("   Gene: MYH9")

try:
    response = requests.get(f"{API_URL}/gene-annotations", params={"gene": "MYH9"}, timeout=30)
    response.raise_for_status()
    data = response.json()
    
    print("   [SUCCESS] Gene annotations retrieved!")
    print(f"   Name: {data.get('name', 'N/A')}")
    print(f"   Chromosome: {data.get('chromosome', 'N/A')}")
    print(f"   Source: {data.get('source', 'N/A')}")
    
    if 'protein_domains' in data:
        print(f"   Protein length: {data.get('protein_length', 'N/A')} aa")
        print(f"   Protein domains: {len(data['protein_domains'])}")
        if data['protein_domains']:
            print(f"   First domain: {data['protein_domains'][0]['name']} ({data['protein_domains'][0]['type']})")
    
    if 'expression' in data:
        print(f"   Expression data: {len(data['expression'])} tissues")
    
except Exception as e:
    print(f"   [ERROR] {e}")

print("\n" + "="*60)
print("Test Complete")
print("="*60)
print("\nNow open http://localhost:8501 in your browser to test the dashboard!")
