import requests

def test_query(q):
    url = "https://myvariant.info/v1/query"
    params = {"q": q, "fields": "hg38", "size": 1}
    print(f"Querying: {q}")
    try:
        resp = requests.get(url, params=params)
        print(f"Status: {resp.status_code}")
        if resp.status_code == 200:
            print("Hits:", len(resp.json().get("hits", [])))
            print(resp.json())
        else:
            print("Error:", resp.text)
    except Exception as e:
        print(f"Exception: {e}")
    print("-" * 20)

# Test cases
test_query("cdk2") 
test_query("dbnsfp.alphamissense") # Check if this field exists
test_query("cadd") # Check for CADD scores
test_query("clinvar") # Check for ClinVar
