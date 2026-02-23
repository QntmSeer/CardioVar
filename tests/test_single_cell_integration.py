
import requests
import json
import os
import sys

# Add parent directory to path so we can import api_integrations
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Assume API is running locally on 8000 for verification or we mock it
# Since we can't easily start the server and keep it running in this environment, 
# we will test the underlying function directly and mocking the API call structure if possible.
# Actually, we can import the app and test via TestClient if fastapi.testclient is available, 
# or just import the function from api_integrations.

try:
    from api_integrations import fetch_single_cell_expression
    
    print("Testing fetch_single_cell_expression directly...")
    
    gene = "MYH9"
    data = fetch_single_cell_expression(gene)
    
    if data:
        print(f"SUCCESS: Fetched data for {gene}")
        print(json.dumps(data, indent=2))
        
        # Verify content
        expected_types = ["Cardiomyocytes", "Fibroblasts"]
        found_types = [d["cell_type"] for d in data]
        
        if all(t in found_types for t in expected_types):
             print("SUCCESS: Found expected cell types.")
        else:
             print(f"FAILURE: Missing expected cell types. Found: {found_types}")
             sys.exit(1)
    else:
        print(f"FAILURE: No data returned for {gene}")
        sys.exit(1)
        
    print("-" * 20)
    
    gene_missing = "NONEXISTENT"
    data_missing = fetch_single_cell_expression(gene_missing)
    if data_missing is None:
        print(f"SUCCESS: Correctly returned None for {gene_missing}")
    else:
        print(f"FAILURE: Returned data for non-existent gene: {data_missing}")
        sys.exit(1)

except ImportError:
    print("Could not import api_integrations. Make sure you are in the root directory.")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)
