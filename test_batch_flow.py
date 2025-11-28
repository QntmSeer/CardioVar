import requests
import pandas as pd
import time
import sys

API_URL = "http://localhost:8000"
CSV_FILE = "batch_smoketest.csv"

def test_batch_flow():
    print(f"Loading {CSV_FILE}...")
    try:
        df = pd.read_csv(CSV_FILE)
        variants = df.to_dict(orient="records")
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    print(f"Submitting batch of {len(variants)} variants...")
    try:
        resp = requests.post(f"{API_URL}/batch-start", json={"variants": variants})
        resp.raise_for_status()
        data = resp.json()
        batch_id = data["batch_id"]
        print(f"Batch started. ID: {batch_id}")
    except Exception as e:
        print(f"Error starting batch: {e}")
        return

    print("Polling for completion...")
    while True:
        try:
            status_resp = requests.get(f"{API_URL}/batch-status/{batch_id}")
            status_resp.raise_for_status()
            status_data = status_resp.json()
            
            state = status_data["status"]
            processed = status_data["processed"]
            total = status_data["total"]
            
            print(f"Status: {state} ({processed}/{total})")
            
            if state in ["completed", "failed"]:
                break
            
            time.sleep(1)
        except Exception as e:
            print(f"Error polling status: {e}")
            break

    if state == "completed":
        print("\nBatch completed successfully!")
        results = status_data["results"]
        res_df = pd.DataFrame(results)
        print("Result columns:", res_df.columns.tolist())
        print(res_df.head())
        print(f"\nTotal results: {len(res_df)}")
    else:
        print(f"\nBatch failed: {status_data.get('error')}")

if __name__ == "__main__":
    test_batch_flow()
