from fastapi.testclient import TestClient
from api import app
import time
import json

client = TestClient(app)

# 1. Test System Status
print("--- Testing System Status ---")
try:
    resp = client.get("/system-status")
    if resp.status_code == 200:
        print("System Status:", json.dumps(resp.json(), indent=2))
    else:
        print("Failed to get system status:", resp.text)
except Exception as e:
    print(f"System status check failed: {e}")

# 2. Test Async Batch Processing
print("\n--- Testing Async Batch Processing ---")
variants = [
    {"chrom": "chr22", "pos": 36191400, "ref": "A", "alt": "C"},
    {"chrom": "chr1", "pos": 155236508, "ref": "A", "alt": "G"},
    {"chrom": "chrX", "pos": 1000, "ref": "A", "alt": "T"}
]

try:
    # Start Batch
    print("Starting batch...")
    start_resp = client.post("/batch-start", json={"variants": variants})
    if start_resp.status_code != 200:
        print("Failed to start batch:", start_resp.text)
        exit(1)
        
    batch_data = start_resp.json()
    batch_id = batch_data["batch_id"]
    print(f"Batch started with ID: {batch_id}")
    
    # Poll Status
    # Note: TestClient runs synchronous code. Background tasks might not run automatically 
    # unless we explicitly trigger them or if TestClient handles them (it usually does run them after response).
    # However, since process_batch_task updates a global dict, we can check it.
    
    # In TestClient, background tasks are executed after the response is sent.
    # So by the time we get the response, the task might have started or finished depending on execution.
    # But since it's synchronous in this context (single thread), it might run sequentially.
    
    # Let's check status immediately
    status_resp = client.get(f"/batch-status/{batch_id}")
    print("Initial Status:", status_resp.json())
    
    # If it's already done (because TestClient runs tasks synchronously), great.
    # If not, we might need to wait, but TestClient usually blocks until tasks are done.
    
    status_data = status_resp.json()
    if status_data["status"] == "completed":
        print("\nFinal Result:")
        print(json.dumps(status_data, indent=2))
    else:
        print("Batch not completed immediately (unexpected for TestClient default behavior)")

except Exception as e:
    print(f"Batch test failed: {e}")
