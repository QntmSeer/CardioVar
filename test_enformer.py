import logging
import sys
import os

# Configure logging
logging.basicConfig(level=logging.INFO)

def test_enformer():
    print("Testing Enformer Integration...")
    
    try:
        import torch
        print(f"Torch version: {torch.__version__}")
    except ImportError:
        print("❌ Torch not installed.")
        return

    try:
        import enformer_pytorch
        print("Enformer-pytorch installed.")
    except ImportError:
        print("❌ enformer-pytorch not installed.")
        return

    try:
        from enformer_wrapper import predict_variant_impact_dl
        
        # Test with a simple variant (PCSK9)
        # chr1:55039974 G>A (rs11591147 - Loss of function)
        print("\nAttempting prediction for PCSK9 variant...")
        result = predict_variant_impact_dl("chr1", 55039974, "G", "A")
        
        if result:
            print("SUCCESS: Enformer Prediction SUCCESS!")
            print(f"Max Impact: {result['max_impact']}")
        else:
            print("FAILURE: Enformer Prediction FAILED (returned None).")
            
    except Exception as e:
        print(f"ERROR during Enformer test: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_enformer()
