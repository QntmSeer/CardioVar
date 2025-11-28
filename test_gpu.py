"""
Test GPU availability and benchmark Enformer performance.
"""
import torch
import time

print("="*60)
print("CardioVar GPU/CUDA Test")
print("="*60)

# 1. Check PyTorch installation
print(f"\nPyTorch version: {torch.__version__}")

# 2. Check CUDA availability
cuda_available = torch.cuda.is_available()
print(f"CUDA available: {cuda_available}")

if cuda_available:
    print(f"CUDA version: {torch.version.cuda}")
    print(f"GPU count: {torch.cuda.device_count()}")
    print(f"Current GPU: {torch.cuda.current_device()}")
    print(f"GPU name: {torch.cuda.get_device_name(0)}")
    
    # 3. Test GPU computation
    print("\nTesting GPU computation...")
    try:
        x = torch.randn(1000, 1000).cuda()
        y = torch.randn(1000, 1000).cuda()
        
        start = time.time()
        z = x @ y
        torch.cuda.synchronize()  # Wait for GPU to finish
        gpu_time = time.time() - start
        
        print(f"GPU matrix multiplication: {gpu_time*1000:.2f}ms")
        print("[SUCCESS] GPU computation works!")
        
    except Exception as e:
        print(f"[ERROR] GPU test failed: {e}")
else:
    print("\n[INFO] No GPU detected. Enformer will run on CPU.")
    print("To enable GPU:")
    print("  1. Install CUDA Toolkit: https://developer.nvidia.com/cuda-downloads")
    print("  2. Reinstall PyTorch: pip install torch --index-url https://download.pytorch.org/whl/cu118")

# 4. Benchmark Enformer (if available)
print("\n" + "="*60)
print("Enformer Performance Test")
print("="*60)

try:
    from variant_engine import compute_variant_impact
    
    print("\nRunning test variant (chr22:36191400 A->C)...")
    print("This will load Enformer and run one prediction.")
    print("Please wait...\n")
    
    start_time = time.time()
    
    result = compute_variant_impact(
        chrom="chr22",
        pos=36191400,
        ref="A",
        alt="C"
    )
    
    total_time = time.time() - start_time
    
    print(f"\n{'='*60}")
    print("Results:")
    print(f"{'='*60}")
    print(f"Total time: {total_time:.2f} seconds")
    print(f"Max Delta: {result['metrics']['max_delta']:.4f}")
    print(f"Gene: {result['metrics']['gene_symbol']}")
    
    if cuda_available:
        print(f"\n[SUCCESS] GPU-accelerated prediction complete!")
        print(f"Expected speedup vs CPU: 10-30x")
    else:
        print(f"\n[INFO] CPU prediction complete.")
        print(f"With GPU, this would take <1 second instead of {total_time:.1f}s")
    
except Exception as e:
    print(f"[ERROR] Enformer test failed: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*60)
print("Test Complete")
print("="*60)
