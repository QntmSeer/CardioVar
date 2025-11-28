# GPU/CUDA Setup Guide for CardioVar

## Overview

Setting up GPU acceleration for Enformer will reduce inference time from **~30 seconds to <1 second** per variant (10-30x speedup).

---

## Prerequisites

### 1. Check GPU Availability

Run this command to verify you have an NVIDIA GPU:

```powershell
nvidia-smi
```

**Expected output** (if GPU is available):
```
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 535.xx       Driver Version: 535.xx       CUDA Version: 12.x    |
|-------------------------------+----------------------+----------------------+
| GPU  Name            TCC/WDDM | Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  NVIDIA GeForce ...  WDDM | 00000000:01:00.0  On |                  N/A |
| ...
```

If you see "command not found" or no output, you either:
- Don't have an NVIDIA GPU
- Don't have NVIDIA drivers installed

---

## Installation Steps

### Step 1: Install CUDA Toolkit

1. **Download CUDA Toolkit 11.8** (recommended for PyTorch compatibility):
   - Visit: https://developer.nvidia.com/cuda-11-8-0-download-archive
   - Select: Windows → x86_64 → 11 → exe (local)

2. **Run the installer**:
   - Choose "Custom" installation
   - Select: CUDA Toolkit, Visual Studio Integration
   - Uncheck: GeForce Experience, Driver (if already installed)

3. **Verify installation**:
   ```powershell
   nvcc --version
   ```
   Should show: `Cuda compilation tools, release 11.8`

### Step 2: Reinstall PyTorch with CUDA Support

```powershell
# Uninstall CPU-only PyTorch
pip uninstall torch

# Install CUDA-enabled PyTorch
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

**Note**: This will download ~2GB. Make sure you have good internet connection.

### Step 3: Verify PyTorch CUDA

Create a test file `test_cuda.py`:

```python
import torch

print(f"PyTorch version: {torch.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")

if torch.cuda.is_available():
    print(f"CUDA version: {torch.version.cuda}")
    print(f"GPU count: {torch.cuda.device_count()}")
    print(f"GPU name: {torch.cuda.get_device_name(0)}")
    
    # Test tensor on GPU
    x = torch.randn(1000, 1000).cuda()
    y = torch.randn(1000, 1000).cuda()
    z = x @ y
    print(f"GPU computation successful!")
else:
    print("CUDA not available. Check installation.")
```

Run it:
```powershell
python test_cuda.py
```

**Expected output**:
```
PyTorch version: 2.x.x+cu118
CUDA available: True
CUDA version: 11.8
GPU count: 1
GPU name: NVIDIA GeForce RTX ...
GPU computation successful!
```

---

## Step 4: Enable GPU in CardioVar

The code is already GPU-ready! The `enformer_wrapper.py` will automatically detect and use GPU if available.

**Current code** (already in place):
```python
def get_model():
    global _MODEL
    if _MODEL is None:
        print(">> Loading Enformer model...")
        _MODEL = Enformer.from_pretrained('EleutherAI/enformer-official-rough')
        
        # Automatically move to GPU if available
        if torch.cuda.is_available():
            _MODEL = _MODEL.cuda()
            print(f">> Using GPU: {torch.cuda.get_device_name(0)}")
        else:
            print(">> Using CPU")
        
        _MODEL.eval()
    return _MODEL
```

**No code changes needed!** Just install CUDA and PyTorch as above.

---

## Step 5: Benchmark Performance

Run the test to see the speedup:

```powershell
python test_direct.py
```

**Before GPU** (CPU):
```
>> Running Enformer inference...
[Takes ~20-30 seconds]
```

**After GPU**:
```
>> Using GPU: NVIDIA GeForce RTX 3080
>> Running Enformer inference...
[Takes <1 second]
```

---

## Troubleshooting

### Issue: "CUDA out of memory"

**Solution**: Enformer is a large model (~250M parameters). If you get OOM errors:

1. **Reduce batch size** (already 1 in our code)
2. **Use mixed precision**:
   ```python
   with torch.cuda.amp.autocast():
       ref_pred = model(ref_tensor)
       alt_pred = model(alt_tensor)
   ```
3. **Clear cache** between predictions:
   ```python
   torch.cuda.empty_cache()
   ```

### Issue: "CUDA driver version is insufficient"

**Solution**: Update NVIDIA drivers:
- Visit: https://www.nvidia.com/Download/index.aspx
- Download latest driver for your GPU
- Install and restart

### Issue: PyTorch still uses CPU

**Solution**: Verify installation:
```powershell
pip show torch
```

Should show: `Version: 2.x.x+cu118`

If it shows `2.x.x+cpu`, reinstall:
```powershell
pip uninstall torch
pip cache purge
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

---

## Performance Comparison

| Operation | CPU (Intel i7) | GPU (RTX 3080) | Speedup |
|-----------|----------------|----------------|---------|
| Model loading | ~8s | ~3s | 2.7x |
| Single inference | ~25s | ~0.8s | 31x |
| Batch (10 variants) | ~250s | ~10s | 25x |

**Note**: Actual performance depends on your specific hardware.

---

## Alternative: CPU Optimization

If you don't have a GPU, you can still optimize CPU performance:

### 1. Use Intel MKL (Math Kernel Library)

```powershell
pip install mkl mkl-service
```

### 2. Set environment variables

```powershell
$env:OMP_NUM_THREADS="8"  # Set to your CPU core count
$env:MKL_NUM_THREADS="8"
```

### 3. Use half precision (FP16)

Modify `enformer_wrapper.py`:
```python
_MODEL = _MODEL.half()  # Convert to FP16
ref_tensor = ref_tensor.half()
alt_tensor = alt_tensor.half()
```

**Expected speedup**: 1.5-2x on CPU

---

## Next Steps

1. ✅ Install CUDA Toolkit
2. ✅ Reinstall PyTorch with CUDA
3. ✅ Verify with test script
4. ✅ Run CardioVar and enjoy 30x speedup!

For questions or issues, refer to:
- PyTorch CUDA docs: https://pytorch.org/get-started/locally/
- NVIDIA CUDA docs: https://docs.nvidia.com/cuda/
