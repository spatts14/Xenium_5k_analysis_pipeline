#!/usr/bin/env python3
"""Check GPU support for PyTorch and scVI."""

print("=== GPU Support Check ===")

# Check PyTorch
try:
    import torch

    print(f"✓ PyTorch version: {torch.__version__}")
    print(f"✓ CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"✓ CUDA version: {torch.version.cuda}")
        print(f"✓ Number of GPUs: {torch.cuda.device_count()}")
        for i in range(torch.cuda.device_count()):
            print(f"  GPU {i}: {torch.cuda.get_device_name(i)}")
    else:
        print("✗ CUDA not available - PyTorch will use CPU only")
except ImportError:
    print("✗ PyTorch not installed")

print()

# Check scVI
try:
    import scvi

    print(f"✓ scVI version: {scvi.__version__}")

    # Check if scVI can detect GPU
    if torch.cuda.is_available():
        print("✓ scVI should be able to use GPU")
        # Test device selection
        device = "cuda" if torch.cuda.is_available() else "cpu"
        print(f"✓ Default device for scVI: {device}")
    else:
        print("! scVI will use CPU (no GPU detected)")

except ImportError:
    print("✗ scVI not installed")

print()

# Check other GPU-related packages
packages_to_check = ["scvi-tools", "jax", "jaxlib"]
for pkg in packages_to_check:
    try:
        __import__(pkg.replace("-", "_"))
        print(f"✓ {pkg} is installed")
    except ImportError:
        print(f"✗ {pkg} not installed")

print("\n=== Test GPU Tensor Operations ===")
if torch.cuda.is_available():
    try:
        # Test basic GPU operations
        x = torch.randn(3, 3).cuda()
        y = torch.randn(3, 3).cuda()
        z = x + y
        print("✓ Basic GPU tensor operations work")

        # Test memory allocation
        memory_allocated = torch.cuda.memory_allocated() / 1024**2  # MB
        print(f"✓ GPU memory allocated: {memory_allocated:.1f} MB")

    except Exception as e:
        print(f"✗ GPU tensor operations failed: {e}")
else:
    print("! Skipping GPU tests (no GPU available)")
