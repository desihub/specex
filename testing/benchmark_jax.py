import os
import sys
import time
import numpy as np
import jax
import jax.numpy as jnp

# Ensure we can import from py/specex
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))

from specex.psf import PSF

def benchmark_psf_jax():
    print(f"JAX devices: {jax.devices()}")
    
    psf = PSF()
    degree = 6
    n_coeff = (degree + 1)**2 - 1
    
    # Large batch of pixels and spots to see GPU advantage
    Np = 10000
    Ns = 100
    
    params = np.random.rand(Ns, 2 + n_coeff)
    params[:, 0:2] += 0.5
    xc = np.random.rand(Ns) + 100
    yc = np.random.rand(Ns) + 100
    xpix = np.random.rand(Np) + 100
    ypix = np.random.rand(Np) + 100
    
    print(f"Benchmarking PSF evaluation for {Np} pixels and {Ns} spots...")
    
    # 1. NumPy version
    t0 = time.time()
    val_numpy = psf.pix_value(xc, yc, xpix, ypix, params, use_jax=False)
    t1 = time.time()
    print(f"NumPy time: {t1 - t0:.4f}s")
    
    # 2. JAX version (first call includes JIT overhead)
    t0 = time.time()
    val_jax = psf.pix_value(xc, yc, xpix, ypix, params, use_jax=True)
    val_jax.block_until_ready()
    t1 = time.time()
    print(f"JAX time (1st call + JIT): {t1 - t0:.4f}s")
    
    # 3. JAX version (2nd call, warm)
    t0 = time.time()
    val_jax = psf.pix_value(xc, yc, xpix, ypix, params, use_jax=True)
    val_jax.block_until_ready()
    t1 = time.time()
    print(f"JAX time (warm): {t1 - t0:.4f}s")
    
    # Verification
    np.testing.assert_allclose(val_numpy, np.array(val_jax), rtol=1e-5, atol=1e-8)
    print("Verification: JAX results match NumPy!")

if __name__ == "__main__":
    benchmark_psf_jax()
