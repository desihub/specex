# Specex Porting Log

## 2026-05-26 14:00 (approx)
### Phase 1: Foundation Completed
- **Legendre Polynomials:** Implemented `Legendre1DPol`, `Legendre2DPol`, and `SparseLegendre2DPol` in `py/specex/math.py`.
- **Gauss-Hermite PSF Core:** Implemented pixel value integration in `py/specex/psf.py`.
- **Verification:** Achieved bit-accurate parity with C++ using direct side-by-side comparison via updated `pybind11` bindings.
- **Environment:** Created `env_setup.sh` and `testing/test_math_psf.py`.

## 2026-05-26 16:30 (approx)
### Phase 2: JAX Integration and Initial Acceleration
- **Vectorization:** Refactored `GaussHermitePSF` to use `jax.numpy` and NumPy broadcasting, eliminating triple-nested loops.
- **Performance:** Achieved **~100x speedup on CPU** via JAX JIT (0.0045s warm vs 0.46s NumPy for 1 million pixel-spot pairs).
- **GPU Readiness:** Verified `jax` and `cupy` imports on Perlmutter compute nodes.

## 2026-05-26 18:00 (approx)
### Major Milestone: Automatic Differentiation and Fitter Core
- **Key Accomplishments:**
    1. **Differentiable PSF Core:** Refactored `GaussHermitePSF` to be fully compatible with JAX's AD engine. `single_pix_value_jnp` now serves as the atomic, differentiable unit.
    2. **Vectorized Derivatives:** Used `jax.jacfwd` and `vmap` in `py/specex/fitter.py` to automatically compute Jacobians for spot centers (xc, yc) and PSF shape parameters (gh_params). This **eliminates thousands of lines** of manual C++ derivative code.
    3. **Efficient Matrix Filling:** Implemented `compute_chi2_ab_full` using `jax.einsum` and `jnp.dot` to accumulate Gauss-Newton matrices directly on the accelerator.
    4. **Brent Line Search:** Ported the robust 1D minimization algorithm from `src/specex_brent.cc` to ensure convergence.
- **Current Status:**
    - Flux fitting: Fully implemented and verified.
    - Position/Shape fitting: AD logic implemented and integrated into the Gauss-Newton loop.
    - End-to-end fit: Initial pipeline operational in Python.

## 2026-05-27 10:00 (approx)
### Current Session: Full Bundle Verification
- **Fixed:** `read_lamp_lines` parser to handle `NAME WAVE` format.
- **Fixed:** `load_python_psf` to call `pyio.load_psf`, correctly bridging traces and 1D parameter models from C++.
- **Blocker (Resolved):** `IndexError` in JAX due to unpopulated PSF parameters (now correctly loaded).
- **Blocker (Current):** `ValueError: operands could not be broadcast together`. 
    - **Cause:** The `fit` loop expected a full parameter update (`delta_P`), but I temporarily used a flux-only solver which returned an empty array for non-linear updates.
    - **Status:** The Python fit successfully initialized and completed a flux-only pass with a stable $\chi^2$, proving the basic I/O and projection logic is sound.
- **Plan:**
    1. Re-enable full AD solver with memory-efficient chunking to avoid OOM on large bundles.
    2. Compare Python vs C++ final fitted parameters.
    3. Benchmark GPU performance on A100.
