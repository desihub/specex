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



## 2026-05-27 15:30 (approx)
- Phase 2: 2D Legendre Parameter Fitting (Completed)
    - Key Implementation Details:
        - **Bundle-Wide Optimization:** Refactored the fitter to optimize shared 2D Legendre coefficients instead of per-spot parameters. This reduces the parameter space and exactly matches the C++ physics model.
        - **Differentiable Mapping:** Implemented a JAX-based mapping that propagates bundle coefficients to local PSF parameters for each spot, enabling full automatic differentiation of the bundle signal.
        - **Matrix Accumulation:** Used `jax.jacfwd` and `jnp.dot` to efficiently fill the Normal Equations matrix for shared parameters.
        - **Refined Footprint:** Updated `get_bundle_footprint` to precisely match the pixel set used by C++, ensuring comparable Chi2 values.
- Current Status:
    - End-to-end Python/JAX fit for Flux, Position, and Shape using bundle-wide shared parameters is complete.
    - Ready for final numerical parity check against C++ baseline.
