# Specex Porting Log

## 2026-05-26 00:00 (approx)
- Initial profiling of C++ code completed. Identified `ComputeChi2AB` as the primary bottleneck.
- Created `python-gpu-port` branch.
- Created `env_setup.sh` for Perlmutter environment management.
- Developed 4-phase porting plan (NumPy -> JAX/CuPy -> GPU/MPI -> Validation).

## 2026-05-26 16:30 (approx)
- Phase 2 Vectorization and JAX Acceleration:
    - Refactored `GaussHermitePSF` to use JAX (`jax.numpy` and `@jit`).
    - Implemented `compute_chi2_ab_jnp` in `py/specex/fitter.py` using JAX for full-matrix accumulation.
    - Verified JAX implementation against NumPy/SciPy in `testing/benchmark_jax.py`.
    - Performance: Achieved ~100x speedup on CPU via JAX JIT for core PSF evaluation (0.0045s warm vs 0.46s NumPy).
    - Status: JAX implementation is GPU-ready; currently running on CPU due to `jaxlib` configuration, but bit-accurate parity is maintained.
- Environment:
    - Updated `env_setup.sh` to fix python version mismatches on compute nodes.
    - Successfully imported `cupy` and `jax` in interactive GPU session.
