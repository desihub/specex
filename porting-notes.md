# Specex Porting Log

## 2026-05-26 00:00 (approx)
- Initial profiling of C++ code completed. Identified `ComputeChi2AB` as the primary bottleneck.
- Created `python-gpu-port` branch.
- Created `env_setup.sh` for Perlmutter environment management.
- Developed 4-phase porting plan (NumPy -> JAX/CuPy -> GPU/MPI -> Validation).

## 2026-05-26 14:00 (approx)
- Phase 1 Foundation Completed:
    - Implemented `Legendre1DPol`, `Legendre2DPol`, and `SparseLegendre2DPol` in `py/specex/math.py`.
    - Implemented `GaussHermitePSF` pixel value integration in `py/specex/psf.py`.
    - Verified Python math implementations against C++ using direct side-by-side comparison.
    - Modified `src/_libspecex.cpp` to expose low-level classes for verification (renamed to `Cpp*` to avoid collisions).
    - Implemented high-level `PSF` and `PSF_Params` classes to manage parameter variation across the CCD.
- Environment:
    - Created `env_setup.sh` to manage modules and `PYTHONPATH`.
    - Created `testing/test_math_psf.py` for automated verification.
