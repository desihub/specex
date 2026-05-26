# Specex C++ Code Summary and Porting Analysis

## Overview
`specex` is a tool for fitting the spectroscopic Point Spread Function (PSF) and wavelength solution for the DESI experiment. It models the PSF as a 2D Gauss-Hermite polynomial where the coefficients themselves vary across the CCD as Legendre polynomials of wavelength and fiber number.

## Code Structure
- **Python Layer (`py/specex/`)**: Provides high-level orchestration, I/O handling (via `fitsio`), and QA. It interfaces with C++ via `pybind11`.
- **C++ Core (`src/`)**:
    - `PSF_Fitter`: Orchestrates the fitting process (Gauss-Newton with Brent line search).
    - `GaussHermitePSF`: Implements the 2D Gauss-Hermite PSF model.
    - `ComputeChi2AB`: The computational heart. It iterates over pixels and spots to fill the Normal Equations matrix (A) and vector (B).
    - `specex_linalg`: Custom and BLAS/LAPACK wrappers for linear algebra.

## Performance Bottlenecks

### 1. Matrix Filling (`ComputeChi2AB`)
This is the most significant bottleneck.
- **Complexity**: $O(N_{iterations} \times N_{pixels} \times N_{spots})$.
- **Details**: For every pixel in a bundle's footprint, the code loops over all spots (emission lines) that might contribute light. It evaluates the PSF and its derivatives with respect to all free parameters (fluxes, trace positions, PSF shapes).
- **Parallelization**: Currently uses OpenMP to split the image into horizontal bands.
- **GPU Opportunity**: This is highly parallelizable. Pixel-wise calculations can be mapped to GPU threads. Accumulating the $A$ matrix ($A = \sum w H H^T$) is a classic reduction problem.

### 2. PSF Evaluation (`PSFValueWithParamsXY`)
- **Details**: Calculates Hermite polynomials and exponentials.
- **GPU Opportunity**: Extremely well-suited for vectorized evaluation in JAX or CuPy.

### 3. Linear System Solving (`cholesky_solve`)
- **Details**: Solves $(A^T W A) \Delta P = A^T W \Delta d$.
- **Scale**: Matrix size is typically $1000 \times 1000$ to $2000 \times 2000$.
- **Status**: Fast on CPU via OpenBLAS/MKL, but can be moved to GPU (e.g., `cupy.linalg.solve`) to avoid CPU-GPU transfers.

## Serial Bottlenecks and Showstoppers

### Serial Logic
- **Brent Line Search**: The outer minimization loop and the Brent step-size optimization are inherently serial but represent a small fraction of the total time.
- **Orchestration**: The logic deciding which parameters to fit in which order (Trace -> PSF -> Flux -> Tail) is serial.

### Showstoppers / Risks
- **Memory Usage**: The $A$ matrix is $N_{par} \times N_{par}$. For $N_{par} \approx 2000$, this is small (~32MB). However, intermediate storage of derivatives for all pixels/spots could be huge. We must use fused kernels or memory-efficient reductions.
- **Complex Parameter Mapping**: The code uses Legendre polynomials to model how PSF parameters vary. Porting this logic exactly is critical for consistency.
- **OpenMP Integration**: The current parallelization is relatively coarse (bands). A naive Python port might be slower if not properly vectorized.

## Porting Plan

### Phase 1: Foundation (Python/NumPy)
- Implement `GaussHermitePSF` and `Legendre` evaluation in pure Python/NumPy.
- Create unit tests comparing Python output against C++ output for specific spots.
- Port I/O and orchestration to minimize reliance on `_libspecex`.

### Phase 2: Vectorization and JAX/CuPy
- Convert `ComputeChi2AB` logic to use vectorized NumPy/CuPy operations.
- Use JAX for automatic differentiation of the PSF model (eliminating manual derivative code).
- Implement the "Matrix Filling" step as a JAX JIT-compiled function.

### Phase 3: GPU Acceleration and MPI
- Optimize the GPU kernels for Perlmutter (A100).
- Integrate MPI to distribute bundles across multiple GPUs/Nodes.
- Target: Fit all 600 bundles of a CCD in the time it currently takes for one.

### Phase 4: Validation and Refinement
- Run full-scale comparisons between C++ and Python outputs.
- Profile and tune memory usage for large footprints.
