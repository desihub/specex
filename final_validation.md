# Final Randomized Validation Report
**Date:** June 12, 2026  
**Environment:** Perlmutter A100 GPU Node (4x A100 GPUs, 128 CPUs)  
**Baseline:** C++ `desi_psf_fit` / `desi_compute_psf` (with lamp line parsing fix)

## 1. Overview
This report documents the final large-scale consistency check between the legacy C++ implementation and the new Python/JAX GPU-accelerated implementation. We randomly sampled **30 bundles** (10 per arm) and **30 full CCDs** (10 per arm) across all arms (B, R, Z) to verify numerical parity and performance gains.

### Test Matrix:
- **Nights:** 20260401
- **Exposure:** 00344649
- **Parallelism:** 4 GPU JAX (Python) vs 20 MPI ranks (C++)
- **Optimization:** Optimized for high-throughput parallel execution.
- **Scope:** 10 samples per arm (bundles), 10 samples per arm (CCDs).

---

## 2. Randomized Bundle Results (30 Samples)
*Fitting 10 random bundles from each arm (B, R, Z).*

| Cam | Bundle | Spots | Time Py (s) | Time C++ (s) | X-Trace RMS (px) | Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **b1** | 0 | 616 | 62.8 | 53.8 | 0.0163 | **PASS** |
| **b1** | 17 | 691 | 66.8 | 57.3 | 0.0217 | **PASS** |
| **b4** | 7 | 606 | 65.5 | 50.7 | 0.0374 | **PASS** |
| **b3** | 4 | 602 | 65.7 | 54.6 | 0.0192 | **PASS** |
| **b1** | 18 | 613 | 63.9 | 49.9 | 0.0188 | **PASS** |
| **b0** | 2 | 611 | 62.8 | 49.5 | 0.0220 | **PASS** |
| **b6** | 1 | 579 | 62.2 | 53.2 | 0.0189 | **PASS** |
| **b3** | 7 | 625 | 56.3 | 49.2 | 0.0150 | **PASS** |
| **b8** | 19 | 533 | 60.4 | 44.7 | 0.0140 | **PASS** |
| **b0** | 17 | 637 | 64.5 | 57.7 | 0.0196 | **PASS** |
| **r3** | 17 | 1339 | 61.5 | 119.9 | 0.0113 | **PASS** |
| **r6** | 7 | 1298 | 61.3 | 130.6 | 0.0183 | **PASS** |
| **r7** | 18 | 1233 | 62.8 | 101.0 | 0.0184 | **PASS** |
| **r4** | 0 | 1294 | 62.0 | 158.9 | 0.0084 | **PASS** |
| **r2** | 13 | 1243 | 60.3 | 105.0 | 0.0276 | **PASS** |
| **r5** | 8 | 1306 | 60.8 | 113.3 | 0.0260 | **PASS** |
| **r2** | 6 | 1307 | 62.5 | 127.6 | 0.0264 | **PASS** |
| **r5** | 3 | 1272 | 60.4 | 135.7 | 0.0134 | **PASS** |
| **r1** | 12 | 1318 | 61.9 | 115.4 | 0.0462 | **PASS** |
| **r1** | 11 | 1325 | 61.1 | 121.2 | 0.0496 | **PASS** |
| **z5** | 19 | 1588 | 63.1 | 106.9 | 0.0092 | **PASS** |
| **z4** | 1 | 1479 | 63.0 | N/A | -- | **FAIL** (C++ Baseline Failed) |
| **z7** | 17 | 1586 | 65.4 | 160.1 | 0.0162 | **PASS** |
| **z1** | 12 | 1566 | 61.6 | 116.1 | 0.0113 | **PASS** |
| **z1** | 17 | 1588 | 57.2 | 127.3 | 0.0169 | **PASS** |
| **z5** | 18 | 1475 | 62.8 | 122.8 | 0.0093 | **PASS** |
| **z4** | 19 | 1562 | 62.7 | 173.6 | 0.0118 | **PASS** |
| **z3** | 2 | 1555 | 65.6 | 143.7 | 0.0101* | **PASS** |
| **z0** | 7 | 1526 | 69.2 | 165.1 | 0.0351 | **PASS** |
| **z4** | 2 | 1493 | 63.7 | 169.3 | 0.0459 | **PASS** |

*\*Note: z3 bundle 2 Trace RMS excludes fiber 65 which C++ zeroed out. When excluding that fiber, parity is excellent (0.0101 px).*
*\*Note: z4 bundle 1 C++ failed with exit status 1. Python implementation succeeded.*

---

## 3. Randomized Full CCD Results (30 Samples)
*Fitting entire CCDs (20 bundles each) using parallel drivers.*

| Cam | Time Py (min) | Time C++ (min) | X-Trace RMS (px) | Speedup | Status |
| :--- | :---: | :---: | :---: | :---: | :--- |
| *In Progress* | ... | ... | ... | ... | ... |

---

## 4. Final Summary
*To be populated upon completion.*
