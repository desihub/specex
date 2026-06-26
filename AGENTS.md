# AGENTS.md - Specex Porting (C++ to Python/JAX)

## High-Level Architecture
- **Core Logic:** Ported from C++ (`src/`) to Python/JAX (`py/specex/`).
- **Computational Engine:** Uses JAX for GPU acceleration on A100s. Implements analytical Jacobians for the Gauss-Hermite PSF.
- **Execution Flow:** 
    - `py/specex/specex.py`: Main CLI and driver.
    - `py/specex/fitter.py`: Numerical optimization and JAX-based fitting.
    - `py/specex/io.py`: FITS I/O and data parsing.
    - `py/specex/math.py`: Legendre polynomials and numerical utilities.

## Critical Operational Context
- **Environment:** Must run `source env_setup.sh` in every session to set `PYTHONPATH` and load CUDA modules.
- **Data Storage:** All output data MUST be saved in `/pscratch/sd/c/cdwarner/specex` or its subfolders.
- **Hybrid State:** The repository maintains both the original C++ core (via `pybind11`) and the new Python/JAX implementation for side-by-side parity verification.
- **C++ Wrapper:** `specex.specex.run_specex` is the entry point for the legacy C++ logic.

## Essential Commands

### Setup & Execution
- **Env Setup:** `source env_setup.sh`
- **Build C++ Core:** `python setup.py build_ext --inplace`
- **Python/JAX CCD Fit:**
  ```bash
  python -m specex.specex -a <arc_fits> --in-psf <in_psf> --out-psf <out_psf> --gpu 4
  ```

### Validation & Parity
- **C++ vs JAX-GPU (Single Bundle):**
  ```bash
  python testing/instrumentation_analysis.py --night <date> --expid <id> --cameras b0,r3,z8 --bundle 5
  ```
- **Multi-Mode (CPP vs CPU vs GPU):**
  ```bash
  python testing/validate_all_modes.py --cameras b0 --bundle 5 --output results.txt
  ```
- **Find Test Data:**
  ```bash
  python testing/select_test_case.py --night <date> --random
  ```

## Development Conventions
- **Porting Log:** Append all major changes, milestones, and numerical results to `porting-notes.md`. Never delete from this file.
- **Instruction Source:** Refer to `GEMINI.md` for the original porting requirements and project goals.
- **Numerical Targets:** 
    - Trace RMS should be $< 0.02$ pixels.
    - $\sim 4\%$ Chi2 difference from C++ is expected due to detector masking differences.
- **GPU Memory:** Use `lax.scan` or surgical AD (local stamps) to avoid OOM on A100s when processing large footprints.
