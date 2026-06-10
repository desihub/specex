# Specex Python/GPU Usage Guide

This guide documents how to run the ported Python/JAX version of Specex and the various validation scripts used to verify numerical parity with the C++ baseline.

## 1. Environment Setup

Before running any scripts, ensure your environment is set up correctly on Perlmutter.

```bash
# From the project root
source env_setup.sh
```

This sets the `PYTHONPATH` to include the local `py` directory and loads necessary modules like `cudatoolkit`.

---

## 2. Finding Test Data

If you want to run on a specific night or exposure, use `select_test_case.py` to find the correct file paths and parameters (like `--broken-fibers`).

```bash
# List all cases for a specific night
python testing/select_test_case.py --night 20260401 --list

# Select a random case for testing
python testing/select_test_case.py --night 20260401 --random
```

---

## 3. Running Parity Comparisons (Recommended)

The primary tool for verifying Python vs C++ metrics is `testing/instrumentation_analysis.py`. It runs both the production C++ code and the new JAX-GPU code on a single bundle and compares results.

### Basic Usage:
```bash
python testing/instrumentation_analysis.py --night 20260401 --expid 00344649 --cameras b0,r3,z8 --bundle 5
```

### Metrics Produced:
The script generates a table (default output: `instrumentation_analysis.txt`) containing:
*   **Time(s):** Wall-clock time for the fit.
*   **Chi2:** Final chi-squared value (lower is generally better).
*   **Spots:** Number of spots identified and used in the fit.
*   **XT RMS / YT RMS:** Root-Mean-Square difference in pixels between the Python-fitted traces and the C++-fitted traces.

---

## 4. Multi-Mode Validation

To compare **CPP** vs **JAX-CPU** vs **JAX-GPU** simultaneously, use `testing/validate_all_modes.py`.

```bash
python testing/validate_all_modes.py --cameras b0 --bundle 5 --output comparison_results.txt
```

This is useful for verifying that the GPU acceleration doesn't introduce numerical divergence from the CPU version of the same code.

---

## 5. Running the Full Python CCD Fit

You can run a full CCD fit directly from the command line using the `specex.specex` module.

```bash
python -m specex.specex \
    -a path/to/preproc.fits \
    --in-psf path/to/input-psf.fits \
    --out-psf output-psf.fits \
    --broken-fibers 473,474 \
    --gpu 4
```

### CLI Arguments:
*   `-a`, `--arc`: Input preprocessed arc image.
*   `--in-psf`: Input PSF file (the "shifted" version).
*   `--out-psf`: Path where the fitted PSF will be saved.
*   `--first-bundle` / `--last-bundle`: Range of bundles to fit (0-19).
*   `--gpu`: Number of GPUs to utilize (default: 4).
*   `--broken-fibers`: List of fiber IDs to exclude from the fit.
*   `--sn-threshold`: Signal-to-Noise threshold for spot selection (default: 3.0).
*   `--h-size-y`: Override the PSF stamp half-size in Y (default: 5).

Alternatively, you can call it from within another Python script:

```python
from specex.specex import fit_ccd_native

fit_ccd_native(
    arc_file='path/to/preproc.fits',
    in_psf_file='path/to/input-psf.fits',
    out_psf_file='output-psf.fits',
    lamp_lines_file='py/specex/data/specex_linelist_desi.txt',
    broken_fibers="473,474"
)
```

---

## 6. Understanding the Metrics

*   **Spots:** If the Python version identifies significantly fewer spots than C++, check the `--sn` (S/N threshold) parameter in the scripts.
*   **Chi2:** A 2-4% difference is currently expected due to "Dead Column Masking" differences in the pre-processor.
*   **Trace Deltas (XT/YT RMS):** We target values < 0.05 pixels. Current results typically show ~0.02 pixels.
*   **Chi2: -1.0:** Usually indicates a crash or a failure to find any spots (often due to the PSF loading bug fixed on June 9th).
