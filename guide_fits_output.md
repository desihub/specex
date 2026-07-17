# Specex FITS Output Data Model Guide

The output `.fits` files produced by Specex (and our Python port) contain the full spectrograph wavelength solution and the Point Spread Function (PSF) model.

## Extension Overview

### Extension 0: Primary
- **Type**: Image (usually empty or containing a copy of the input).
- **Purpose**: Holds global metadata in the header (MJD, Exposure ID, Camera, etc.).

### Extension 1: XTRACE
- **Type**: Image Array (shape: `[500, 7]`)
- **Purpose**: **CCD X-coordinate vs. Wavelength.**
- **Details**: 500 fibers. Each row contains 7 Legendre coefficients (degree 6). 
- **Usage**: To find the X-pixel position for a given wavelength $\lambda$, evaluate the Legendre polynomial $P(\lambda)$ using these coefficients.

### Extension 2: YTRACE
- **Type**: Image Array (shape: `[500, 7]`)
- **Purpose**: **CCD Y-coordinate vs. Wavelength.**
- **Details**: Same structure as XTRACE.
- **Usage**: Maps wavelength to the Y-pixel (row) on the detector.

### Extension 3: PSF
- **Type**: Binary Table
- **Purpose**: **Gauss-Hermite PSF Shape Parameters.**
- **Columns**: Usually 59 columns.
- **Rows**: Each row corresponds to a specific PSF parameter (e.g., `GHSIGX`, `GHSIGY`, `GH-0-0`, `GH-1-1`, `TAILAMP`, etc.).
- **Data**: Each cell contains an array of Legendre coefficients mapping that parameter's variation along the fiber and wavelength.
  - In production Specex, these are often `[500, 4]` or similar (4 coeffs per fiber).
  - In our port, we use bundle-wide 2D Legendres.

### Extension 4: EXTOFF (external wavelength offset), when present
- **Type**: Binary Table
- **Purpose**: **Inherited calibration metadata, not a PSF-fit output.**
- **Fields**: `WAVE`, `DWAVE`, `DWAVE_ERR` -- an external-reference wavelength-offset table written by `desispec.trace_shifts.write_traces_in_psf()`, a separate downstream pipeline stage (`desi_compute_trace_shifts`) that runs *after* PSF fitting and compares extracted spectra to a reference (e.g. sky) spectrum.
- **Provenance**: `desi_compute_psf`'s bundle merge step (`merge_psf()` in `desispec/scripts/specex.py`) only overwrites the `XTRACE`/`YTRACE`/`PSF` HDUs of the `--input-psf` template; any other extension already present in that template (like `EXTOFF`, or `INTOFF` for internal/fiber-to-fiber offsets) is carried through to the output byte-for-byte, unchanged. Neither the C++ fitter nor `merge_psf()` computes these values -- they're whatever was already baked into the input PSF file from an earlier calibration run.
- **Port implication**: `write_python_psf` (`py/specex/io.py`) reproduces this by copying any extension present in `input_template` other than `XTRACE`/`YTRACE`/`PSF` straight into the output, with no computation. Do not try to recompute `DWAVE`/`EXTOFF` values from the fit -- that's not what C++ does either.
