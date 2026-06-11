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

### Extension 4: Wavelength Residuals
- **Type**: Binary Table
- **Purpose**: **Spot-by-Spot Calibration Diagnostics.**
- **Fields**:
  - `WAVE`: The reference vacuum wavelength (Angstrom) of the arc lamp emission line used for fitting.
  - `DWAVE`: The **Wavelength Residual** ($\Delta\lambda = \lambda_{measured} - \lambda_{model}$). It represents the sub-pixel shift required to perfectly center the PSF on the measured spot relative to the smooth Legendre solution.
  - `DWAVE_ERR`: The 1-sigma statistical uncertainty of the `DWAVE` measurement.

## Interpreting the Corrections
If `DWAVE` is consistently positive or negative in a specific region of the CCD, it indicates a "local" perturbation in the wavelength solution that the smooth degree-6 Legendre polynomial could not capture. In the DESI pipeline, these residuals are often used for high-precision checks of the spectrograph stability.
