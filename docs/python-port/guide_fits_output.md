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

## Practical: Converting Between (fiber, wavelength) and (X, Y) Pixel Position

This is the recipe both `bundle_parity_suite.py` and the new `testing/wave_xy_convert.py` use. A ready-to-run implementation lives in `testing/wave_xy_convert.py` (`PSFTrace` class + CLI) -- this section explains what it's doing and why, so you can extend or debug it independently (e.g. for the shared ~0.5A wavelength-vs-line-list offset investigation).

**Where the numbers actually live:**
- `XTRACE` and `YTRACE` are each a plain `(n_fibers, ncoeff)` float image array. Row `fib` is the Legendre coefficient vector for fiber `fib` -- **fiber indexing is absolute** (row 130 is fiber 130, CCD-wide 0-499), not bundle-relative. This exact absolute-vs-relative distinction was the root cause of bug #5 (`gh_params` indexing the wrong fiber's PSF shape) -- worth remembering any time you touch per-fiber arrays in this codebase.
- `ncoeff = ` the trace's Legendre degree `+ 1`. Don't hardcode this (the guide above says "7" because that's what our standard z8/00344649 test file happens to use, degree 6) -- read it off the array shape (`xtrace.shape[1] - 1`) since it can vary by input PSF file. **This is a completely different, unrelated degree from the small `wdeg` (1 or 3, band-dependent) used for the within-bundle joint-fit trace/PSF-shape *correction* discussed at length in `porting-notes.md` (task 18/19)** -- that correction is added on top of this per-fiber trace curve during fitting, but the *stored, final* `XTRACE`/`YTRACE` in the output file is always this single flat per-fiber Legendre-in-wavelength curve, regardless of what degree was used internally to fit it. Don't confuse the two when reading code or debugging.
- Both HDU headers carry `WAVEMIN`/`WAVEMAX` (should be identical between XTRACE and YTRACE in a well-formed file -- `PSFTrace.__init__` asserts this). This is the domain the Legendre basis is normalized over.

**The normalization convention** (standard Legendre polynomials are only defined/orthogonal on `[-1, 1]`, so wavelength gets rescaled first):
```
rx(wave) = 2 * (wave - WAVEMIN) / (WAVEMAX - WAVEMIN) - 1
x(wave)  = sum_i  XTRACE[fiber, i] * P_i(rx(wave))
y(wave)  = sum_i  YTRACE[fiber, i] * P_i(rx(wave))
```
where `P_i` is the standard degree-`i` Legendre polynomial (`specex.math.legendre_pol`). This is exactly `Legendre1DPol.value()` in `py/specex/math.py`.

**Forward direction (fiber, wavelength) -> (x, y):** direct evaluation, as above. Both X and Y are well-defined and cheap for any wavelength.

**Inverse direction (fiber, x, y) -> wavelength:** **only Y is usable for this.** Y is the dispersion axis and (by construction, physically) monotonic in wavelength across the fiber's range, so inverting `Y_vs_W` at a measured `y` gives a unique, well-posed wavelength. X is the cross-dispersion axis -- it varies only slowly and non-monotonically with wavelength (see the whole z0:5/degree-1 X-trace discussion in porting-notes.md), so "invert X to get wavelength" is not a meaningful operation and isn't provided. In practice you don't need X for the inversion at all; it's supplied to the CLI/description for symmetry but the code only ever uses `y`. Inversion itself (`Legendre1DPol.invert()`) is a robust fine-grid lookup (1000-point `linspace` over `[WAVEMIN, WAVEMAX]`, then `np.interp`) rather than an analytic polynomial inverse -- fine for anything below ~1e-3 pixel precision, which is well beyond what matters here.

**Usage:**
```bash
# fiber, wave(s) -> x, y
python testing/wave_xy_convert.py fiber_wave_to_xy  fit-psf-z8-00344649.fits 130 8670.33

# fiber, y (measured centroid) -> wave
python testing/wave_xy_convert.py xy_to_wave         fit-psf-z8-00344649.fits 130 2015.46
```
or from Python:
```python
from wave_xy_convert import PSFTrace
t = PSFTrace('fit-psf-z8-00344649.fits')
x, y = t.wave_to_xy(fiber=130, wave=8670.33)   # -> (1128.5726, 2013.1869)
wave  = t.y_to_wave(fiber=130, y=2015.46)       # -> 8671.7137 (round-trip residual is real fit scatter, not a bug)
```
Both methods accept a scalar or a `np.array` for the wavelength/y argument. To compare C++ vs Python (or either vs the line list) at a specific spot, build a `PSFTrace` on each output file and evaluate at the same `(fiber, wave)` or `(fiber, y)` -- this is exactly the `wave_residual_stats()` methodology already used throughout the b/r and z-band campaigns in `porting-notes.md`.
