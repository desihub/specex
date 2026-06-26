# State of the Project: Specex Porting (C++ to Python/JAX)

## Goal
Achieve exact spot selection parity and low centroid RMS between the Python/JAX implementation and the C++ baseline for the NIR (z8) band.

## Current Status
- **Divergence:** Python finds 1533 spots, C++ finds 1523 spots for z8 Bundle 5.
- **C++ Baseline:** The C++ implementation is currently failing to write the final "relaxed" spot selection to disk, resulting in filenames like `cppspots_pass4.txt` having unexpectedly low counts (e.g., 703 or 959) instead of the expected 1523.
- **Hypothesis:** The `select_spots` function is called multiple times in `FitEverything`. Since `selection_pass_count` is static, it increments on every call, but the `output_filename` is only passed in a few of those calls, making the resulting filenames (`passN.txt`) misleading and gaps in the exported data.

## Progress
- [x] Fixed C++ build and formatting.
- [x] Implemented high-precision (15 decimal places) exports for candidates and selected spots.
- [x] Verified that both implementations generate exactly 1700 raw candidates.
- [x] Aligned Python S/N threshold (3.0) and distance check logic with C++.
- [x] Implemented "Checkpoint Exports" (`cpp_cp1...`, `cpp_cp2...`) in both implementations.
- [x] Standardized precision to `.15f` / `std::setprecision(15)`.

## Key Decisions
- **Precision:** Use 15 decimal places for all wavelength and coordinate exports to eliminate matching failures.
- **Diagnostics:** Use indexed checkpoint files to pinpoint where the C++ selection logic diverges from Python.
- **Consistency:** Pass the `output_filename` to every call of `select_spots` in C++ to ensure a complete audit trail of selection passes.

## Next Steps (For New Thread)
1. **Fix C++ Export Logic:** Modify `src/specex_psf_fitter.cc` to pass `psf->output_psf_filename` to every `select_spots` call in `FitEverything`.
2. **Verify Spot Count:** Run the C++ fit and confirm that the final pass export reaches exactly 1523 spots.
3. **Isolate the 10-Spot Difference:** Identify the specific 10 spots that Python includes but C++ excludes (1533 vs 1523).
4. **Centroid Parity:** Calculate the Relative Centroid RMS for the overlapping 1523 spots.

## Critical Context
- **Case:** z8 Bundle 5 (z8-00344649).
- **C++ Fit Command:** `desi_psf_fit -a /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz --in-psf /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits --lamp-lines /global/cfs/cdirs/desi/users/cdwarner/code/specex/py/specex/data/specex_linelist_desi.txt --out-psf /pscratch/sd/c/cdwarner/specex/cppfit-psf-z8-00344649_05.fits --first-bundle 5 --last-bundle 5 --first-fiber 125 --last-fiber 149 --legendre-deg-wave 3 --fit-continuum --broken-fibers 473,474`
- **Python Fit Command:** `python -m specex.specex -a /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz --in-psf /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits --lamp-lines /global/cfs/cdirs/desi/users/cdwarner/code/specex/py/specex/data/specex_linelist_desi.txt --out-psf /pscratch/sd/c/cdwarner/specex/pyfit-psf-z8-00344649_05.fits --first-bundle 5 --last-bundle 5 --first-fiber 125 --last-fiber 149 --legendre-deg-wave 3 --fit-continuum --broken-fibers 473,474 --gpu 4`
