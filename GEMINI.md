# Project: specex

## Description
- This is part of the overall DESI data pipeline

- This code fits the spectcograph wavelength solution and point spread function using arc-lamp spectra with emission lines at known wavelengths.

- It models the PSF as 2D Gauss-Hermite polynomials

- The spot PSFs overlap between neighboring fibers, so the solution is fit simultaneously in bundles of 25 fibers.  Between the bundles there is an extra gap on the CCD which allows us to fit the bundles independently.

## The Project
- We would like to port this to python, with GPU acceleration where possible with a goal of the same speed or hopefully faster and the same output.

- Currently the code is OpenMP parallelized in C++ with outer wrappers in python to distribute the work; we expect that we can achieve similar performance with Python+MPI+numpy+numba+cupy or jax.

- Each spectrograph CCD has 500 fibers = 20 bundles/CCD of 25 fibers/bundle.  Total 30 CCDs * 20 bundles/CCD = 600 bundles.  Currently these are fit using jobs with 3 CPU nodes each; ideally we'll get this to 1 GPU node.

- We are using Perlmutter at Berkeley national lab with A100 GPUs.

- The algorithm is described in section 4.3 of /global/cfs/cdirs/desi/users/cdwarner/code/specex/2209.14482v2.pdf


## Examples:
# fit a single bundle of 25 fibers (~10 seconds on a login node)
```
desi_psf_fit \
  -a /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz \
  --in-psf /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits \
  --lamp-lines $SPECEX/py/specex/data/specex_linelist_desi.txt \
  --out-psf $SCRATCH/fit-psf-z8-00344649_05.fits \
  --first-bundle 5 --last-bundle 5 --first-fiber 125 --last-fiber 149 --legendre-deg-wave 3 --fit-continuum --broken-fibers 473,474
```

# Fit a full CCD PSF (20 bundles) (~3 minutes on a CPU batch node)
```
srun -n 20 desi_compute_psf --mpi --broken-fibers 473,474 \
  --input-image /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz \
  --input-psf /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits \
  --output-psf $SCRATCH/fit-psf-z8-00344649.fits
```

## General Instructions
- Ensure all new functions and classes have appropriate comments

- We should do a git commit on the python branch after every major code change with a meaningful commit message summarizing the changes

- Coding style should avoid multiple commands on a line for better readability.  We should also avoid lines like x, y, z = True, False, True when defining variables for readability.

- We should use best practices for python coding

- We can use mpi, cupy, and jax

- We are open to rewrites that produce similar if not identical results if a literal translation won't be very performance effective.

- As we develop, please document changes by writing summary notes to porting-notes.md so we have a record of the progress and changes.  We should never delete from porting-notes.md - it should be a running log with only additions.

- We should ensure a clear separation of I/O and methods that "do things" in the code and not port over legacy intermediate data products that we don't want to rely on having to be there for an expandable architecture that can flexibly be added to a larger data pipeline for instance.
