# specex : DESI PSF fit using HARP

## Introduction

Package for PSF measurement in fiber-fed spectrograph (BOSS data, DESI simulations).
This package is intended to be used with [HARP](https://github.com/tskisner/HARP) and/or [specter](https://github.com/desihub/specter) extraction codes.
The installation procedure is described in the [INSTALL file](INSTALL.md), the code depends on HARP.

## Using specex for DESI

### specex usage for sparse fiber slit (spectrograph tests)

The code needs an input PSF file which can contain only trace coordinates.

Full fit :
```
desi_fit_psf --arc preproc-arc-lamp-image.fits --in-psf input-psf.fits --out-psf output-psf.fits 
```

Only trace coordinates fit :
```
desi_fit_psf --arc preproc-arc-lamp-image.fits --in-psf input-psf.fits --out-psf output-psf.fits --only-trace-fit
```

Averaging a set of PSF :
```
specex_mean_psf.py -i psf-1.fits psf-2.fits ... -o psf-average.fits
```

Many more options are available, see ```desi_fit_psf --help```

### Parallel computing

specex is not MPI, but it uses internally `openmp` to speed up the computation (note also specex runs significantly faster with the configuration option `--optimize=3`). The number of cpu used is given by the value of the environment variable `OMP_NUM_THREADS`. A convenient way to do parallel computing is to fit the PSF per block of fibers (or bundle, as called in the code), and then merge the results.

Example (in `bash`):

THIS EXAMPLE IS DEPRECATED, THIS HAS TO BE MODIFIED.

```
for BUNDLE in `seq 0 19`; do
  desi_psf_fit ... --first_bundle $BUNDLE --last_bundle $BUNDLE --out-psf-xml psf-of-bundle-${BUNDLE}.xml &
done
wait
specex_merge_psf psf-of-bundle-*.xml --out-fits psf.fits --out-xml  psf.xml
```

## License

License information needed.
