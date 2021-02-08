
# specex

This repository contains code for PSF measurement in fiber-fed spectrograph for DESI.

## Introduction

Package for PSF measurement in fiber-fed spectrograph (BOSS data, DESI simulations).
This package is intended to be used with [HARP](https://github.com/tskisner/HARP) and
/or [specter](https://github.com/desihub/specter) extraction codes.
The installation procedure is described in the [INSTALL file](INSTALL.md), as well as below. The code depends on HARP and uses pybind11 (2.2.0)

## Installation

```
# uncomment below for consistency with current harpconfig compiler (gcc)
# used to compile specex master branch running at nersc
# this should give identical results for the psf
# export CC=gcc
# export CXX=gcc

git clone --single-branch --branch io_refactor https://github.com/desihub/specex
cd specex
python setup.py	install	--prefix .

```

## Using specex for DESI

Access to specex in python is through a shared object library extension, `specex._libspecex`:
```
import specex
from specex._libspecex import run_specex

com = ['desi_psf_fit']
com.extend(['-a',
            '/global/cfs/cdirs/desi/spectro/redux/blanc/preproc/20201216/00068217/preproc-b1-00068217.fits'])
com.extend(['--in-psf', '/global/cfs/cdirs/desi/spectro/redux/blanc/exposures/20201216/00068217/shifted-input-psf-b1-00068217.fits'])
com.extend(['--out-psf', './fit-psf-b1-00068217-00.fits'])
com.extend(['--first-bundle', '0'])
com.extend(['--last-bundle', '0'])
com.extend(['--first-fiber', '0'])
com.extend(['--last-fiber', '24'])
com.extend(['--legendre-deg-wave', '1'])

retval = run_specex(com)
```
This should produce a file `fit-psf-b1-00068217-00.fits` in the same directory.

