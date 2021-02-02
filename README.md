
# specex

This repository contains code for PSF measurement in fiber-fed spectrograph for DESI.

## Introduction

Package for PSF measurement in fiber-fed spectrograph (BOSS data, DESI simulations).
This package is intended to be used with [HARP](https://github.com/tskisner/HARP) and
/or [specter](https://github.com/desihub/specter) extraction codes.
The installation procedure is described in the [INSTALL file](INSTALL.md), as well as below. The code depends on HARP.

## Installation

```

# uncomment below for consistency with current harpconfig compiler (gcc)
# used to compile specex master branch running at nersc
# this should give identical results for the psf
# export CC=gcc
# export CXX=gcc

git clone --single-branch --branch io_refactor https://github.com/desihub/specex
cd specex
mkdir build
cd build
cmake ../
make -j8
export PYTHONPATH=$PWD:$PYTHONPATH

```

## Using specex for DESI

Access to specex in python is through a shard object library and desi_compute_psf
See [desispec.scripts.specex](https://github.com/desihub/desispec/blob/spx_io_refactor/py/desispec/scripts/specex.py) for more details on how to run specex for multiple bundles in serial or parallel and then merge. Here is an example:
```
from desispec.pybindspecex import run_specex

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

