
# specex

This repository contains code for PSF measurement in fiber-fed spectrograph for DESI.

## Installation
```

# uncomment below for consistency with current harpconfig compiler (gcc)
# used to compile specex master branch running at nersc
# this should give identical results for the psf
# export CC=gcc
# export CXX=gcc

git clone --single-branch --branch io_refactor https://github.com/desihub/specex
cd specex
python setup.py install --prefix .

```
