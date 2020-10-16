
# specex

This repository contains code for PSF measurement in fiber-fed spectrograph (BOSS data, DESI simulations).

## pip installation
```
git clone --single-branch --branch io_refactor https://github.com/marcelo-alvarez/specex
pip -v install ./specex/specex_pybind
```

## cmake build
```
git clone --single-branch --branch io_refactor https://github.com/marcelo-alvarez/specex
cd specex/specex_pybind
mkdir build
cd build
cmake ../
make -j8
export PYTHONPATH=$PWD:$PYTHONPATH
```
