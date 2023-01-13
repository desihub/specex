# specex Change Log

## 0.8.6 (unreleased)

* No changes yet

## 0.8.5 (2023-01-12)

* I/O updates (PR #65, #67)
* pybind11 STL container fix (PR #64)
* MKL BLAS update (PR #62)

## 0.8.4 (2022-01-04)

* Fatal error scope and prefix for log messages (PR #61)

## 0.8.3 (2021-11-19)

* Don't write merged PSF if SPECEX_ERROR is called (PR #57)

## 0.8.2 (2021-11-10)

* Use OpenBLAS for BLAS by default, instead of MKL (PR #55)

## 0.8.1 (2021-10-07)

* Use input PSF if its properties match those expected (PR #52)
* Fix OpenMP bug (PR #51)
* Fix silent PSF fit error bug (PR #54)

## 0.8.0 (2021-07-06)

* Major refactor to remove boost and harp dependencies.

## 0.7.0 (2021-02-10)

* Major refactor: move all I/O to python and change build method (PR #36)

## 0.6.9 (2021-02-10)

* Comment out blended and faint lines; add missing lines (PR #37)

## 0.6.8 (2020-12-11)

* Fixed floating point exception crash (PR #33)
* Adjusted lines based upon inspection of arc lamp data (0e48b41 and 809c554)
* Added timeout option (PR #34)

## 0.6.7 (2020-04-15)

* Added Xe line.

## 0.6.6 (2020-04-07)

* Don't raise an error when spot is outside of image frame.

## 0.6.5 (2019-12-20)

* Add broken fibers option.
* Do not choose a reference fiber that is off.

## 0.6.4 (2019-10-31)

* Do not fail on binary image HDUs.
* Set SPECEXDATA in module file.

## 0.6.3 (2019-05-30)

* Bring back xml I/O functionalities (were broken)

## 0.6.2 (2018-09-26)

* A bug fix (when a spot is masked)

## 0.6.1 (2018-03-29)

* Versioning utilities

## 0.6.0 (2017-11-10)

* Add a prior on high degree of traces in same bundle + several minor improvements in the trace fit

## 0.5.0 (2017-09-30)

* Major refactor of fitting and formats to be faster and more robust by
  accepting any previous PSF as a starting point for the optimization (PR #21).

## 0.4.5 (2017-07-19)

* Change build of libraries to differentiate between
  internal libraries and loadable modules.  Build
  proper bundles on OS X.

## 0.4.4 (2017-06-21)

* On OS X, disable harp plugin.
* Add portable fenv.h wrapper to fix compile on OS X

## 0.4.3 (2017-03-02)

* Install data files along with code.

## 0.4.2 (2016-11-09)

* fix etc/specex.module to work with desiInstall

## 0.4.1 (2016-11-09)

* Support desiInstall.
* Minor changes to documentation files: README.md, INSTALL.md, etc.
* Better handling of missing variables in the top-level Makefile.

## 0.4.0 (2016-07-01)

* Adds ctypes wrapper functionality

## Some earlier versions

* Uses output of desi_bootcalib.py
* Option to increase trace Legendre polynomial degree
* Other changes we didn't document here...

## 0.3.4 (2015-01-12)

As used for the DESI data challenge 2. See [DESI wiki page (restricted access)](https://desi.lbl.gov/trac/wiki/Pipeline/DataChallenges/2015-01).
