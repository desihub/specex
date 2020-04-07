# specex Change Log

## 0.6.7 (unreleased)

* No changes yet.

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
