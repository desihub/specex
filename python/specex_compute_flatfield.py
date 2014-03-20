#!/usr/bin/env python

import pyfits,sys,json,pylab,string,numpy,os
from math import *


if len(sys.argv)<3 :
    print sys.argv[0],"inspec.fits outspec.fits (meanflat.fits)"
    sys.exit(12);

infilename=sys.argv[1]
outfilename=sys.argv[2]
meanflatfilename=""
if len(sys.argv)>3 :
    meanflatfilename=sys.argv[3]

hdulist=pyfits.open(infilename)
spectra=hdulist[0].data
invar=hdulist[1].data
wave   =hdulist[2].data
median_flat=numpy.median(spectra, axis=0)
spectra /= median_flat
invar *= median_flat*median_flat

if os.path.isfile(outfilename) :
    os.unlink(outfilename)

hdulist.writeto(outfilename)

if not median_flatfilename == "" :
    median_flat_array=numpy.zeros(spectra.shape)
    median_flat_array[0]=median_flat
    median_flat_invar=numpy.ones(spectra.shape)
    hdu0 = pyfits.PrimaryHDU(median_flat_array)
    hdu1 = pyfits.ImageHDU(median_flat_invar)
    hdu2 = pyfits.ImageHDU(wave)

    if os.path.isfile(median_flatfilename) :
        os.unlink(median_flatfilename)

    pyfits.HDUList([hdu0,hdu1,hdu2]).writeto(median_flatfilename)

if False :
    pylab.plot(wave,median_flat)
    pylab.figure()
    pylab.plot(wave,spectra[1,:])
    pylab.show()


