#!/usr/bin/env python

import pyfits,sys,json,pylab,string,numpy,os,scipy.interpolate
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

hdulist.writeto(outfilename,clobber=True)

if not meanflatfilename == "" :
    median_flat_array=numpy.zeros((1,wave.shape[0]))
    median_flat_array[0]=median_flat
    median_flat_invar=numpy.zeros(median_flat_array.shape)
    hdu0 = pyfits.PrimaryHDU(median_flat_array)
    hdu1 = pyfits.ImageHDU(median_flat_invar)
    hdu2 = pyfits.ImageHDU(wave)

    pyfits.HDUList([hdu0,hdu1,hdu2]).writeto(meanflatfilename,clobber=True)

if False :
    pylab.plot(wave,median_flat)
    pylab.figure()
    pylab.plot(wave,spectra[1,:])
    pylab.show()


