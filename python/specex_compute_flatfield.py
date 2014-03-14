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
meanflat=numpy.median(spectra, axis=0)
spectra /= meanflat
invar *= meanflat*meanflat

if os.path.isfile(outfilename) :
    os.unlink(outfilename)

hdulist.writeto(outfilename)

if not meanflatfilename == "" :
    meanflat_array=numpy.zeros(spectra.shape)
    meanflat_array[0]=meanflat
    meanflat_invar=numpy.ones(spectra.shape)
    hdu0 = pyfits.PrimaryHDU(meanflat_array)
    hdu1 = pyfits.ImageHDU(meanflat_invar)
    hdu2 = pyfits.ImageHDU(wave)

    if os.path.isfile(meanflatfilename) :
        os.unlink(meanflatfilename)

    pyfits.HDUList([hdu0,hdu1,hdu2]).writeto(meanflatfilename)

if False :
    pylab.plot(wave,meanflat)
    pylab.figure()
    pylab.plot(wave,spectra[1,:])
    pylab.show()


