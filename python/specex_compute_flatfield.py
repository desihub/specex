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

# we do a spline fit here because we know the relative transmission
# from one fiber to another is a smooth function of wavelength
step1=10 # A, this is the chosen spline resolution
knots1=wave[0]+step1/2+step1*numpy.arange(int((wave[-1]-wave[0])/step1))


# we "expect" the transmission to be much smoother than observed, must be a issue with
# spectro-perf or psf
# compute spline on a much smoother grid to evaluate rms and assign error
step2=200 # A, this is the chosen spline resolution
knots2=wave[0]+step2/2+step2*numpy.arange(int((wave[-1]-wave[0])/step2))


for fiber in range(spectra.shape[0]) :
    tck1=scipy.interpolate.splrep(wave,spectra[fiber],w=invar[fiber],task=-1,t=knots1,k=1) # linear
    spectra[fiber]=scipy.interpolate.splev(wave,tck1,der=0)
    tck2=scipy.interpolate.splrep(wave,spectra[fiber],w=invar[fiber],task=-1,t=knots2,k=3) # 3rd order
    smooth=scipy.interpolate.splev(wave,tck2,der=0)
    chi2pdf=numpy.sum(invar[fiber]*(spectra[fiber]-smooth)**2)/(wave.shape[0]-knots2.shape[0]+2)
    invar[fiber]/=chi2pdf
    print "fiber",fiber,"mean flat error=",1./sqrt(numpy.mean(invar[fiber]))
    # impossible to compute easily the variance !!!

hdulist.writeto(outfilename,clobber=True)

if not meanflatfilename == "" :
    median_flat_array=numpy.zeros(spectra.shape)
    median_flat_array[0]=median_flat
    median_flat_invar=numpy.ones(spectra.shape)
    hdu0 = pyfits.PrimaryHDU(median_flat_array)
    hdu1 = pyfits.ImageHDU(median_flat_invar)
    hdu2 = pyfits.ImageHDU(wave)

    pyfits.HDUList([hdu0,hdu1,hdu2]).writeto(meanflatfilename,clobber=True)

if False :
    pylab.plot(wave,median_flat)
    pylab.figure()
    pylab.plot(wave,spectra[1,:])
    pylab.show()


