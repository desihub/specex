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
mask=hdulist["FMASK"].data
valid_fibers=numpy.where(mask==0)[0]
median_flat=numpy.median(spectra[valid_fibers,:], axis=0)
spectra /= median_flat
invar *= median_flat*median_flat

# remove outliers due to cosmics
wstep=50
knots=wave[0]+wstep/2+wstep*numpy.arange(int((wave[-1]-wave[0])/wstep))
wstep=200
loose_knots=wave[0]+wstep/2+wstep*numpy.arange(int((wave[-1]-wave[0])/wstep))
for fiber in valid_fibers :
    spec=spectra[fiber]
    specivar=invar[fiber]
    indices=numpy.where((spec<0)|(spec>1.5))[0]
    specivar[indices]=0
    for loop in range(100) :
        tck = scipy.interpolate.splrep(wave,spec,task=-1,t=knots,w=specivar)
        mod = scipy.interpolate.splev(wave,tck,der=0)
        res = spec-mod
        rms = numpy.std(res[numpy.where(specivar>0)[0]])
        # detect 5 sigma outliers of fit and set ivar=0
        indices=numpy.where((res**2>(5*rms)**2)&(specivar>0))[0]
        print "fiber %d iter %d rms %f outliers %d"%(fiber,loop,rms,len(indices))
        if len(indices)==0 :
            break
        specivar[indices]=0
        
    # replace by interpolation
    tck = scipy.interpolate.splrep(wave,spec,task=-1,t=loose_knots,w=specivar)
    mod = scipy.interpolate.splev(wave,tck,der=0)
    indices=numpy.where(specivar>0)
    meaninvar=numpy.mean(specivar[indices])
    indices=numpy.where(specivar==0)
    if len(indices)>0 :
        spec[indices]=mod[indices]
        specivar[indices]=meaninvar
    
    spectra[fiber]=spec
    invar[fiber]=specivar
    
    
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


