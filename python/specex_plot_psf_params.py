#!/usr/bin/env python

import sys,numpy,pylab
import astropy.io.fits as pyfits

from numpy.polynomial.legendre import legval

def eval(pname,fiber,wavelength) :
    wr=2*(wavelength-(wavemin+wavemax)/2)/(wavemax-wavemin)
    return legval(wr,fiber_coeffs[fiber])

if len(sys.argv)<3 :
    print sys.argv[0],"psf.fits param"
    sys.exit(0)


hdulist=pyfits.open(sys.argv[1])
param=sys.argv[2]

print hdulist[1].data["PARAM"]
paramindex=numpy.where(hdulist[1].data["PARAM"]==param)[0][0]
wavemin=hdulist[1].data["WAVEMIN"][paramindex]
wavemax=hdulist[1].data["WAVEMAX"][paramindex]
wave=numpy.arange(wavemin,wavemax,10)
fibermin=int(hdulist[1].header["FIBERMIN"])
fibermax=int(hdulist[1].header["FIBERMAX"])
 
legdeg=hdulist[1].header["LEGDEG"]
table=hdulist[1].data
fiber_coeffs=table["COEFF"][paramindex]

if 0 :
    fibers=numpy.arange(fibermin,fibermax+1,10)
    for fiber in fibers :
        x=eval(param,fiber-fibermin,wave)
        pylab.plot(wave,x)

    pylab.xlabel("wavelength")
    pylab.ylabel(param)

pylab.figure()
wave=(wavemin+wavemax)/2.
x=numpy.zeros((fibermax+1-fibermin))
fibers=numpy.arange(fibermin,fibermax+1)
for fiber in fibers :
    x[fiber-fibermin]=eval(param,fiber-fibermin,wave)
pylab.plot(fibers,x,"o-")
pylab.xlabel("fiber #")
pylab.ylabel(param)

pylab.show()
