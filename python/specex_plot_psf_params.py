#!/usr/bin/env python

import sys,numpy,pylab
import astropy.io.fits as pyfits

from numpy.polynomial.legendre import legval

def eval(pname,fiber,wavelength,verbose=False) :
    legdeg=hdulist[1].header["LEGDEG"]
    table=hdulist[1].data
    paramindex=numpy.where(table["PARAM"]==pname)[0][0]
    wavemin=table["WAVEMIN"][paramindex]
    wavemax=table["WAVEMAX"][paramindex]
    fiber_coeffs=table["COEFF"][paramindex][fiber]
    
    
    wr=2*(wavelength-(wavemin+wavemax)/2)/(wavemax-wavemin)
    val=legval(wr,fiber_coeffs)
    if verbose :
        print pname,"fiber=",fiber,"wave=",wavelength,"wrange=",wavemin,wavemax,"coefs=",fiber_coeffs,"val=",val
    return val

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

#fibers=numpy.arange(0,20,1)
fibers=numpy.arange(fibermin,fibermax+1,10)
for fiber in fibers :
    print fiber
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
pylab.plot(fibers,x)
pylab.plot(fibers,x,"o")
pylab.xlabel("fiber #")
pylab.ylabel(param)

pylab.show()
