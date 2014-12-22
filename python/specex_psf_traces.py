#!/usr/bin/env python

import pyfits,sys,numpy
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

if len(sys.argv)<2 :
    print sys.argv[0],"psf.fits"
    sys.exit(0)

hdulist=pyfits.open(sys.argv[1])
#print hdulist[1].header.tostring
fiber=0
for wave in numpy.arange(3600,4500,50) :
    x=eval("X",fiber,wave,verbose=(wave==3600))
    y=eval("Y",fiber,wave)
    print fiber,wave,x,y
