#!/usr/bin/env python

import sys,numpy,pylab
import astropy.io.fits as pyfits
import argparse
from numpy.polynomial.legendre import legval
import numpy as np

def eval(fiber,wavelength) :
    wr=2*(wavelength-(wavemin+wavemax)/2)/(wavemax-wavemin)
    return legval(wr,fiber_coeffs[fiber])

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--psf', type = str, default = None, required=True,
                        help = 'path of DESI PSF fits file')
parser.add_argument('--param', type = str, default = None, required=True,
                        help = 'parameter (GHSIGX,GHSIGY,GH-2-0,...)')
parser.add_argument('--wavelength', type = float, default = None, required=False,
                        help = 'for this wavelength, as a function of fiber number')
parser.add_argument('--fiber', type = float, default = None, required=False,
                        help = 'for this fiber, as a function of wavelength')

args = parser.parse_args()


hdulist=pyfits.open(args.psf)

print hdulist[1].data["PARAM"]
paramindex=numpy.where(hdulist[1].data["PARAM"]==args.param)[0][0]
wavemin=hdulist[1].data["WAVEMIN"][paramindex]
wavemax=hdulist[1].data["WAVEMAX"][paramindex]
fibermin=int(hdulist[1].header["FIBERMIN"])
fibermax=int(hdulist[1].header["FIBERMAX"])
 
print "fibermin,fibermax =",fibermin,fibermax
print "wavemin,wavemax   =",wavemin,wavemax



legdeg=hdulist[1].header["LEGDEG"]
table=hdulist[1].data
fiber_coeffs=table["COEFF"][paramindex]


if args.wavelength is not None :

    if args.wavelength < wavemin or args.wavelength > wavemax :
        print "ERROR wavelength %f outside of valid range [%f,%f]"%(args.wavelength,wavemin,wavemax)
        sys.exit(12)
    
    pylab.figure()
    x=numpy.zeros((fibermax+1-fibermin))
    fibers=numpy.arange(fibermin,fibermax+1)
    for fiber in fibers :
        x[fiber-fibermin]=eval(fiber-fibermin,args.wavelength)
    pylab.plot(fibers,x,"o-")
    pylab.xlabel("fiber #")
    pylab.ylabel("%s @ %dA"%(args.param,int(args.wavelength)))

if args.fiber is not None :

    if args.fiber < fibermin or args.fiber > fibermax :
        print "ERROR fiber %d outside of valid range [%d,%d]"%(args.fiber,fibermin,fibermax)
        sys.exit(12)
    
    pylab.figure()
    wave=np.linspace(wavemin,wavemax,500)
    x=np.zeros(wave.shape)
    for i in range(x.size) :
        x[i]=eval(args.fiber-fibermin,wave[i])
    pylab.plot(wave,x,"-")
    pylab.xlabel("wavelength [A]")
    pylab.ylabel("%s for fiber #%d"%(args.param,args.fiber))

pylab.show()
