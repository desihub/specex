#!/usr/bin/env python

import pyfits,sys,json,pylab,string,numpy,os
from math import *


if len(sys.argv)<4 :
    print sys.argv[0],"inspec.fits inflat.fits outspec.fits inflat.fits"
    sys.exit(12);

inspecfilename=sys.argv[1]
inflatfilename=sys.argv[2]
outspecfilename=sys.argv[3]

hdulist=pyfits.open(inspecfilename)
spectra=hdulist[0].data
invar=hdulist[1].data
hdulist2=pyfits.open(inflatfilename)
flat=hdulist2[0].data
invarflat=hdulist2[1].data
spectra /= flat
# scale invar and add noise of flat field itself (ignore uncertainty on the mean flat continuum)
# var_1 = var_0/flat/flat + ((sigma(flat)/flat**2)*spectra_0)**2
#       = var_0/flat/flat + ((sigma(flat)/flat)*spectra_1)**2
#       = var_0/flat/flat + var_flat/(flat*flat)*(spectra_1*spectra_1)
#       = 1/(invvar_0*flat*flat) + 1/(invarflat*flat*flat)*(spectra_1*spectra_1)
 
invvar = 1/(1/(invar*flat*flat)+ 1/(invarflat*flat*flat)*spectra*spectra)

if os.path.isfile(outspecfilename) :
    os.unlink(outspecfilename)

hdulist.writeto(outspecfilename)
