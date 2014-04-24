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
mask=hdulist["FMASK"].data
valid_fibers=numpy.where(mask==0)[0]

print "%d valid fibers in %s"%(valid_fibers.shape[0],inspecfilename)

hdulist2=pyfits.open(inflatfilename)
flat=hdulist2[0].data
invarflat=hdulist2[1].data
spectra[valid_fibers] /= flat[valid_fibers]
# scale invar and add noise of flat field itself (ignore uncertainty on the mean flat continuum)
# var_1 = var_0/flat/flat + ((sigma(flat)/flat**2)*spectra_0)**2
#       = var_0/flat/flat + ((sigma(flat)/flat)*spectra_1)**2
#       = var_0/flat/flat + var_flat/(flat*flat)*(spectra_1*spectra_1)
#       = 1/(invvar_0*flat*flat) + 1/(invarflat*flat*flat)*(spectra_1*spectra_1)
 
invar[valid_fibers] = 1/(1/(invar[valid_fibers]*flat[valid_fibers]*flat[valid_fibers])+ 1/(invarflat[valid_fibers]*flat[valid_fibers]*flat[valid_fibers])*spectra[valid_fibers]*spectra[valid_fibers])

if os.path.isfile(outspecfilename) :
    os.unlink(outspecfilename)

hdulist.writeto(outspecfilename)
