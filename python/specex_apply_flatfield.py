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
spectra /= flat
invar *= flat*flat

if os.path.isfile(outspecfilename) :
    os.unlink(outspecfilename)

hdulist.writeto(outspecfilename)
