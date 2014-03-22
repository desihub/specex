#!/usr/bin/env python

import pyfits,sys
from math import *
from specex_air_to_vacuum import *


if len(sys.argv)<3 :
    print sys.argv[0],"inspec.fits outspec.fits"
    sys.exit(12);

infilename    = sys.argv[1]
outfilename   = sys.argv[2]

hdulist       = pyfits.open(infilename)
wave          = hdulist[2].data # in air
wave          = convert_air_to_vacuum(wave) # now in vacuum

print "writing result to",outfilename
hdulist.writeto(outfilename,clobber=True)
