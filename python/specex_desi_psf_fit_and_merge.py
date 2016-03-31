#!/usr/bin/env python

import sys,numpy,pylab
import astropy.io.fits as pyfits
import argparse
import numpy as np
import os.path
import
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--psfboot', type = str, default = None, required=True,
                        help = 'path of DESI PSF boot fits file')
parser.add_argument('--arc', type = str, default = None, required=True,
                        help = 'path of DESI arc image fits file')
args = parser.parse_args()
hdulist=pyfits.open(args.arc)
print hdulist[0].data
cam=hdulist[0].header["CAMERA"].strip()
try :
    expid=hdulist[0].header["EXPID"].strip()
except KeyError :
    expid=0

print expid,cam

params="--gauss_hermite_deg 2 --half_size_x 4 --half_size_y 4 -v"

for bundle in range(20) :
    psfxml="psf-specex-exp%08d-%s-bun%02d.xml"%(expid,cam,bundle)
    psffits="psf-specex-exp%08d-%s-bun%02d.fits"%(expid,cam,bundle)
    spotxml="spots-exp%08d-%s-bun%02d.xml"%(expid,cam,bundle)
    log="specex_desi_psf_fit-exp%08d-%s-bun%02d.log"%(expid,cam,bundle)
    cmd="specex_desi_psf_fit -a %s --first_bundle %d --last_bundle %d --xcoord-file %s --xcoord-hdu 1 --ycoord-file %s --ycoord-hdu 2 --out_xml %s --out_fits %s --out_spots %s %s"%(args.arc,bundle,bundle,args.psfboot,args.psfboot,psfxml,psffits,spotxml,params)
    if os.path.isfile(psfxml) :
        print "skip existing",psfxml
        continue
    print "%s >& %s"%(cmd,log)
    system("%s >& %s"%(cmd,log))



    
