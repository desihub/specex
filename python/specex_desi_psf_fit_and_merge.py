#!/usr/bin/env python

import sys,numpy,pylab
import astropy.io.fits as pyfits
import argparse
import numpy as np
import os,os.path
import subprocess

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--psfboot', type = str, default = None, required=True,
                        help = 'path of DESI PSF boot fits file')
parser.add_argument('--arc', type = str, default = None, required=True,
                        help = 'path of DESI arc image fits file')
parser.add_argument('--nthreads', type = int, default = 8, required=False,
                        help = 'number of threads per process')
parser.add_argument('--nproc', type = int, default = 4, required=False,
                        help = 'number of processes')
parser.add_argument('--out_xml', type = str, default = None, required=True,
                        help = 'output PSF file in xml')
parser.add_argument('--out_fits', type = str, default = None, required=True,
                        help = 'output PSF file in fits')

args = parser.parse_args()
hdulist=pyfits.open(args.arc)
cam=hdulist[0].header["CAMERA"].strip()
try :
    expid=hdulist[0].header["EXPID"]
except KeyError :
    expid=0

print "EXPID=%d CAM=%s"%(expid,cam)
print "NPROC=%d"%args.nproc
print "NTHREADS=%d"%args.nthreads



params="--gauss_hermite_deg 2 --half_size_x 4 --half_size_y 4 -v"

pids=[]
for bundle in range(20) :
    psfxml="psf-specex-exp%08d-%s-bun%02d.xml"%(expid,cam,bundle)
    psffits="psf-specex-exp%08d-%s-bun%02d.fits"%(expid,cam,bundle)
    spotxml="spots-exp%08d-%s-bun%02d.xml"%(expid,cam,bundle)
    log="specex_desi_psf_fit-exp%08d-%s-bun%02d.log"%(expid,cam,bundle)
    cmd="export OMP_NUM_THREADS=%d ; specex_desi_psf_fit -a %s --first_bundle %d --last_bundle %d --xcoord-file %s --xcoord-hdu 1 --ycoord-file %s --ycoord-hdu 2 --out_xml %s --out_fits %s --out_spots %s %s"%(args.nthreads,args.arc,bundle,bundle,args.psfboot,args.psfboot,psfxml,psffits,spotxml,params)
    if os.path.isfile(psfxml) :
        print "skip existing",psfxml
        continue
    print "%s >& %s"%(cmd,log)
    p = subprocess.Popen("%s >& %s"%(cmd,log), stdin=None, stdout=None, shell=True)
    pids.append(p)
    print "nproc=%d"%len(pids)
    if len(pids) >= args.nproc :
        print "Now waiting"
         
        for i in xrange(len(pids)) :
            pids[i].wait()
        pids=[]

# wait for the last ones to end
for i in xrange(len(pids)) :
    pids[i].wait()

print "merging results in %s and %s"%(args.out_xml,args.out_fits)
cmd="specex_merge_psf --out-fits %s --out-xml %s"%(args.out_fits,args.out_xml)
for bundle in range(20) :
    bundlefile="psf-specex-exp%08d-%s-bun%02d.xml"%(expid,cam,bundle)
    if os.path.isfile(bundlefile) :
        cmd += " %s"%bundlefile
    else :
        print "WARNING : missing '%s'"%bundlefile

print cmd
os.system(cmd)
print "done"
        


    
