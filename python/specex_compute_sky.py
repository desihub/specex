#!/usr/bin/env python

import pyfits,sys,json,pylab,string,numpy,os
from math import *


if len(sys.argv)<3 :
    print sys.argv[0],"inspec.fits plPlugMapM.par outspec.fits"
    sys.exit(12);

infilename=sys.argv[1]
plplgmap=sys.argv[2]
outfilename=sys.argv[3]

# get spectrograph id, hardcoded for now, will be read in fits
specid=1

# find sky fibers 
skyfibers=[]
file=open(plplgmap)
for line in file.readlines() :
    if line.find("PLUGMAPOBJ") != 0 : 
        continue
    
    
    vals=string.split(line," ")
    holetype=vals[8]
    if holetype != "OBJECT" :
        continue
    objType=vals[21]
    if objType != "SKY" :
        continue
    spectrographId=string.atoi(vals[24])
    if spectrographId != specid :
        continue
    fiberId=string.atoi(vals[25])
    print line
    print objType,spectrographId,fiberId
    myfiberid=fiberId-1
    if specid==2 :
        myfiberid-=500
    skyfibers.append(myfiberid)
file.close()

print "skyfibers (now starting at 0)=",skyfibers


hdulist=pyfits.open(infilename)
spectra=hdulist[0].data
invar=hdulist[1].data
wave=hdulist[2].data

skyspectra=spectra[skyfibers,:]
skyinvar=invar[skyfibers,:]

if True :
    #for i in skyfibers :
    #    pylab.plot(wave,spectra[i,:])
    #pylab.figure()
    for i in range(skyspectra.shape[0]) :
        pylab.plot(wave,skyspectra[i,:])
    pylab.show()

