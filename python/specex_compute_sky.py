#!/usr/bin/env python

import pyfits,sys,json,pylab,string,numpy,os,scipy,scipy.sparse,scipy.linalg
from scipy.sparse.linalg import spsolve
from math import *

from specex_cholesky import *


if len(sys.argv)<3 :
    print sys.argv[0],"inspec.fits plPlugMapM.par outspec.fits (sky.fit)"
    sys.exit(12);

infilename=sys.argv[1]
plplgmap=sys.argv[2]
outfilename=sys.argv[3]
skyfilename=""
if(len(sys.argv)>3) :
    skyfilename=sys.argv[4]

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
    #print line
    #print objType,spectrographId,fiberId
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
Rdata=hdulist[3].data

skyspectra=spectra[skyfibers,:]
skyinvar=invar[skyfibers,:]

nskyfibers=len(skyfibers)
nfibers=Rdata.shape[0]
d=Rdata.shape[1]/2
nwave=Rdata.shape[2]
offsets = range(d,-d-1,-1)

print "solving for the mean deconvolved sky"
print "filling A and B"

A=numpy.matrix(numpy.zeros((nwave,nwave))) # dense because additions of band matrices not implemented
B=numpy.zeros((1,nwave))
for fiber in skyfibers:
    
    R=scipy.sparse.dia_matrix((Rdata[fiber],offsets),(nwave,nwave))
    Ninv=scipy.sparse.dia_matrix((invar[fiber,:],[0]),(nwave,nwave))
    tmp=invar[fiber,:]*spectra[fiber,:]
    tmp2=R.transpose()*Ninv*R
    A+=tmp2.todense()
    B+=R.transpose().dot(tmp)

print "done"

print "solving"
deconvolvedsky,dskycovmat=cholesky_solve_and_invert(A,B[0])
print "done"

# compute only once the sky variance because expensive and in any case approximate because we only keep the diagonal
# most conservative is to evaluate it at the highest resolution (most variance)
# also the *mean* sky statistical uncertainty is negligible wrt to the Poisson noise of the subtracted sky of each
# fiber that is already included in the invar of each spectrum
# last point, the sky error is certainly dominated by sky variations in field of view that we have neglected 
R=scipy.sparse.dia_matrix((Rdata[nfibers/2],offsets),(nwave,nwave))
Rt=R.transpose()
sky=numpy.dot(R.toarray(),deconvolvedsky)
print "computing covmat" 
skycovmat=Rt.dot(Rt.dot(dskycovmat).transpose())
skyvar=numpy.diag(skycovmat)
print "done" 

if skyfilename != "" :
    print "writing skymodel to",skyfilename  
    skyinvar=1/numpy.diag(skycovmat)
    sky_array=numpy.zeros((1,sky.shape[0]))
    sky_array[0]=sky
    skyinvar_array=numpy.zeros((1,skyinvar.shape[0]))
    skyinvar_array[0]=skyinvar    
    pyfits.HDUList([pyfits.PrimaryHDU(sky_array),pyfits.ImageHDU(skyinvar_array),pyfits.ImageHDU(wave)]).writeto(skyfilename,clobber=True)
    #pyfits.HDUList([pyfits.PrimaryHDU(skycovmat)]).writeto("skycovmat.fits",clobber=True)

print "subtracting sky to all fibers"
for fiber in range(nfibers) :
    R=scipy.sparse.dia_matrix((Rdata[fiber],offsets),(nwave,nwave))
    Rt=R.transpose()
    sky=numpy.dot(R.toarray(),deconvolvedsky) # it is a numpy.matrix that has to be converted to a numpy.array
    
    spectra[fiber] -= sky
    invar[fiber] = 1/( 1/invar[fiber] + skyvar )

print "done"

print "writing result to",outfilename
hdulist.writeto(outfilename,clobber=True)

sys.exit(0)
