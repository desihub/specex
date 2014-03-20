#!/usr/bin/env python

import pyfits,sys,json,pylab,string,numpy,os,scipy,scipy.sparse,scipy.linalg
from scipy.sparse.linalg import spsolve
from math import *

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

#A=scipy.sparse.dia_matrix((nwave,nwave)) # 
A=numpy.matrix(numpy.zeros((nwave,nwave))) # dense because additions of band matrices not implemented
B=numpy.zeros((1,nwave))
for fiber in skyfibers:
    
    R=scipy.sparse.dia_matrix((Rdata[fiber],offsets),(nwave,nwave))
    Ninv=scipy.sparse.dia_matrix((invar[fiber,:],[0]),(nwave,nwave))
    tmp=invar[fiber,:]*spectra[fiber,:]
    tmp2=R.transpose()*Ninv*R
    A+=tmp2.todense()
    B+=R.transpose().dot(tmp)

#convert A to sparse for solving
Ad=d*2
Aoffsets = range(Ad,-Ad-1,-1)
Adata =  numpy.zeros((len(Aoffsets),nwave))  
for i in range(len(Aoffsets)) :
    diagonal=numpy.diag(A,Aoffsets[i])
    off=max(0,Aoffsets[i])
    Adata[i,off:len(diagonal)+off]=diagonal
As=scipy.sparse.dia_matrix( (Adata,Aoffsets), shape=(nwave,nwave))

print "done"

print "solving"
As = As.tocsr()
deconvolvedsky=spsolve(As,B)
print "deconvolvedsky",deconvolvedsky.shape
print "done"

if skyfilename != "" :
    print "writing skymodel to",skyfilename
    R=scipy.sparse.dia_matrix((Rdata[nfibers/2],offsets),(nwave,nwave))
#convolvedsky=R.dot(deconvolvedsky) #error of byte swap I don't understand here 
    sky=numpy.dot(R.todense(),deconvolvedsky) 
    skyinvar=numpy.zeros((1,nwave))
    skyinvar[0,:]=numpy.diag(A).copy()
    pyfits.HDUList([pyfits.PrimaryHDU(sky),pyfits.ImageHDU(skyinvar),pyfits.ImageHDU(wave)]).writeto(skyfilename,clobber=True)
   

print "subtracting sky to all fibers"
for fiber in range(nfibers) :
    R=scipy.sparse.dia_matrix((Rdata[fiber],offsets),(nwave,nwave))
    sky=numpy.dot(R.toarray(),deconvolvedsky) # it is a numpy.matrix that has to be converted to a numpy.array
    spectra[fiber] -= sky
print "done"

print "writing result to",outfilename
hdulist.writeto(outfilename,clobber=True)

sys.exit(0)
