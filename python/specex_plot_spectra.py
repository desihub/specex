#!/usr/bin/env python

import pyfits,sys,json,pylab,string,numpy
from math import *

#pylab.ion() #??

class dataset :
    
    def __init__(self,filename,fibers=None) :
        
        self.filename = filename
        self.fibers = fibers
        
        hdulist=pyfits.open(self.filename)
        hdulist.info()

        self.spectra=hdulist[0].data
        (nspec,nwave) = self.spectra.shape

        if len(fibers)==0 :
            fibers=range(nspec)

        self.fibers = fibers

        try :
            self.wave=hdulist["WAVELENGTH"].data
        except :
            for hdu in range(1,len(hdulist)) :
                data=hdulist[hdu].data
                if len(data.shape) == 1 :
                    self.wave=hdulist[hdu].data
                    break
        
        
        ivar=None
        try :
            ivar=hdulist["IVAR"].data
        except :
            for hdu in range(1,len(hdulist)) :
                data=hdulist[hdu].data
                if len(data.shape) == 2 :
                    ivar=hdulist[hdu].data
                    break
        
        self.errors=[]
        if ivar!=None and ivar.shape==self.spectra.shape :
            self.errors=ivar
            for i in range(self.errors.shape[0]) :
                for j in range(self.errors.shape[1]) :
                    if ivar[i,j]>0 :
                        self.errors[i,j]=1./sqrt(ivar[i,j])
                    else :
                        self.errors[i,j]=0
                    
    def get_wave(self,fiber) :
        if len(self.wave.shape)==1 :
            return self.wave
        else :
            return self.wave[fiber]

def usage() :
    print sys.argv[0],"spec.fits fiber1 fiber2 fiber3:fiber4 spec2.fits fiber5 fiber6 ... (-e)"
    sys.exit(12);

if len(sys.argv)<2 :
    usage()

specfilenames=[]
datasets=[]

filename=""
fibers=[]

show_errors=False
show_log=False

for arg in sys.argv[1:] :
    
    if arg[0]=='-' :
        if arg=="-e" :
            show_errors=True
            continue
        elif arg=="-h" :
            usage()
        else :
            print "unknown option:",arg
            usage()

    isfiber=False
    try :
        i=string.atoi(arg[0])
        isfiber=True
    except :
        pass
    
    if not isfiber :
        
        specfilenames.append(arg)
        
    else :
        
        vals=string.split(arg,":")
        if len(vals)==1 :
            fibers.append(string.atoi(arg))
        elif len(vals)==2 :
            fiber1=string.atoi(vals[0])
            fiber2=string.atoi(vals[1])
            for fiber in range(fiber1,fiber2) :
                fibers.append(fiber)
    

for filename in specfilenames :
    toto=dataset(filename,fibers)
    datasets.append(toto)


for dset in datasets :
    print dset.filename,dset.fibers

    if dset.filename==datasets[-1].filename :
        color="r"
    else :
        color="Gray"
        
    for fiber in dset.fibers :
        if show_errors and len(dset.errors) :
            pylab.errorbar(dset.get_wave(fiber),dset.spectra[fiber,:],yerr=dset.errors[fiber,:],color=color)
        elif show_log :
            pylab.plot(dset.get_wave(fiber),numpy.log(dset.spectra[fiber,:]),color=color)
        else :
            pylab.plot(dset.get_wave(fiber),dset.spectra[fiber,:],color=color)

pylab.show() # don't need cause "ion"

#import IPython
#IPython.embed()


