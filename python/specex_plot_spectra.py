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

        self.wave=[]
        if len(hdulist)==2 :
            self.wave =  hdulist[1].data[0]
        elif len(hdulist[2].data.shape)==1 :
            self.wave =  hdulist[2].data
        elif len(hdulist[2].data.shape)==2 and hdulist[2].data.shape[0] == 1   :
            self.wave =  hdulist[2].data[0]
            
        if len(self.wave)==0 :
            print "error couldn't find wave data in file"
            sys.exit(12)

        if len(self.wave) != nwave :
            print "error"
            sys.exit(12)

        self.errors=[]
        ivar=hdulist[1].data
        if ivar.shape==self.spectra.shape :
            self.errors=ivar
            for i in range(self.errors.shape[0]) :
                for j in range(self.errors.shape[1]) :
                    if ivar[i,j]>0 :
                        self.errors[i,j]=1./sqrt(ivar[i,j])
                    else :
                        self.errors[i,j]=0

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
    for fiber in dset.fibers :
        if show_errors and len(dset.errors) :
            pylab.errorbar(dset.wave,dset.spectra[fiber,:],yerr=dset.errors[fiber,:])
        elif show_log :
            pylab.plot(dset.wave,numpy.log(dset.spectra[fiber,:]))
        else :
            pylab.plot(dset.wave,dset.spectra[fiber,:])

pylab.show() # don't need cause "ion"

#import IPython
#IPython.embed()


