#!/usr/bin/python




import pyfits,sys,json,pylab,string
from math import *

#pylab.ion() #??

class dataset :
    
    def __init__(self,filename,fibers) :
        
        self.filename = filename
        self.fibers = fibers
        
        hdulist=pyfits.open(self.filename)
        hdulist.info()

        self.spectra=hdulist[0].data
        (nspec,nwave) = self.spectra.shape

        self.wave=[]
        if hdulist[1].data.shape[0] == 1 :
            self.wave =  hdulist[1].data[0]
        elif len(hdulist[2].data.shape)==1 :
            self.wave =  hdulist[2].data

        if len(self.wave)==0 :
            print "error"
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
                    self.errors[i,j]=1./sqrt(ivar[i,j])


if len(sys.argv)<3 :
    print sys.argv[0],"spec.fits fiber1 fiber2 fiber3:fiber4 spec2.fits fiber5 fiber6 ..."
    sys.exit(12);


datasets=[]

filename=""
fibers=[]

for arg in sys.argv[1:] :
    
    isfiber=False
    try :
        i=string.atoi(arg[0])
        isfiber=True
    except :
        pass
    
    
    

    if isfiber and filename=="" :
        print "arg error" 
        sys.exit(12)
    
    if not isfiber :
        
        if filename != "" :
            toto=dataset(filename,fibers)
            datasets.append(toto)
            fibers=[]

        filename = arg

    else :
        
        vals=string.split(arg,":")
        if len(vals)==1 :
            fibers.append(string.atoi(arg))
        elif len(vals)==2 :
            fiber1=string.atoi(vals[0])
            fiber2=string.atoi(vals[1])
            for fiber in range(fiber1,fiber2) :
                fibers.append(fiber)
        

toto=dataset(filename,fibers)
datasets.append(toto)


for dset in datasets :
    print dset.filename,dset.fibers



for dset in datasets :
    for fiber in dset.fibers :
        if len(dset.errors) :
            pylab.errorbar(dset.wave,dset.spectra[fiber,:],yerr=dset.errors[fiber,:])
        else :
            pylab.plot(dset.wave,dset.spectra[fiber,:])

pylab.show() # don't need cause "ion"

import IPython
#IPython.embed()


