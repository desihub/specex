#!/usr/bin/env python


import pyfits,sys,json,pylab,string,numpy,os,scipy.interpolate
from scipy.stats import norm
from math import *

def rebin(wave,spec,width) :
    n1=wave.shape[0]
    n2=int((wave[-1]-wave[0])/width)
    n2=n1/(n1/n2)
    owave    = wave[0:n1-n1%n2].reshape((n2,n1/n2)).mean(-1)
    ospec    = spec[:,0:n1-n1%n2].reshape((spec.shape[0],n2,n1/n2)).mean(-1)
    return owave,ospec

if len(sys.argv)<3 :
    print sys.argv[0],"flat1.fits flat2.fits"
    sys.exit(12);

hdulist1=pyfits.open(sys.argv[1]);
hdulist2=pyfits.open(sys.argv[2]);

delta=hdulist1[0].data/hdulist2[0].data-1
wave=hdulist1[2].data

pyfits.HDUList([pyfits.PrimaryHDU(delta)]).writeto("delta.fits",clobber=True)

owave,odelta=rebin(wave,delta,500)



save=True
if True :
    fiber=101
    fig0=pylab.figure("flat_fiber%d_vs_wavelength"%(fiber+1))
    pylab.plot(wave,hdulist1[0].data[fiber])
    pylab.plot(wave,hdulist2[0].data[fiber])
    pylab.xlabel('Wavelength ($\AA$)')
    pylab.ylabel('flat-field correction')
    pylab.title(r'fiber %d, exp. 104768 and 104774'%(fiber+1))
    if save :
        fig0.savefig("flat_fiber101_vs_wavelength.pdf")
if True :
    fiber=101
    fig=pylab.figure("toto")
    toto=hdulist1[0].data[fiber]/hdulist2[0].data[fiber]
    var=1/(hdulist1[1].data[fiber]*(hdulist1[0].data[fiber])**2)+1/(hdulist2[1].data[fiber]*(hdulist2[0].data[fiber])**2)
    
    
    n,bins,p= pylab.hist(toto-1,bins=50,range=[-0.03,0.03],histtype='step',normed=1)
    indices=numpy.where(abs(toto-1)<0.015)[0]
    (mu, sigma) = norm.fit(toto[indices]-1)
    print mu,sigma
    gx=numpy.arange(-0.015,0.015,0.0005)
    gy=pylab.normpdf(gx,mu,sigma)
    pylab.plot(gx,gy,'k')
    chi2=numpy.sum((toto[indices]-1-mu)**2/var[indices])
    n=toto.shape[0]
    print "mu=",mu," sigma=",sigma," sigma/sqrt(2)=",sigma/sqrt(2.)
    print "chi2/pdf=",chi2/n
    
if True :
    fig1=pylab.figure("flat_vs_wavelength")
    for fiber in range(delta.shape[0]) :
        #pylab.plot(wave,delta[fiber])
        pylab.plot(owave,odelta[fiber])
        
    pylab.xlabel('Wavelength ($\AA$)')
    pylab.ylabel('$\Delta \log(flat)$')
    if save :
        fig1.savefig("flat_vs_wavelength.pdf")

if True :
    fig2=pylab.figure("flat_vs_fiber")
    pylab.plot(odelta[:,owave.shape[0]/2])
    pylab.xlabel('fiber')
    pylab.ylabel('$\Delta \log(flat)$')
    if save : 
        fig2.savefig("flat_vs_fiber.pdf")
    
if True :
    fig3=pylab.figure("flat_histogram")

    x=odelta.reshape((odelta.shape[0]*odelta.shape[1]))
    #x=odelta[:,owave.shape[0]/2]
    n,bins,p= pylab.hist(x,bins=30,range=[-0.015,0.015],histtype='step',normed=1)

    (mu, sigma) = norm.fit(x[numpy.where(abs(x)<0.008)[0]])
    print mu,sigma
    gx=numpy.arange(-0.015,0.015,0.0005)
    gy=pylab.normpdf(gx,mu,sigma)
    print gx.shape,gy.shape
    pylab.plot(gx,gy,'k')
    pylab.xlabel('$\Delta \log(flat)$')
    pylab.title(r'$\mathrm{Histogram\ of\ }\Delta \log(flat):\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
    if save : 
        fig3.savefig("flat_histogram.pdf")

pylab.show()

