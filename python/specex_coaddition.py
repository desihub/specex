#!/usr/bin/env python


import pyfits,sys,json,pylab,string,numpy,os,scipy.interpolate,scipy.linalg
from scipy.sparse.linalg import spsolve
from scipy.stats import norm
from math import *
from specex_cholesky import *
from math import *


def compute_model(flux,invar,calibcoeff,calibderivatives) :
    #print "compute mean model"
    #sys.stdout.flush()   
    
    nexpo=flux.shape[0]
    nwave=flux.shape[1]
        
    a = numpy.zeros(nwave)
    b = numpy.zeros(nwave)    
        
    calibcorr=numpy.ones((nexpo,nwave))
    calt=calibderivatives.transpose()
    for expo in range(nexpo) :
        calibcorr[expo,:] += calt[:].dot(calibcoeff[expo])
        #print calibcorr[expo]
    
    a = numpy.sum(invar*calibcorr**2,axis=0)
    b = numpy.sum(invar*calibcorr*flux,axis=0)
    
    model=b/a
    
    
    #print "done compute mean model"
    #sys.stdout.flush()  
    
    return model

def compute_model_invar(invar,calibderivatives) :
    #print "compute mean model"
    #sys.stdout.flush()   
    
    nexpo=flux.shape[0]
    nwave=flux.shape[1]
        
    a = numpy.zeros(nwave)
    b = numpy.zeros(nwave)    
        
    calibcorr=numpy.ones((nexpo,nwave))
    calt=calibderivatives.transpose()
    for expo in range(nexpo) :
        calibcorr[expo,:] += calt[:].dot(calibcoeff[expo])
        #print calibcorr[expo]
    
    a = numpy.sum(invar*calibcorr**2,axis=0)
    return 1/a
        
    
    #print "done compute mean model"
    #sys.stdout.flush()  
    
    return model

def fit_calibration(wave,flux,invar,model,calibcoeff,calibderivatives,verbose=True) :
    
    #print "model",model.shape
    #print "calibmodel",calibmodel.shape
    #print "calibderivatives",calibderivatives.shape
    
    # flux = model*(1+calibcorrections*calibderivatives)
    # flux = model*calibval
    # flux = calibmodel
    
    nexpo=flux.shape[0]
    nwave=wave.shape[0]
    ncoeff=calibderivatives.shape[0]

    
    calt=calibderivatives.transpose()
    val=numpy.zeros(nwave)
    
    for expo in range(nexpo) :

        
        val[:]=calt[:].dot(calibcoeff[expo])
        calibmodel = model*(1.+val)

        A=numpy.zeros((ncoeff,ncoeff))
        B=numpy.zeros((ncoeff))
        for i in range(ncoeff) :
            ider=calibderivatives[i]
            B[i]=numpy.sum(invar[expo]*(flux[expo]-calibmodel)*ider*model)
            for j in range(i+1) :
                jder=calibderivatives[j]
                A[i,j]=numpy.sum(invar[expo]*ider*jder*model**2)
                A[j,i]=A[i,j]
        
        # new calib corrections :
        calibcoeff[expo] += cholesky_solve(A,B)
    
    
    
    # center this
    for i in range(ncoeff) :
        calibcoeff[:,i] -= numpy.mean(calibcoeff[:,i])
        
    
    # apply to model
    # calt=calibderivatives.transpose()
    # for expo in range(nexpo) :
    #     val=numpy.zeros(nwave)
    #     val[:]=calt[:].dot(calibcoeff[expo])
    
    return calibcoeff
    
    
    
 
def outlier_clipping(wave,flux,invar,model,calibcoeff,calibderivatives,calibmodel,nsig=3.,wave_bin=200.) :
    
    step=(wave[-1]-wave[0])/int((wave[-1]-wave[0])/wave_bin)
    binwave=numpy.arange(wave[0],wave[-1],step)
    
    nexpo=flux.shape[0]
    nwave=wave.shape[0]
    
    # leave a minimum of (nexpo/2+1) in any case
    nexpo_min=int(nexpo/2)+1
    nexpo_per_wave=numpy.sum(invar>0,axis=0) # shape is (nwave)
    
    windices_with_enough_exposures=numpy.where(nexpo_per_wave>nexpo_min)[0] # so we can rm one
    
    #update calib model
    calt=calibderivatives.transpose()
    val=numpy.zeros(nwave)
    for expo in range(nexpo) :
        val[:]=calt[:].dot(calibcoeff[expo])
        calibmodel[expo] = model*(1.+val)

    # compute residuals
    dchi2=invar*(flux-calibmodel)**2 # shape is (nexpo,nwave)
    
    chi2pdf=0
    nout=0
    
    for bw1 in binwave :
        bw2=bw1+step
        
        # compute chi2
        windices=numpy.where((wave>=bw1)&(wave<bw2))[0]
        chi2=numpy.sum(dchi2[:,windices])
        ndata=numpy.sum(invar[:,windices]>0)
        npar=len(windices)
        ndf=ndata-npar
        bin_chi2pdf=chi2/ndf
        

        # remove at max one outlier entry per wave, and leave a minimum of (nexpo/2+1) in any case
        
        # find wave bins where there are outlier and at least nexpomin+1
        windices=numpy.intersect1d(windices,windices_with_enough_exposures)
        
        tmp=numpy.zeros(dchi2.shape[1])
        tmp[windices]=numpy.sum((dchi2[:,windices]>(nsig*bin_chi2pdf)),axis=0)
        
        windices=numpy.intersect1d(windices,numpy.where(tmp>0)[0])
        
        bin_nout=len(windices)
        
        #print bw1,bw2,"chi2pdf=",bin_chi2pdf,ndata,npar,"nout=",bin_nout
        
        for w in windices :
            # find max outlier
            expo=numpy.argmax(dchi2[:,w])
            # set its weight to zero
            invar[expo,w]=0
        
        nout += bin_nout
        chi2pdf += bin_chi2pdf
        
    chi2pdf/=binwave.shape[0] # mean
    
    return chi2pdf,nout
    

if len(sys.argv)<2 :
    print sys.argv[0],"spXvfsc-*.fits"
    sys.exit(12);


print "check and load data"
print "count data"
sys.stdout.flush()
n={'b1':0,'r1':0,'b2':0,'r2':0}
nfibers=0
input_wave={}
input_filename={}
plateid=None
mjd=None

for c in range(2,len(sys.argv)) :
    filename=sys.argv[c]
    print "inspecting",filename
    hdulist=pyfits.open(filename)
    band=hdulist[0].header["CAMERAS"]
    if not n.has_key(band)  :
        print "ERROR unknown band",band,"in",filename
        sys.exit(12)
    n[band] += 1
    

    if plateid == None :
        plateid = hdulist[0].header["PLATEID"]
    else :
        if hdulist[0].header["PLATEID"] != plateid :
            print "ERROR not same plate ",plateid,hdulist[0].header["PLATEID"],"in",filename
            sys.exit(12)
    
    if mjd == None :
        mjd = hdulist[0].header["MJD"]

    if not input_wave.has_key(band) :
        input_wave[band]=hdulist[2].data.copy()
    else :
        # check compatibility
        if numpy.any(input_wave[band] != hdulist[2].data) :
            print "ERROR not wave wavelength grid for band ",band,"in",filename
            sys.exit(12)
    if not input_filename.has_key(band) :
        input_filename[band]=[]
    input_filename[band].append(filename)
    if nfibers==0 :
        nfibers=hdulist[0].data.shape[0]
    else :
        if nfibers != hdulist[0].data.shape[0] :
            print "ERROR not same number of fibers",nfibers,hdulist[0].data.shape[0],"in",filename
            sys.exit(12)
    hdulist.close()

input_flux={}
input_invar={}
model_flux={}
model_invar={}
model_calibcoeff={}
calibrated_model_flux={}
ncoeff=2 # umber of calibration correction coefficients per input spectrum, first order correction
calibration_derivatives={}

bands=input_wave.keys()
for band in bands :
    nexpo=n[band]
    nwave=input_wave[band].shape[0]
    input_flux[band]=numpy.zeros((nexpo,nfibers,nwave))
    input_invar[band]=numpy.zeros((nexpo,nfibers,nwave))
    model_flux[band]=numpy.zeros((nfibers,nwave))
    model_invar[band]=numpy.zeros((nfibers,nwave))
    model_calibcoeff[band]=numpy.zeros((nexpo,nfibers,ncoeff))
    calibrated_model_flux[band]=numpy.zeros((nexpo,nfibers,nwave))
    calibration_derivatives[band]=numpy.zeros((ncoeff,nwave))
    
print "load data"
sys.stdout.flush()

for band in bands :
    n[band]=0

for c in range(2,len(sys.argv)) :
    filename=sys.argv[c]
    print "loading",filename
    hdulist=pyfits.open(filename)
    band=hdulist[0].header["CAMERAS"]
    input_flux[band][n[band]]=hdulist[0].data.copy()
    input_invar[band][n[band]]=hdulist[1].data.copy()
    hdulist.close()
    n[band] += 1
print "done loading data"
sys.stdout.flush()   




print "compute calibration derivatives"
sys.stdout.flush() 
  
for band in bands :
    #nwave=input_wave[band].shape[0]
    mean_wave=(input_wave[band][0]+input_wave[band][-1])/2.
    range_of_wave=(input_wave[band][-1]-input_wave[band][0])/2.
    rwave=(input_wave[band]-mean_wave)/range_of_wave
    der=calibration_derivatives[band]
    for i in range(ncoeff) :
        der[i]=rwave**i
print "done calibration derivatives"
sys.stdout.flush()   

for band in bands :

    #fibers=[178,179]
    fibers=numpy.arange(input_flux[band].shape[1])
    
    for fiber in fibers :

        wave=input_wave[band]
        flux=input_flux[band][:,fiber]
        invar=input_invar[band][:,fiber]
        #model=model_flux[band][fiber]
        calibcoeff=model_calibcoeff[band][:,fiber,:]
        calibmodel=calibrated_model_flux[band][:,fiber,:]
        calibderivatives=calibration_derivatives[band]

        for loop in range(20) :
                        
            # compute model
            model = compute_model(flux=flux,invar=invar,calibcoeff=calibcoeff,calibderivatives=calibderivatives)
            
            # outlier rejection
            chi2pdf,nout=outlier_clipping(wave=wave,flux=flux,invar=invar,model=model,calibcoeff=calibcoeff,calibderivatives=calibderivatives,calibmodel=calibmodel,nsig=4.,wave_bin=500.)
            
            # fit calibration
            calibcoeff=fit_calibration(wave=wave,flux=flux,invar=invar,model=model,calibcoeff=calibcoeff,calibderivatives=calibderivatives,verbose=True)
                        
            line = "#%d %s f=%03d chi2pdf=%4.3f nout=%03d coef="%(loop,band,fiber,chi2pdf,nout)
            for expo in range(nexpo) :
                line += " %f"%calibcoeff[expo,0]
            
            print line
            
            
            if nout==0 :
                break
        

        model_flux[band][fiber]=model
        model_invar[band][fiber]=compute_model_invar(invar=invar,calibderivatives=calibderivatives)
        model_calibcoeff[band][:,fiber,:]=calibcoeff
        
        sys.stdout.flush()
    
        
if False :
    print "write recalibrated data"
    
    for band in bands :

        calibderivatives=calibration_derivatives[band]
        calt=calibderivatives.transpose()

        print calt.shape

        expo=-1
        for ifilename in input_filename[band] :
            expo+=1
            ofilename=string.replace(ifilename,"spXvfsc","spXvfscc")
            print ifilename,"->",ofilename
            sys.stdout.flush()
            
            hdulist=pyfits.open(ifilename);

            print hdulist[0].data.shape

            nfibers=hdulist[0].data.shape[0]
            nwave=hdulist[0].data.shape[1]
            
            for fiber in range(nfibers) :
                calibcoeff=model_calibcoeff[band][expo,fiber,:]
                calibcorr=numpy.ones((nwave))
                calibcorr[:] += calt[:].dot(calibcoeff)
                hdulist[0].data[fiber]/=calibcorr  # flux
                hdulist[1].data[fiber]*=(calibcorr**2)  # invar
                
            hdulist.writeto(ofilename,clobber=True)
            print "wrote",ofilename
            sys.stdout.flush()

            hdulist.close()

if True :
    print "write coadded data per band"
    
    for band in bands :

        print "do the average of the resolution matrices"
        R=None
        for ifilename in input_filename[band] :
            hdulist=pyfits.open(ifilename);
            if R==None :
                R=hdulist[3].data.copy()
            else :
                R+=hdulist[3].data
            hdulist.close()
        R/=len(input_filename[band])
        print "done averaging the resolution matrices"
        
        ofilename = "coadd-%s-%s-%s.fits"%(band,str(plateid),str(mjd))
        
        wave=input_wave[band]
        model=model_flux[band]
        invar=model_invar[band]
        output_hdulist=pyfits.HDUList([pyfits.PrimaryHDU(model),pyfits.ImageHDU(invar,name="IVAR"),pyfits.ImageHDU(wave,name="WAVELENGTH"),pyfits.ImageHDU(R,name="RESOLUTION")])
        ncoef=model_calibcoeff[band].shape[2]
        for i in range(ncoef) :
            output_hdulist.append(pyfits.ImageHDU(model_calibcoeff[band][:,:,i],name="COEF%d"%i))
        
        # add some keys
        output_header=output_hdulist[0].header
        output_header.update("NCOEF",ncoef,"number of recalib. coefficients")
        
        first=True
        index=-1
        for ifilename in input_filename[band] :
            index+=1
            hdulist=pyfits.open(ifilename);
            input_header=hdulist[0].header
            if first :
                for k in input_header.keys() :
                    try :
                        output_header.update(k,input_header[k],"from %s"%ifilename)
                    except :
                        pass
            first=False
            output_header.update("FILE%02d"%index,filename,"used in coadd")
            for k in ["EXPOSURE","MJD","AZ","ALT"] :
                tk=k
                if len(tk)>6 :
                    tk=k[:6]
                nk="%s%02d"%(tk,index)
                output_header.update(nk,input_header[k],"%s key of FILE%02d"%(k,index))
            hdulist.close()

        output_hdulist.writeto(ofilename,clobber=True)
        print "wrote",ofilename
        sys.stdout.flush()
        hdulist.close()

sys.exit(0)
    
