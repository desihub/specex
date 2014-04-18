#!/usr/bin/env python


import pyfits,sys,json,pylab,string,numpy,os,scipy.interpolate,scipy.linalg
from scipy.sparse.linalg import spsolve
from scipy.stats import norm
from math import *
from specex_cholesky import *
from math import *




# global variables
# ----------------------------------------------------
# input data
input_waves=[]
input_fluxes=[]
input_invars=[]

# needed to define reduced wavelength
mean_wave=[]
range_of_wave=[]

# output data
output_wave=None # will be 1D a numpy array
output_flux=None # will be 1D a numpy array
covmat=None      # will be 2D a numpy array
lscale=0         # step in log(wave), NOT log10(wave)
calib_coeffs=None # will be 2D a numpy array

# derivatives
interp_derivatives=[] # derivative wrt to uncalibrated model flux, list per input spec of numpy 2D array
interp_derivatives_indices=[] # first index in output flux of derivative wrt to uncalibrated model flux, list of numpy 1D array
calibration_derivatives=[] # derivative wrt to calibration coefficents

# model in original wave grid
flux_model_input_waves=[]

# parameters 
ncoefperspec  = 2 # umber of calibration correction coefficients per input spectrum, first order correction


# ----------------------------------------------------

def set_output_wave() :
    # define an output wavelength grid, choose log for now
    # from file spec-3647-55181-0180.fits , same as BOSS pipeline
    log10bin=9.9897385e-05
    
    # check this using bin size of spectra
    nspec = len(input_waves)

    for s in range(nspec) :
        slog10bin=log10(input_waves[s][1]/input_waves[s][0])
        if slog10bin>log10bin :
            print "need to increase bin size to %f instead of %f"%(slog10bin,log10bin)
            log10bin=slog10bin

    scale=10**log10bin
    global lscale
    lscale=log(scale)
    nbins=int((log10(iwave[-1])-log10(iwave[0]))/log10bin)+1
    print "nbins=",nbins
    global output_wave
    output_wave=numpy.zeros((nbins))
    output_wave[0]=iwave[0]
    for i in range(1,nbins) :
        output_wave[i]=output_wave[i-1]*scale
        
    


def compute_derivatives() :
    # precompute non variable derivatives
    nspec = len(input_waves)

    for s in range(nspec) :

        # wave
        iwave=input_waves[s]
        nw=iwave.shape[0]

        # reduced wavelength
        mean_wave.append((iwave[0]+iwave[-1])/2.)
        range_of_wave.append(iwave[-1]-iwave[0])
        
        rwave=2.*(iwave-mean_wave[s])/range_of_wave[s] # range is now -1,1
        
        # monomials of a polynomial of reduced wavelength
        spec_calibration_derivatives=numpy.zeros((iwave.shape[0],ncoefperspec))
        for i in range(nw) :
            spec_calibration_derivatives[i]=rwave[i]**numpy.arange(ncoefperspec) 
        calibration_derivatives.append(spec_calibration_derivatives)

        # derivatives of model flux at this wave wrt to model fluxes on nodes
        spec_interp_derivatives_indices = numpy.zeros((iwave.shape[0])).astype(int) # only first index is stored
        spec_interp_derivatives         = numpy.zeros((iwave.shape[0],2)) # two derivatives per flux point

        for i in range(nw) :
            oi=int(log(iwave[i]/output_wave[0])/lscale)
            spec_interp_derivatives_indices[i]=oi
            if oi>=output_wave.shape[0]-1 :
                continue # will fix this later
            # derivatives of interpolation
            spec_interp_derivatives[i,0]=(output_wave[oi+1]-iwave[i])/(output_wave[oi+1]-output_wave[oi])
            spec_interp_derivatives[i,1]=(iwave[i]-output_wave[oi])/(output_wave[oi+1]-output_wave[oi])
        interp_derivatives_indices.append(spec_interp_derivatives_indices)
        interp_derivatives.append(spec_interp_derivatives)

    global calib_coeffs
    calib_coeffs = numpy.zeros((nspec,ncoefperspec))
    
    # ----------------------------------------------------------------------     

def fit_model(verbose=True,covar=False) :
    if verbose :
        print "filling A B ... "
    npar=output_wave.shape[0]
    A=numpy.zeros((npar,npar))
    B=numpy.zeros((npar))

    nspec=len(input_waves)

    for s in range(nspec) :
        nw=input_waves[s].shape[0]
        for i in range(nw) :
            oi=interp_derivatives_indices[s][i]
            if oi>=output_wave.shape[0]-1 :
                continue # will fix this later
            
            cder=calibration_derivatives[s][i]
            calibcorrection=(1+cder.dot(calib_coeffs[s]))
            der=calibcorrection*interp_derivatives[s][i]
            
            #print "der.shape",der.shape
            B[oi:oi+2] += (input_invars[s][i]*input_fluxes[s][i])*der[:]
            A[oi:oi+2,oi:oi+2] += input_invars[s][i]*numpy.outer(der,der)
    
    global output_flux
    global covmat

    if not covar : 
        if verbose :
            print "conversion to band matrix ..."
        d=3
        offsets = range(d,-d-1,-1)
        data =  numpy.zeros((len(offsets),npar))  
        for i in range(len(offsets)) :
            diagonal=numpy.diag(A,offsets[i])
            off=max(0,offsets[i])
            data[i,off:len(diagonal)+off]=diagonal
        As=scipy.sparse.dia_matrix((data,offsets),shape=(npar,npar))
        As = As.tocsr()
        if verbose :
            print "solve sparse ..."     
        output_flux=scipy.sparse.linalg.spsolve(As,B)
    else :
        if verbose :
            print "solve and invert with cholesky ..."  
        output_flux,covmat=cholesky_solve_and_invert(A,B)
    
    if verbose :
        print "done"

def fit_calibration(verbose=True) :
        
    nspec=len(input_waves)
    
    
    

    for s in range(nspec) :
        nw=input_waves[s].shape[0]
        npar=ncoefperspec
        A=numpy.zeros((npar,npar))
        B=numpy.zeros((npar))
        
        for i in range(nw) :
            oi=interp_derivatives_indices[s][i]
            if oi>=output_wave.shape[0]-1 :
                continue # will fix this later

            der=interp_derivatives[s][i]
            model=der.dot(output_flux[oi:oi+2])    # uncalibrated model value    
            #print "model",model
            cder=calibration_derivatives[s][i]
            #print "cder",cder
            cmodel=model*(1+cder.dot(calib_coeffs[s])) # calibrated model value
            #print "cmodel",cmodel
            res=input_fluxes[s][i]-cmodel
            #print "cres",res
            cdermodel=model*cder
            #print "cdermodel",cdermodel
            
            B += input_invars[s][i]*res*cdermodel
            A += input_invars[s][i]*numpy.outer(cdermodel,cdermodel)
            
            
        # solving
        calcoeff=cholesky_solve(A,B)
        calib_coeffs[s]+=calcoeff
        if verbose :
            print "spec",s,"cal. coeffs",calib_coeffs[s]

def update_model_per_spec() :
    
    global flux_model_input_waves
    flux_model_input_waves=[]
    
    nspec=len(input_waves)
    for s in range(nspec) :
        nw=input_waves[s].shape[0]
        spec_flux_model_input_waves=numpy.zeros((nw))
        for i in range(nw) :
            oi=interp_derivatives_indices[s][i]
            if oi>=output_wave.shape[0]-1 :
                continue # will fix this later
            der=interp_derivatives[s][i]
            model=der.dot(output_flux[oi:oi+2])    # uncalibrated model value    
            cder=calibration_derivatives[s][i]
            spec_flux_model_input_waves[i]=model*(1+cder.dot(calib_coeffs[s])) # calibrated model value

        flux_model_input_waves.append(spec_flux_model_input_waves)
    

# returns mean chi2/ndf,nout,ndata
def outlier_clipping(nsig=3.,wave_bin=200.) :
    
    step=(output_wave[-1]-output_wave[0])/int((output_wave[-1]-output_wave[0])/wave_bin)
    binwave=numpy.arange(output_wave[0],output_wave[-1],step)
    

    ndata_tot=0
    nout_tot=0
    mean_chi2pdf=0

    for bw1 in binwave :
        bw2=bw1+step
        
        # compute chi2
        chi2=0
        ndata=0
        for s in range(nspec) :
            indices=numpy.where((input_waves[s]>=bw1)&(input_waves[s]<bw2)&(input_invars[s]>0))[0]
            n=len(indices)
            if n>0:
                chi2+=numpy.sum(input_invars[s][indices]*(input_fluxes[s][indices]-flux_model_input_waves[s][indices])**2)
                ndata+=n
        npar=len(numpy.where((output_wave>=bw1)&(output_wave<bw2))[0])
        ndf=ndata-npar
        chi2pdf=chi2/ndf
        #print bw1,bw2,chi2pdf
        
        # discard nsig outliers
        nout=0
        for s in range(nspec) :
            indices=numpy.where((input_waves[s]>=bw1)&(input_waves[s]<bw2)&(input_invars[s]>0)&(input_invars[s]*(input_fluxes[s]-flux_model_input_waves[s])**2>nsig*chi2pdf))[0]
            n=len(indices)
            #print "wavebin",bw1,bw2,"spec=",s,"nout=",n,"ntot=",w.shape[0]
            if n==0 :
                continue
            input_invars[s][indices]=0 # set to zero weight of input var
            nout+=n
        #print  "wavebin",bw1,bw2,"nout=",nout,"ndata=",ndata
        nout_tot+=nout
        ndata_tot+=ndata
        mean_chi2pdf+=chi2pdf

    mean_chi2pdf/=binwave.shape[0]

    return mean_chi2pdf,nout_tot,ndata_tot


if len(sys.argv)<3 :
    print sys.argv[0],"output.fits spXvfsc-*.fits"
    sys.exit(12);


output_filename=sys.argv[1]



fibers=[178,179]
output_waves=[]
output_fluxes=[]
output_invars=[]


for fiber in fibers :
    print "Adding spectra of fiber",fiber
    
    input_waves=[]
    input_fluxes=[]
    input_invars=[]
    mean_wave=[]
    range_of_wave=[]
    interp_derivatives=[]
    interp_derivatives_indices=[]
    calibration_derivatives=[]
    flux_model_input_waves=[]
    output_wave=None
    output_flux=None
    covmat=None  
    calib_coeffs=None
    
    for c in range(2,len(sys.argv)) :
        filename=sys.argv[c]
        hdulist=pyfits.open(filename);
        band=hdulist[0].header["CAMERAS"][0]
        print band,filename

        input_waves.append(hdulist[2].data)
        input_fluxes.append(hdulist[0].data[fiber,:])
        input_invars.append(hdulist[1].data[fiber,:])
        hdulist.close()
    
    nspec=len(input_waves)
    # first create a merged list of wave
    iwave=numpy.zeros((0))
    for s in range(nspec) :
        iwave=numpy.union1d(iwave,input_waves[s])
    print iwave.shape,iwave[0],iwave[-1]

    
    
    set_output_wave()    
    compute_derivatives()
    
    # now we are going to fit the data assuming the model is a piece-wise spectral density
    # this computes output_flux
    for loop in range(10) :
        fit_model(verbose=False)
        update_model_per_spec()
        fit_calibration(verbose=False)
        chi2pdf,nout,ndata=outlier_clipping(nsig=5.,wave_bin=200.)
        print "loop %d chi2pdf=%f nout=%d/%d"%(loop,chi2pdf,nout,ndata) 
        if nout==0 :
            break
    
    # compute_covariance
    fit_model(verbose=True,covar=True)
    
    output_invar=numpy.zeros((output_wave.shape[0]))
    if not covmat==None :
        for i in range(output_wave.shape[0]) :
            output_invar[i]=1/covmat[i,i]
        
    output_waves.append(output_wave)
    output_fluxes.append(output_flux)
    output_invars.append(output_invar)
    


# saving results
nfibers=len(fibers)
# we cannot keep a single wavelenght, but we need to keep the same wavelength range
nwave=0
for fiber in range(nfibers) :
    nwave=max(nwave,output_waves[fiber].shape[0])

wave=numpy.zeros((nfibers,nwave))
flux=numpy.zeros((nfibers,nwave))
invar=numpy.zeros((nfibers,nwave))

for fiber in range(nfibers) :
    nw=output_waves[fiber].shape[0]
    wave[fiber,:nw]=output_waves[fiber]
    flux[fiber,:nw]=output_fluxes[fiber]
    invar[fiber,:nw]=output_invars[fiber]

first_filename=sys.argv[2]
print "copying header from",first_filename
first_hdulist=pyfits.open(first_filename);
input_header=first_hdulist[0].header

output_hdulist=pyfits.HDUList([pyfits.PrimaryHDU(flux),pyfits.ImageHDU(invar,name="IVAR"),pyfits.ImageHDU(wave,name="WAVELENGTH")])


print input_header
print input_header["TELESCOP"]
output_header=output_hdulist[0].header
print output_header

# copy header of first hdulist
for k in input_header.keys() :
   output_header.update(k,input_header[k],"from %s"%first_filename)
   

keys=["CAMERAS","EXPOSURE","MJD","AZ","ALT"]

for c in range(2,len(sys.argv)) :
    index=c-2
    filename=sys.argv[c]
    hdulist=pyfits.open(filename)
    header=hdulist[0].header
    output_header.update("FILE%d"%index,filename,"used in coadd")
    
    for k in  keys :
        tk=k
        if len(tk)>7 :
            tk=k[:7]
        nk="%s%d"%(tk,index)
        output_header.update(nk,header[k],"%s key of FILE%d"%(k,index))

    hdulist.close()   
    


output_hdulist.writeto(output_filename,clobber=True)

print "wrote",output_filename

sys.exit(0)
    
