#!/usr/bin/env python


import pyfits,sys,json,pylab,string,numpy,os,scipy.interpolate,scipy.linalg
from scipy.sparse.linalg import spsolve
from scipy.stats import norm
from math import *
from specex_cholesky import *
from math import *

def compute_wave_split(input_wave,input_invar) :
    # compute a common wavelength grid
    # switch from b to r resolution at the point where inverse variance are roughly equal
    # because we want to minimize resampling at this stage

    mean_invar_b=0.5*(numpy.mean(input_invar["b1"],axis=0)+numpy.mean(input_invar["b2"],axis=0))
    mean_invar_r=0.5*(numpy.mean(input_invar["r1"],axis=0)+numpy.mean(input_invar["r2"],axis=0))

    wave_b=input_wave["b1"]
    wave_r=input_wave["r1"]

    # interpolate mean_invar_b on r wave
    intersect_indices_r=numpy.where(wave_r<wave_b[-1])[0]
    intersect_indices_b=numpy.where(wave_b>wave_r[0])[0]
    
    # fit both with smooth func to avoid problems due to very low invar at some wavelength
    wmin=wave_r[intersect_indices_r][0]
    wmax=wave_r[intersect_indices_r][-1]
    wstep=50
    knots=wmin+wstep/2+wstep*numpy.arange(int((wmax-wmin)/wstep))

    bfunc=scipy.interpolate.splrep(wave_b[intersect_indices_b],mean_invar_b[intersect_indices_b],task=-1,t=knots)
    rfunc=scipy.interpolate.splrep(wave_r[intersect_indices_r],mean_invar_r[intersect_indices_r],task=-1,t=knots)

    wstep=0.5
    knots=wmin+wstep/2+wstep*numpy.arange(int((wmax-wmin)/wstep))
    mean_invar_b_smooth=scipy.interpolate.splev(knots,bfunc,der=0)
    mean_invar_r_smooth=scipy.interpolate.splev(knots,rfunc,der=0)
    wavesplit=knots[numpy.where(mean_invar_r_smooth>mean_invar_b_smooth)[0][0]]
    return wavesplit

def compute_interpolation_derivatives(wave_data,wave_model) :
    n=wave_data.shape[0]
    
    der=numpy.zeros((n,2))
    indices=numpy.zeros((n)).astype(int)
    
    for i in range(1,n-1) :
        i2=numpy.where(wave_model>wave_data[i])[0][0]
        i1=i2-1
        indices[i]=i1
        stepinv=1/(wave_model[i2]-wave_model[i1])
        der[i,0]=stepinv*(wave_model[i2]-wave_data[i])
        der[i,1]=stepinv*(wave_data[i]-wave_model[i1])
    der[0,0]=1
    der[n-1,0]=1
    indices[0]=0
    indices[n-1]=wave_model.shape[0]-2
    
    #print indices

    return der,indices



if len(sys.argv)<5 :
    print sys.argv[0],"coadd-b1-xxxx-xxxxx.fits coadd-b2-xxxx-xxxxx.fits coadd-r1-xxxx-xxxxx.fits coadd-r2-xxxx-xxxxx.fits"
    sys.exit(12);


print "check and load data"
print "count data"
sys.stdout.flush()
n={'b1':0,'r1':0,'b2':0,'r2':0}
nfibers=0
input_wave={}
input_flux={}
input_invar={}
input_filename={}
plateid=None
mjd=None

for c in range(1,len(sys.argv)) :
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
    else :
        if hdulist[0].header["MJD"] != mjd :
            print "ERROR not same mjd ",plateid,hdulist[0].header["MJD"],"in",filename
            sys.exit(12)
    
    if input_wave.has_key(band) :
        print "ERROR, expect only one input file per band"
        sys.exit(12)

    input_wave[band]=hdulist["WAVELENGTH"].data.copy()
    input_flux[band]=hdulist[0].data.copy()
    input_invar[band]=hdulist["IVAR"].data.copy()
    input_filename[band]=filename
    hdulist.close()

# check bands
for band in ["b1","b2","r1","r2"] :
    if not input_wave.has_key(band) :
        print "ERROR missing band",band
        sys.exit(12)

# assume wavelength grid of both spectrographs is the same, we check this here
if numpy.any(input_wave["b1"] != input_wave["b2"]) :
    print "ERROR not same wavelength grid for b1 and b2, this is an assumption of this code"
    sys.exit(12)
if numpy.any(input_wave["r1"] != input_wave["r2"]) :
    print "ERROR not same wavelength grid for r1 and r2, this is an assumption of this code"
    sys.exit(12)


# compute a common wavelength grid
# switch from b to r resolution at the point where inverse variance are roughly equal
# because we want to minimize resampling at this stage
wavesplit = compute_wave_split(input_wave,input_invar)
wave_b = input_wave["b1"]
wave_r = input_wave["r1"]

b_operlap_indices=numpy.where(wave_b>wave_r[0])[0]
r_operlap_indices=numpy.where(wave_r<wave_b[-1])[0]
#owave=wave[g_operlap_indices]
owave_b=wave_b[b_operlap_indices]
owave_r=wave_r[r_operlap_indices]
nwb=b_operlap_indices.shape[0]
nwr=r_operlap_indices.shape[0]
nfibers=input_flux["b1"].shape[0]


print "compute derivatives"

# calib_b(lambda)=1+coef_b*((lambda-lambda_min)/(lambda_max-lambda_min))**2 for b
# calib_r(lambda)=1+coef_r*((lambda-lambda_max)/(lambda_min-lambda_max))**2 for r

wmin=owave_b[0]
wmax=owave_r[-1]
wstep=10
knots=wmin+wstep/2+wstep*numpy.arange(int((wmax-wmin)/wstep))
n=knots.shape[0]

# in the r-band wave grid for both
calder_b=((knots-wmin)/(wmax-wmin))**2
calder_r=((wmax-knots)/(wmax-wmin))**2

# in original wave grid
calder_b_original=numpy.zeros(wave_b.shape)
calder_b_original[b_operlap_indices]=((wave_b[b_operlap_indices]-wmin)/(wmax-wmin))**2
cal_b_original=numpy.zeros(wave_b.shape)

calder_r_original=numpy.zeros(wave_r.shape)
calder_r_original[r_operlap_indices]=((wave_r[r_operlap_indices]-wmin)/(wmax-wmin))**2
cal_r_original=numpy.zeros(wave_r.shape)


print calder_b
print calder_r

split_index=numpy.where(knots>wavesplit)[0][0]
calder_b_split=calder_b[split_index]
calder_r_split=calder_r[split_index]
print "calib. derivatives at split index=",calder_b_split,calder_r_split

# cross-calibrate at the overlap region
# we just want to calculate the calibration coeff. here so we brutally interpolate data on the same grid, that of r
# which has lower resolution

coef={}
for band in ["b1","b2","r1","r2"] :
    coef[band]=numpy.zeros((nfibers))


for spectro in range(1,3) :
    for fiber in range(nfibers) :
        
        flux_b_original=input_flux["b%d"%spectro][fiber,b_operlap_indices]
        invar_b_original=input_invar["b%d"%spectro][fiber,b_operlap_indices]

        if numpy.sum(invar_b_original)==0 :
            print "ignore NULL fiber",fiber
            input_flux["b%d"%spectro][fiber]=0
            input_invar["b%d"%spectro][fiber]=0
            input_flux["r%d"%spectro][fiber]=0
            input_invar["r%d"%spectro][fiber]=0
            continue
        
        tck=scipy.interpolate.splrep(owave_b,flux_b_original,w=invar_b_original,task=-1,t=knots,k=1)
        flux_b=scipy.interpolate.splev(knots,tck,der=0)
        tck=scipy.interpolate.splrep(owave_b,invar_b_original,w=invar_b_original,task=-1,t=knots,k=1)
        invar_b=(nwb/n)*scipy.interpolate.splev(knots,tck,der=0)
        
        flux_r_original=input_flux["r%d"%spectro][fiber,r_operlap_indices]
        invar_r_original=input_invar["r%d"%spectro][fiber,r_operlap_indices]

        if numpy.sum(invar_r_original)==0 :
            print "ignore NULL fiber",fiber
            input_flux["b%d"%spectro][fiber]=0
            input_invar["b%d"%spectro][fiber]=0
            input_flux["r%d"%spectro][fiber]=0
            input_invar["r%d"%spectro][fiber]=0
            continue

        tck=scipy.interpolate.splrep(owave_r,flux_r_original,w=invar_r_original,task=-1,t=knots,k=1)
        flux_r=scipy.interpolate.splev(knots,tck,der=0)
        tck=scipy.interpolate.splrep(owave_r,invar_r_original,w=invar_r_original,task=-1,t=knots,k=1)
        invar_r=(nwr/n)*scipy.interpolate.splev(knots,tck,der=0)
        
        
        if False :
            pylab.errorbar(knots,flux_b,1/numpy.sqrt(invar_b))
            pylab.errorbar(knots,flux_r,1/numpy.sqrt(invar_r))
            #pylab.errorbar(owave_b,flux_b_original,1/numpy.sqrt(invar_b_original))
            #pylab.errorbar(owave_r,flux_r_original,1/numpy.sqrt(invar_r_original))
            pylab.show()
            sys.exit(12)
            
        coef_b=0.
        coef_r=0.
        model=numpy.zeros((n))
        
        saved_coef_b=12
        saved_coef_r=12
        
        
        for loop in range(20) :

            # compute model

            cal_b=numpy.exp(coef_b*calder_b)
            model_b=cal_b*model
            res_b=flux_b-model_b
            a = invar_b*cal_b**2
            b = invar_b*cal_b*res_b
            chi2_b = numpy.sum(invar_b*res_b**2)

            cal_r=numpy.exp(coef_r*calder_r)
            model_r=cal_r*model
            res_r=flux_r-model_r
            a += invar_r*cal_r**2
            b += invar_r*cal_r*res_r
            chi2_r = numpy.sum(invar_r*res_r**2)
            
            model+=b/a
            model_b=cal_b*model
            model_r=cal_r*model
            
            # compute calibration
            a = numpy.sum(invar_b*(calder_b*model_b)**2)
            b = numpy.sum(invar_b*calder_b*model_b*(flux_b-model_b))
            coef_b+=b/a
            
            a = numpy.sum(invar_r*(calder_r*model_r)**2)
            b = numpy.sum(invar_r*calder_r*model_r*(flux_r-model_r))
            coef_r+=b/a
            
            maxval=2.
            coef_b=max(-maxval,min(maxval,coef_b))
            coef_r=max(-maxval,min(maxval,coef_r))
            
            
            # center coefs, one must have cal_b(wavesplit)*cal_r(wavesplit) = 1
            # coef_b*calder_b_split + coef_r*calder_r_split = 0
            delta=(coef_b*calder_b_split+coef_r*calder_r_split)/(calder_b_split+calder_r_split)
            coef_b-=delta
            coef_r-=delta
            
            #print "iter",loop,"fiber",fiber,"chi2pdf=",(chi2_b+chi2_r)/(2*nwr-2),"coefs",cal_b[-1],cal_r[0]
            
            if abs(coef_b-saved_coef_b)<0.005 and abs(coef_b-saved_coef_b)<0.005 :
                break
            saved_coef_b=coef_b
            saved_coef_r=coef_r
            
            if False :
                pylab.plot(knots,cal_b)
                pylab.plot(knots,cal_r)
                
                pylab.errorbar(knots,flux_b,1/numpy.sqrt(invar_b))
                pylab.errorbar(knots,flux_r,1/numpy.sqrt(invar_r))
                pylab.plot(knots,model,lw=2)
                pylab.show()
            
            

        coef["b%d"%spectro][fiber]=coef_b
        coef["r%d"%spectro][fiber]=coef_r
        print "spec",spectro,"fiber",fiber,"coef=",coef_b,coef_r
        
        # apply to data before merge
        cal_b_original[b_operlap_indices]=numpy.exp(coef_b*calder_b_original[b_operlap_indices])
        cal_r_original[r_operlap_indices]=numpy.exp(coef_r*calder_r_original[b_operlap_indices])
        
        input_flux["b%d"%spectro][fiber,b_operlap_indices]/=cal_b_original[b_operlap_indices]
        input_flux["r%d"%spectro][fiber,r_operlap_indices]/=cal_r_original[r_operlap_indices]
        input_invar["b%d"%spectro][fiber,b_operlap_indices]*=(cal_b_original[b_operlap_indices])**2
        input_invar["r%d"%spectro][fiber,r_operlap_indices]*=(cal_r_original[r_operlap_indices])**2
        
        
    

# merge data        
b_indices=numpy.where(wave_b<=wavesplit)[0]
r_indices=numpy.where(wave_r>wavesplit)[0]
wave = numpy.append(wave_b[b_indices],wave_r[r_indices])
nwave=wave.shape[0]
nfibers1=input_flux["b1"].shape[0]
nfibers2=input_flux["b2"].shape[0]
flux  = numpy.zeros((nfibers1+nfibers2,nwave))
invar = numpy.zeros((nfibers1+nfibers2,nwave))

nb=b_indices.shape[0]
flux[:nfibers1,:nb]=input_flux["b1"][:,b_indices]
invar[:nfibers1,:nb]=input_invar["b1"][:,b_indices]
flux[nfibers1:,:nb]=input_flux["b2"][:,b_indices]
invar[nfibers1:,:nb]=input_invar["b2"][:,b_indices]
indices=r_indices
flux[:nfibers1,nb:]=input_flux["r1"][:,r_indices]
invar[:nfibers1,nb:]=input_invar["r1"][:,r_indices]
flux[nfibers1:,nb:]=input_flux["r2"][:,r_indices]
invar[nfibers1:,nb:]=input_invar["r2"][:,r_indices]

# write file

ofilename = "coadd-merged-%s-%s.fits"%(str(plateid),str(mjd))
output_hdulist=pyfits.HDUList([pyfits.PrimaryHDU(flux),pyfits.ImageHDU(invar,name="IVAR"),\
                                   pyfits.ImageHDU(wave,name="WAVELENGTH")]) #,pyfits.ImageHDU(R,name="RESOLUTION")])

# add some keys
output_header=output_hdulist[0].header

if True :
    ifilename=input_filename["b1"]
    hdulist=pyfits.open(ifilename);
    input_header=hdulist[0].header 
    odico={}
    for k in output_header.keys() :
        odico[k]=output_header[k]
       
    for k in input_header.keys() :
        if not odico.has_key(k) :
            try :                
                output_header.update(k,input_header[k],"from %s"%ifilename)
            except :
                pass
    hdulist.close()

output_hdulist.writeto(ofilename,clobber=True)
print "wrote",ofilename
sys.exit(0)
    
