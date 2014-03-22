#!/usr/bin/env python


# based on spflux_v5.pro
# we do not refit the kurucz parameters or redshift
# we simply load the correct kurucz model with correct redshift and apply
# magnitude correction and reddenning

import pyfits,sys,json,pylab,string,numpy,os,scipy,scipy.sparse,scipy.linalg
from scipy.sparse.linalg import spsolve
import scipy.interpolate
from math import *
from specex_air_to_vacuum import *

hc = 1.986438e-8 # in units of (ergs.A)^{-1} (hc=2*pi* hbar*c=2*pi*197.326 eV nm = 6.28318*197.326*1.60218e-12*10)



# from idlutils/pro/dust/ext_ccm.pro
def ccm_dust_extinction(wave,Rv=3.1) :
    
    xx = 10000./wave
    indices_LO  = numpy.where(xx>8.0)[0]                                                 # No data, lambda < 1250 Ang
    indices_FUV = numpy.intersect1d(numpy.where(xx>5.9)[0],numpy.where(xx<=8.0)[0])      # UV + FUV
    indices_NUV = numpy.intersect1d(numpy.where(xx>3.3)[0],numpy.where(xx<=5.9)[0])      # UV + FUV
    indices_OPT = numpy.intersect1d(numpy.where(xx>1.1)[0],numpy.where(xx<=3.3)[0])      # Optical/NIR
    indices_IR  = numpy.intersect1d(numpy.where(xx>0.3)[0],numpy.where(xx<=1.1)[0])      # IR
    indices_HI  = numpy.where(xx<=0.3)[0]                                                # No data, lambda > 33,333 Ang
    
    extinction = numpy.zeros(wave.shape)
    #tmp 
    yy    = numpy.zeros(wave.shape)
    afac  = numpy.zeros(wave.shape)
    bfac  = numpy.zeros(wave.shape)
    
    extinction[indices_LO]=5.0
    
    afac[indices_FUV] = 1.752 - 0.316*xx[indices_FUV] - 0.104 / ( (xx[indices_FUV]-4.67)**2 + 0.341 ) - 0.04473*(xx[indices_FUV]-5.9)**2 - 0.009779*(xx[indices_FUV]-5.9)**3
    bfac[indices_FUV] = -3.090 + 1.825*xx[indices_FUV] + 1.206 / ( (xx[indices_FUV]-4.62)**2 + 0.263 ) + 0.2130*(xx[indices_FUV]-5.9)**2 + 0.1207*(xx[indices_FUV]-5.9)**3
    
    afac[indices_NUV] = 1.752 - 0.316*xx[indices_NUV] - 0.104 / ( (xx[indices_NUV]-4.67)**2 + 0.341 ) 
    bfac[indices_NUV] = -3.090 + 1.825*xx[indices_NUV] + 1.206 / ( (xx[indices_NUV]-4.62)**2 + 0.263 ) 
    
    yy[indices_OPT] = xx[indices_OPT] - 1.82
    afac[indices_OPT] = 1.0 + 0.17699*yy[indices_OPT] \
        - 0.50447*yy[indices_OPT]**2 - 0.02427*yy[indices_OPT]**3 \
        + 0.72085*yy[indices_OPT]**4 + 0.01979*yy[indices_OPT]**5 \
        - 0.77530*yy[indices_OPT]**6 + 0.32999*yy[indices_OPT]**7
    bfac[indices_OPT] = 1.41338*yy[indices_OPT] \
        + 2.28305*yy[indices_OPT]**2 + 1.07233*yy[indices_OPT]**3 \
        - 5.38434*yy[indices_OPT]**4 - 0.62251*yy[indices_OPT]**5 \
        + 5.30260*yy[indices_OPT]**6 - 2.09002*yy[indices_OPT]**7
    
    yy[indices_IR] = xx[indices_IR]**1.61
    afac[indices_IR] = 0.574*yy[indices_IR]
    bfac[indices_IR] = -0.527*yy[indices_IR]
    
    yy[indices_HI] = xx[indices_HI]**1.61
    afac[indices_HI] = 0.574*yy[indices_HI]
    bfac[indices_HI] = -0.527*yy[indices_HI]
    
    extinction = afac + bfac / Rv
    
    return extinction

# from idlutils/pro/dust/ext_odonnell.pro
def odonnell_dust_extinction(wave,Rv=3.1) :
    xx = 10000./wave
    
    indices_optical = numpy.intersect1d(numpy.where(xx>=1.1)[0],numpy.where(xx<=3.3)[0])
    indices_other   = numpy.union1d(numpy.where(xx<1.1)[0],numpy.where(xx>3.3)[0])

    extinction = numpy.zeros(wave.shape)
    # tmp
    yy    = numpy.zeros(wave.shape)
    afac  = numpy.zeros(wave.shape)
    bfac  = numpy.zeros(wave.shape)
    
    yy[indices_optical] = xx[indices_optical] - 1.82

    afac[indices_optical] = 1.0 + 0.104*yy[indices_optical] \
        - 0.609*yy[indices_optical]**2 + 0.701*yy[indices_optical]**3 \
        + 1.137*yy[indices_optical]**4 - 1.718*yy[indices_optical]**5 \
        - 0.827*yy[indices_optical]**6 + 1.647*yy[indices_optical]**7 \
        - 0.505*yy[indices_optical]**8

    bfac[indices_optical] = 1.952*yy[indices_optical] \
        + 2.908*yy[indices_optical]**2 - 3.989*yy[indices_optical]**3 \
        - 7.985*yy[indices_optical]**4 + 11.102*yy[indices_optical]**5 \
        + 5.491*yy[indices_optical]**6 - 10.805*yy[indices_optical]**7 \
        + 3.347*yy[indices_optical]**8
    
    extinction[indices_optical] = afac[indices_optical] + bfac[indices_optical] / Rv
    
    if len(indices_other)>0 :
        extinction[indices_other] = ccm_dust_extinction(wave[indices_other],Rv)
    
    return extinction
    
    

# returns a spectrum based on model parameters for each fiber
# spec_parameters is a dictionnary of of dictionnaries (primary key is the fiber)
# target_wave is an array of wavelength 
def specex_read_kurucz(modelfilename, spec_parameters, target_wave=None) :

    print "reading kurucz spectra in",modelfilename
    
    hdulist=pyfits.open(modelfilename)

    table=hdulist[1].data
    data=hdulist[0].data
    n_model_spec=data.shape[0]
    n_model_wave=data.shape[1]
    
    spec_indices=[]

    # find corresponding models
    spec=0
    for fiber in spec_parameters :
        #print fiber,spec_parameters[fiber]
        tmp=numpy.where(table.field("MODEL")==spec_parameters[fiber]["MODEL"])[0]
        if len(tmp)!=1 :
            print "error in specex_read_kurucz cannot find model",par["MODEL"]
            sys.exit(12)
        spec_indices.append(tmp[0])
        spec+=1
        
    print "kurucz spectra used = ",spec_indices
    
   
    
    # load wavelength array
    # wave number (CRPIX1-1) is given by CRVAL1, wave step is CD1_1 (-1 because in fits standard first pix index=1)
    crpix1=hdulist[0].header['CRPIX1']
    crval1=hdulist[0].header['CRVAL1']
    cd1_1=hdulist[0].header['CD1_1']
    model_wave_step   = cd1_1
    model_wave_offset = (crval1-cd1_1*(crpix1-1))
    model_air_wave=model_wave_step*numpy.arange(n_model_wave) + model_wave_offset
    
    # convert air to vacuum
    model_wave_at_redshift0 = convert_air_to_vacuum(model_air_wave)

    # need to duplicate here for each spectrum
    model_wave=numpy.zeros((len(spec_indices),len(model_wave_at_redshift0)))
    
    # apply redshift (norm doesn't matter, fixed afterwards)
    spec=0
    for fiber in spec_parameters :
        redshift=spec_parameters[fiber]["Z"]
        model_wave[spec] = (1.+redshift)*model_wave_at_redshift0
        spec+=1
    
    # model fluxes 
    model_flux = data[spec_indices]
    
    # convert erg/s/cm2/A to photons/s/cm2/A 
    # number of photons = energy/(h*nu) = energy * wl/(2*pi*hbar*c)
    # 2*pi* hbar*c = 2* pi * 197 eV nm = 6.28318*197.326*1.60218e-12*10 = 1.986438e-8 = 1/5.034135e7 ergs.A    
    # (but we don't care of the norm. anyway here)
    model_flux *= 5.034135e7*model_wave
    
    
    if target_wave == None :
        
        # no rebinning requested, we simply keep data in generous range 2500,12000A
        #waveindexmin=n_model_wave
        #waveindexmax=0
        #for spec in range(model_wave.shape[0]) :
        #    waveindexmin=min(waveindexmin,numpy.where(model_wave[spec]>2500)[0][0]-1)
            #waveindexmax=max(waveindexmax,numpy.where(model_wave[spec]>12000)[0][0])
        #model_wave = model_wave[:,waveindexmin:]
        #model_flux = model_flux[:,waveindexmin:]
        
        hdulist.close() 
        return model_wave,model_flux
    
    n_target_wave=len(target_wave)
    
    # keep only data in wave range of interest before rebinning 
    
    wavemin=target_wave[0]-(target_wave[1]-target_wave[0])
    wavemax=target_wave[-1]+(target_wave[-1]-target_wave[-2])
    waveindexmin=n_model_wave
    waveindexmax=0
    for spec in range(model_wave.shape[0]) :
        waveindexmin=min(waveindexmin,numpy.where(model_wave[spec]>wavemin)[0][0]-1)
        waveindexmax=max(waveindexmax,numpy.where(model_wave[spec]>wavemax)[0][0]+2)
    model_wave = model_wave[:,waveindexmin:waveindexmax]
    model_flux = model_flux[:,waveindexmin:waveindexmax]
    
    #print model_flux.shape
    #print model_wave[0,:20]
    #print model_flux[0,:20]    
    #print target_wave[0],target_wave[-1]
    #print wavemin,wavemax
    #print waveindexmin,waveindexmax,len(model_wave),len(model_air_wave)
    
    
    
    # now integrate models in bins

    # compute bounds of bins (ok for any monotonous binning)    
    wavebinbounds=numpy.zeros((target_wave.shape[0]+1))
    wavebinbounds[1:-1]=(target_wave[:-1]+target_wave[1:])/2
    wavebinbounds[0]=target_wave[0]-0.5*(target_wave[1]-target_wave[0])
    wavebinbounds[-1]=target_wave[-1]+0.5*(target_wave[-1]-target_wave[-2])
    
    # mean in bins, quite slow (few seconds), there must be a faster version of this
    flux=numpy.zeros((len(spec_indices),n_target_wave))
    for s in range(model_wave.shape[0]) :
        i1=numpy.where(model_wave[s]>wavebinbounds[0])[0][0]
        for w in range(n_target_wave):
            i2=numpy.where(model_wave[s]>wavebinbounds[w+1])[0][0]
            flux[s,w]=numpy.mean(model_flux[s,i1:i2])
            #print w,i1,i2,wavebinbounds[w],wavebinbounds[w+1]
            i1=i2

    hdulist.close()
    return target_wave,flux

# read spFluxcalib fits image
def get_kurucz_parameters(spfluxcalibfilename,starfibers) :
    print "reading kurucz parameters in",spfluxcalibfilename
    hdulist=pyfits.open(spfluxcalibfilename)
    table=hdulist[2].data
    columns=hdulist[2].columns
    modelparams={}
    for fiber in starfibers :
        row=numpy.where(table.field("FIBERID")==fiber+1)[0][0]
        #print fiber,row
        dico={}
        for k in columns.names :
            dico[k]=table.field(k)[row]
        modelparams[fiber]=dico
    
    hdulist.close()
    return modelparams

def read_spframe(spframefilename) :
    print "reading spFrame table in",spframefilename
    hdulist=pyfits.open(spframefilename)
    table=hdulist[5].data
    columns=hdulist[5].columns
    keys=columns.names
    dico={}
    fibers=table.field("FIBERID")
    for fiberid in fibers :
        myfiberid=fiberid-1
        row=numpy.where(table.field("FIBERID")==fiberid)[0][0]
        entry={}
        for k in keys :
            entry[k]=table.field(k)[row]
        dico[myfiberid]=entry
    
    hdulist.close()
    return dico

# input fluxes have to be in (photons/cm2/s/A), NOT (ergs/cm2/s/A)
def compute_ab_mags(wave,flux,filter_filenames) :
    
    nbands=len(filter_filenames)
    nwave=wave.shape[0]

    #print nbands,wave.shape,flux.shape
    
    filter_transmissions=numpy.zeros((nbands,nwave))
    
    band=0
    for filename in filter_filenames :
        #print "reading",filename

        fwave=[]
        ftrans=[]

        file=open(filename)
        for line in file.readlines() :
            if line[0]=="#" :
                continue
            line=line.strip()
            if len(line)==0 :
                continue
            tmp=line.split(" ")
            vals=[]
            for t in tmp :
                if len(t) == 0 :
                    continue
                vals.append(string.atof(t))
            fwave.append(vals[0])
            ftrans.append(vals[1])
        file.close()
        
        fwave=numpy.array(fwave)
        ftrans=numpy.array(ftrans)
        #print band,len(wave),len(fwave),len(ftrans)
        filter_transmissions[band]=numpy.interp(wave,fwave,ftrans)
        band+=1
    
    #print filter_transmissions
    
    # compute integrated flux of input spectrum in electrons
    # input spectrum has to be propto photons/cm2/s/A, not ergs
    # because unit of transmission is electron/photon (it includes a quantum efficiency in electron/photon, the rest of the terms are adimentional)
    integrated_flux=numpy.zeros((nbands))
    for b in range(nbands) :
        integrated_flux[b]=numpy.dot(filter_transmissions[b],flux)
    
        
    # compute AB flux in photons/s/cm2/A , NOT ergs, 
    # we have to be precise here because it defines the relation between spec. and imaging magnitudes
    # Fukugita 1996 :
    # Mag_AB = -2.5 log10( flux (ergs/cm2/s/Hz) ) -48.60
    # We want :
    # Mag_AB = -2.5 log10( flux (photons/cm2/s/A) / AB_ref )
    # so,
    # AB_ref = (photons/ergs) c/lambda**2 (in Hz/A) * 10**(-48.6/2.5)
    #        = (photons/ergs) 2.99792458 * 10**(18-48.6/2.5) / lambda**2 (with lambda in A, c from PDG, no harm to be precise)
    # (photons/energy) = lambda/(2*pi*hbar*c) , with 2*pi* hbar*c = 2* pi * 197.326 eV nm = 6.28318*197.326*1.60218e-12*10 = 1.986438e-8 = 1/5.034135e7 ergs.A 
    # AB_ref = (5.034135e7*lambda) * 2.99792458 * 10**(18-48.6/2.5) / lambda**2 (with lambda in A)
    # AB_ref = 5.479558e6/lambda photons/cm2/s/A , with lambda in A

    # because we use hc in other places to convert back to ergs we use a global variable hc = 1.986438e-8 (ergs.A)^{-1}
    # (we use the speed of light only here)
        
    ab_spectrum = 2.99792458 * 10**(18-48.6/2.5)/hc/wave

    # compute integrated AB flux
    integrated_ab_flux=numpy.zeros((nbands))
    for b in range(nbands) :
        integrated_ab_flux[b]=numpy.dot(filter_transmissions[b],ab_spectrum)
    
    # compute magnitudes
    ab_mags=numpy.zeros((nbands))
    for b in range(nbands) :
        if integrated_flux[b]<=0 :
            ab_mags[b]=99.
        else :
            ab_mags[b]=-2.5*log10(integrated_flux[b]/integrated_ab_flux[b])

    return ab_mags


def evaluate_relative_noise(x,y,step) :
    knots=x[0]+step/2+step*numpy.arange(int((x[-1]-x[0])/step))
    func=scipy.interpolate.splrep(x,y,task=-1,t=knots)
    smooth=scipy.interpolate.splev(x,func,der=0)
    diff=y/smooth-1
    rms=numpy.std(diff)
    return rms

"""
def specex_clip(x,y,iw,step,nsig=3,replace=False,verbose=False) :
    knots=x[0]+step/2+step*numpy.arange(int((x[-1]-x[0])/step))
    
    w=iw
    if w==None :
        w=numpy.ones(y.shape)
    
    saved_rms=0
    for loop in range(10) : 
        func=scipy.interpolate.splrep(x,y,w=w,task=-1,t=knots)
        smooth=scipy.interpolate.splev(x,func,der=0)
        diff=y-smooth
        rms=numpy.std(diff[numpy.where(w>0)[0]])
        if verbose :
            print loop,rms
        if(rms==saved_rms) :
            break
        saved_rms=rms
        indices=numpy.where(numpy.absolute(diff)>nsig*rms)[0]
        w[indices]=0
        if(replace) :
            y[indices]=smooth[indices]

    if not iw==None :
        iw=w
"""

# main

# test
if False :
    wave=1000.+10*numpy.arange(2000)
    #xx=10000/wave
    dust_extinction1=odonnell_dust_extinction(wave,3.1)
    dust_extinction2=ccm_dust_extinction(wave,3.1)
    pylab.plot(wave,dust_extinction1)
    pylab.plot(wave,dust_extinction2)
    pylab.show()
    sys.exit(12)


    


if len(sys.argv)<6 :
    print sys.argv[0],"inspec.fits plPlugMapM.par spFluxCalib.fits spFrame.fits outspec.fits (throughput.fits)"
    sys.exit(12);

infilename=sys.argv[1]
plplgmap=sys.argv[2]
spfluxcalibfilename=sys.argv[3]
spframefilename=sys.argv[4]
outfilename=sys.argv[5]
throughputfilename=""
if len(sys.argv)>6 :
    throughputfilename=sys.argv[6]

specexdatadir=""

try :
    specexdatadir=os.environ["SPECEXDATA"]
except :
    print "error:  you need to set the environment variable SPECEXDATA because this code needs various data"
    sys.exit(12)

modelfilename=specexdatadir+"/kurucz_stds_raw_v5.fits"
if not os.path.isfile(modelfilename) :
    print "error, cannot find",modelfilename
    sys.exit(12);
imaging_filters=[]
for band in ["u","g","r","i","z"] :
    filename="%s/sdss_jun2001_%s_atm.dat"%(specexdatadir,band)
    if not os.path.isfile(filename) :
        print "error, cannot find",filename
        sys.exit(12);
    imaging_filters.append(filename)

# get spectrograph id, hardcoded for now, will be read in fits
specid=1

# find  fibers with std stars
starfibers=[]
file=open(plplgmap)
for line in file.readlines() :
    if line.find("PLUGMAPOBJ") != 0 : 
        continue
    
    
    vals=string.split(line," ")
    holetype=vals[8]
    if holetype != "OBJECT" :
        continue
    objType=vals[21]
    if objType != "SPECTROPHOTO_STD" :
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
    starfibers.append(myfiberid)
file.close()

print "std stars fibers (now starting at 0)=",starfibers

hdulist=pyfits.open(infilename)
spectra=hdulist[0].data
invar=hdulist[1].data
wave=hdulist[2].data
Rdata=hdulist[3].data

#print "DEBUG: do only few stars"
#starfibers = starfibers[1:5]

nstarfibers=len(starfibers)
nfibers=Rdata.shape[0]
d=Rdata.shape[1]/2
nwave=Rdata.shape[2]
offsets = range(d,-d-1,-1)

# read kurucz model parameters for each fiber
model_params   = get_kurucz_parameters(spfluxcalibfilename,starfibers)


# read ebv and calibrated ab magnitudes, I know in advance there are 5 bands
ebv = numpy.zeros((nstarfibers))
measured_ab_mags = numpy.zeros((nstarfibers,5))

spframe_params = read_spframe(spframefilename)
for spec in range(nstarfibers) :
    fiber=starfibers[spec]
    ebv[spec]=spframe_params[fiber]["SFD_EBV"]
    
    # MAG is not what we are interested in : 
    # david: "The MAG values in those files are from the FIBER2FLUXes.
    # Specifically, FIBER2FLUX is the SDSS-calibrated (not exactly AB) of a 2 arc sec
    # diameter aperture placed on the SDSS image after convolving
    # it to 2 arc sec seeing.
    # CALIBFLUX gives the AB-calibrated PSF fluxes for point sources.
    # It's actually ~2.1X the fiber fluxes with an AB correction.  This reduces
    # to PSF fluxes for point sources, and carries along extended objects
    # with the same calibrations.
    #
    #measured_ab_mags[spec]=spframe_params[fiber]["MAG"] 
    calibfluxes=spframe_params[fiber]["CALIBFLUX"]
    for b in range(5) :
        measured_ab_mags[spec,b]=22.5-2.5*log10(calibfluxes[b]) # formula from spflux_v5.pro
    

# load model spectra at native resolution to compute ab mags and compare with measured_ab_mags,
unbinned_wave,unbinned_model_fluxes=specex_read_kurucz(modelfilename,model_params,None)

#print unbinned_wave.shape

# apply extinction to models
for spec in range(nstarfibers) :
    unbinned_model_fluxes[spec] *= 10.**(-(3.1 * ebv[spec] / 2.5) * odonnell_dust_extinction(unbinned_wave[spec],3.1) )

# compute model magnitudes to apply scaling
model_ab_mags = numpy.zeros((nstarfibers,5))
for spec in range(nstarfibers) :
    model_ab_mags[spec] = compute_ab_mags(unbinned_wave[spec],unbinned_model_fluxes[spec],imaging_filters)

print "model AB mags in r =",model_ab_mags[:,2]

# compute scale to apply to models
# as in spflux_v5, used only band 2(of id) = 'r' = band 2(of python) ????? to define the norm 
scale_to_apply_to_models = numpy.zeros((nstarfibers))
scale_to_apply_to_models[:] = 10**(-0.4*(measured_ab_mags[:,2]-model_ab_mags[:,2]))
 
print "scale to apply to models =",scale_to_apply_to_models

# now compute model flux at the wavelength of our data
# fluxes are redshifted and resampled but not dust extinguised nor flux calibrated
unused,model_photon_fluxes=specex_read_kurucz(modelfilename,model_params,wave)

# apply extinction to models
for spec in range(nstarfibers) :
    model_photon_fluxes[spec] *= 10.**(-(3.1 * ebv[spec] / 2.5) * odonnell_dust_extinction(wave,3.1) )

# apply pre-computed scale to models
for spec in range(nstarfibers) :
    model_photon_fluxes[spec] *=  scale_to_apply_to_models[spec]

# models are in electrons/s/cm2/A
# data are electrons/A
# so the ratio data/model is in cm2*s



# now we dump the ratio data/model to see how it goes, this is only for monitoring

# first convolve models to data resolution
convolved_model_photon_fluxes=numpy.zeros(model_photon_fluxes.shape)
for spec in range(nstarfibers) :
    fiber=starfibers[spec]
    R=scipy.sparse.dia_matrix((Rdata[fiber],offsets),(nwave,nwave))
    convolved_model_photon_fluxes[spec]=numpy.dot(R.todense(),model_photon_fluxes[spec])

if True :
    pyfits.HDUList([pyfits.PrimaryHDU(convolved_model_photon_fluxes),pyfits.ImageHDU(numpy.zeros(convolved_model_photon_fluxes.shape)),pyfits.ImageHDU(wave)]).writeto("models.fits",clobber=True)

# compute the ratio
# calib_from_phot_to_elec is the conversion factor from photons/cm2/s/A to electrons/A
# its unit is cm2*s
calib_from_phot_to_elec=numpy.zeros(convolved_model_photon_fluxes.shape)
calib_from_phot_to_elec_invar=numpy.zeros(convolved_model_photon_fluxes.shape)
calib_from_phot_to_elec[:]=spectra[starfibers[:]]/convolved_model_photon_fluxes[:]
calib_from_phot_to_elec_invar[:]=invar[starfibers[:]]*(convolved_model_photon_fluxes[:])**2

# now we start the fitting procedure :

# 1) median of calib_from_phot_to_elec   
median_calib_from_phot_to_elec=numpy.median(calib_from_phot_to_elec,axis=0)

# 2.1) iterative clipping with spline fit of calib_from_phot_to_elec 
# due to incorrect line models (either the model, or the redshift, or worse = the resolution)
if True :
    print "iterative clipping"
    wstep=5.
    for spec in range(nstarfibers) :
        mean=numpy.median(calib_from_phot_to_elec[spec])
        knots=wave[0]+wstep/2+wstep*numpy.arange(int((wave[-1]-wave[0])/wstep))
        saved_rms=0
        for loop in range(10) :
            func=scipy.interpolate.splrep(wave,calib_from_phot_to_elec[spec],w=calib_from_phot_to_elec_invar[spec],task=-1,t=knots)
            smooth_calib_from_phot_to_elec=scipy.interpolate.splev(wave,func,der=0)
            diff=calib_from_phot_to_elec[spec]-smooth_calib_from_phot_to_elec
            rms=numpy.std(diff[numpy.where(calib_from_phot_to_elec_invar[spec]>0)[0]])
            print "spec=",spec,"rms/mean=",rms/mean
            if(rms==saved_rms) :
                break
            saved_rms=rms
            calib_from_phot_to_elec_invar[spec,numpy.where(numpy.absolute(diff)>3*rms)[0]]=0

    print "done"
    #median_calib_from_phot_to_elec=numpy.median(smooth_calib_from_phot_to_elec,axis=0)

# 2.2) compute an achromatic offset of the ratio of each star to the median
# this will absorb at zeroth order differences of fiber aperture corrections
# probably due primarily to mis-alignements (themselves due to uncertainties
# in the plate production + atmospheric differential refraction)
# model is calib_from_phot_to_elec = corr*median_calib_from_phot_to_elec
# fit is corr = (sum_i w_i calib_from_phot_to_elec_i median_calib_from_phot_to_elec_i)/(sum_i w_i median_calib_from_phot_to_elec_i**2)
fiber_aperture_correction=numpy.zeros((nstarfibers))
for spec in range(nstarfibers) :
    wmc = numpy.dot(calib_from_phot_to_elec_invar[spec]*calib_from_phot_to_elec[spec],median_calib_from_phot_to_elec)
    wcc = numpy.dot(calib_from_phot_to_elec_invar[spec]*median_calib_from_phot_to_elec,median_calib_from_phot_to_elec)
    fiber_aperture_correction[spec]=wmc/wcc
    
print "fiber aperture corrections     =",fiber_aperture_correction
print "fiber aperture corrections rms =",numpy.std(fiber_aperture_correction)
# 3) apply achromatic offset
for spec in range(nstarfibers) :
    calib_from_phot_to_elec[spec]/=fiber_aperture_correction[spec]
    calib_from_phot_to_elec_invar[spec]*=(fiber_aperture_correction[spec])**2

# 4) refit median
# this is one way to get the calibration we will record it, but also test other methods
# need to do the median only on data with invar>0
median_calib_from_phot_to_elec=numpy.ma.median(numpy.ma.array(calib_from_phot_to_elec,mask=(calib_from_phot_to_elec_invar==0)),axis=0)
#median_calib_from_phot_to_elec=numpy.median(calib_from_phot_to_elec,axis=0)

if False :
    for spec in range(nstarfibers) :
        pylab.plot(wave,calib_from_phot_to_elec[spec],color='r')
        indices=numpy.where(calib_from_phot_to_elec_invar[spec]>0)[0]
        pylab.plot(wave[indices],calib_from_phot_to_elec[spec,indices],color='k')
    pylab.plot(wave,median_calib_from_phot_to_elec, color='b')
    pylab.show()    

# 5) direct fit
# 5.1) apply achromatic offset to models (inverse of applying it to data)
# for spec in range(nstarfibers) :
#     model_photon_fluxes[spec]*=fiber_aperture_correction[spec]

# 5.2) fit
# compare data to  R*calib*model = R*diag(model)*calib
# print "filling A and B"
# A=numpy.matrix(numpy.zeros((nwave,nwave))) # dense because additions of band matrices not implemented
# B=numpy.zeros((1,nwave))
# for spec in range(nstarfibers) :
#     fiber=starfibers[spec]
#     R=scipy.sparse.dia_matrix((Rdata[fiber],offsets),(nwave,nwave)) # resolution matrix
#     M=scipy.sparse.dia_matrix((model_photon_fluxes[spec,:],[0]),(nwave,nwave)) # model in the form of a diagonal matrix
#     Ninv=scipy.sparse.dia_matrix((invar[fiber,:],[0]),(nwave,nwave)) # inverse variance in the form of a diagonal matrix
#     RM=R*M
#     NinvD=invar[fiber,:]*spectra[fiber,:]
#     tmp2=RM.transpose()*Ninv*RM
#     A+=tmp2.todense()
#     B+=RM.transpose().dot(NinvD)

#convert A to sparse for solving
# Ad=d*2
# Aoffsets = range(Ad,-Ad-1,-1)
# Adata =  numpy.zeros((len(Aoffsets),nwave))  
# for i in range(len(Aoffsets)) :
#     diagonal=numpy.diag(A,Aoffsets[i])
#     off=max(0,Aoffsets[i])
#     Adata[i,off:len(diagonal)+off]=diagonal
# As=scipy.sparse.dia_matrix( (Adata,Aoffsets), shape=(nwave,nwave))
# print "done"
# print "solving"
# As = As.tocsr()
# deconvolved_calib_from_phot_to_elec=spsolve(As,B)
# print "done"

# 5.3) reconvolve at resolution of central fiber for display
# Rcentral=scipy.sparse.dia_matrix((Rdata[nfibers/2],offsets),(nwave,nwave))
# convolved_calib_from_phot_to_elec=numpy.dot(Rcentral.toarray(),deconvolved_calib_from_phot_to_elec)


#6) last method, simple fit in convolved space

# 6.1) apply achromatic offset to models (inverse of applying it to data)
for spec in range(nstarfibers) :
    convolved_model_photon_fluxes[spec]*=fiber_aperture_correction[spec]


# 6.2) fit
# compare data to (R*model)*calib 
# iterative fit with clipping

# copy inverse variance of star data to modify it in the clipping
starinvar=numpy.zeros( convolved_model_photon_fluxes.shape)
for spec in range(nstarfibers) :
    fiber=starfibers[spec]
    starinvar[spec]=invar[fiber]



    
# fit
A=numpy.zeros((nwave)) # it's diagonal
B=numpy.zeros((nwave))
for spec in range(nstarfibers) :
    fiber=starfibers[spec]
    A += convolved_model_photon_fluxes[spec]*convolved_model_photon_fluxes[spec]*starinvar[spec]
    B += convolved_model_photon_fluxes[spec]*spectra[fiber]*starinvar[spec]

convolved_calib_from_phot_to_elec=B/A

dchi2=numpy.zeros( convolved_model_photon_fluxes.shape)
ndata=0
for spec in range(nstarfibers) :
    fiber=starfibers[spec]
    dchi2[spec]=starinvar[spec]*(spectra[fiber]-convolved_calib_from_phot_to_elec*convolved_model_photon_fluxes[spec])**2
    ndata+=len(numpy.where(starinvar[spec]>0)[0])

chi2pdf=numpy.sum(dchi2)/(ndata-nwave)
print "chi2pdf=",chi2pdf

for loop in range(50) :

    # collect outliers
    wave_indices=[]
    for spec in range(nstarfibers) :
        indices=numpy.where(dchi2[spec]>4*chi2pdf)[0]
        wave_indices=numpy.union1d(wave_indices,indices)
    wave_indices = wave_indices.astype(int)
    #print wave_indices
    
    print "number of waves with 4 sigma outliers=",len(wave_indices)
    if len(wave_indices)==0 :
        break
    
    # rm largest outlier
    for w in wave_indices :
        spec=numpy.argmax(dchi2[:,w])
        #print "max outlier at w=",w,"is spec",spec
        starinvar[spec,w]=0
        ndata -= 1
    
    # refit only those wavelength
    A[wave_indices]*=0
    B[wave_indices]*=0
    for spec in range(nstarfibers) :
        fiber=starfibers[spec]
        A[wave_indices] += convolved_model_photon_fluxes[spec,wave_indices]*convolved_model_photon_fluxes[spec,wave_indices]*starinvar[spec,wave_indices]
        B[wave_indices] += convolved_model_photon_fluxes[spec,wave_indices]*spectra[fiber,wave_indices]*starinvar[spec,wave_indices]
    
    convolved_calib_from_phot_to_elec[wave_indices]=B[wave_indices]/A[wave_indices]
    
    # recompute chi2pdf
    for spec in range(nstarfibers) :
        fiber=starfibers[spec]
        dchi2[spec,wave_indices]=starinvar[spec,wave_indices]*(spectra[fiber,wave_indices]-convolved_calib_from_phot_to_elec[wave_indices]*convolved_model_photon_fluxes[spec,wave_indices])**2
        
    chi2pdf=numpy.sum(dchi2)/(ndata-nwave)
    print "loop",loop,"chi2pdf=",chi2pdf
    
calib_from_phot_to_elec_variance = 1/A

if throughputfilename != "" :
    # apply a scaling with exposure time and telescope aperture (approximatly) to compute throughput
    telescope_aperture=3.1415*( ((2.5e2)/2)**2 - ((1.5e2)/2)**2 ) # cm2
    exptime=900. #s
    throughput=numpy.zeros((3,median_calib_from_phot_to_elec.shape[0]))
    throughput[0]=(1./(telescope_aperture*exptime))*median_calib_from_phot_to_elec # s cm2 -> 1
    throughput[1]=(1./(telescope_aperture*exptime))*convolved_calib_from_phot_to_elec # s cm2 -> 1
    pyfits.HDUList([pyfits.PrimaryHDU(throughput),pyfits.ImageHDU(numpy.zeros(throughput.shape)),pyfits.ImageHDU(wave)]).writeto(throughputfilename,clobber=True)

print "method 1 rms (median)    = ",evaluate_relative_noise(wave,median_calib_from_phot_to_elec,5)
print "method 2 rms (conv. fit)  = ",evaluate_relative_noise(wave,convolved_calib_from_phot_to_elec,5)

# apply calibration to data
# here I choose the most robust = median
calib_from_phot_to_elec = convolved_calib_from_phot_to_elec


# we want ergs/cm2/s/A not photons/cm2/s/A so we have to multiply back the calibration
# (photons/energy) = lambda/(2*pi*hbar*c) , with 2*pi* hbar*c = 2* pi * 197.326 eV nm = 6.28318*197.326*1.60218e-12*10 = 1.986438e-8 = 1/5.034135e7 ergs.A 
#calib_from_elec_to_phot = 1./calib_from_phot_to_elec  # unit is photons/cm2/s/electron
calib_from_elec_to_ergs = hc/wave/calib_from_phot_to_elec  # unit is ergs/cm2/s/electron

calib_from_elec_to_ergs_variance = (hc/wave)**2*calib_from_phot_to_elec_variance/calib_from_phot_to_elec**4


# add 1.e17 , this is a choice of norm
calib_from_elec_to_ergs *= 1e17

calib_from_elec_to_ergs_variance *= (1e17)**2 

# data : electron/A  -> 10^-17 ergs/cm2/s/A

spectra[:] *= calib_from_elec_to_ergs
invar[:]   /= calib_from_elec_to_ergs**2

# add noise of calibration
invar[:] = 1/( 1/invar[:] + (spectra[:]/calib_from_elec_to_ergs)**2 * calib_from_elec_to_ergs_variance )

print "writing result to",outfilename
hdulist.writeto(outfilename,clobber=True)

