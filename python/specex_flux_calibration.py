#!/usr/bin/env python

# convert air to vacuum, this is IDL routine airtovac for instance :
# http://idlastro.gsfc.nasa.gov/ftp/pro/astro/airtovac.pro
def convert_air_to_vacuum(air_wave) :
    # idl code :
    # for iter=0, 1 do begin
    # sigma2 = (1d4/double(wave_vac[g]) )^2.     ;Convert to wavenumber squared
    # ; Compute conversion factor
    # fact = 1.D +  5.792105D-2/(238.0185D0 - sigma2) + $
    #                        1.67917D-3/( 57.362D0 - sigma2)
    # wave_vac[g] = wave_air[g]*fact              ;Convert Wavelength
    # endfor
    
    
    sigma2 = (1e4/air_wave)**2
    fact = 1. +  5.792105e-2/(238.0185 - sigma2) +  1.67917e-3/( 57.362 - sigma2)
    vacuum_wave = air_wave*fact

    # comparison with http://www.sdss.org/dr7/products/spectra/vacwavelength.html
    # where : AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4)
    # air_wave=numpy.array([4861.363,4958.911,5006.843,6548.05,6562.801,6583.45,6716.44,6730.82])
    # expected_vacuum_wave=numpy.array([4862.721,4960.295,5008.239,6549.86,6564.614,6585.27,6718.29,6732.68])
    # test ok
    return vacuum_wave


# returns a spectrum based on model parameters
def specex_read_kurucz(modelfilename, spec_parameters, target_wave) :

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
    
    n_target_wave=len(target_wave)
    
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
    
    # apply redshift (I should in principle convert energy flux as well here, but I think it's going to be hard to match BOSS norm anyway)
    spec=0
    for fiber in spec_parameters :
        redshift=spec_parameters[fiber]["Z"]
        model_wave[spec] = (1.+redshift)*model_wave_at_redshift0
        spec+=1
    
    # model fluxes 
    model_flux = data[spec_indices]
    
    if True :
        # keep only data in wave range of interest    
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
    
    # convert erg/s/A to photons/s/A I NEED TO CHECK THE ORIGINAL UNIT !!!!!!!!!
    # number of photons = energy/(h*nu) = energy * wl/(2*pi*hbar*c)
    # 2*pi* hbar*c = 2* pi * 197 eV nm = 6.28318*197.326*1.60218e-12*10 = 1.986438e-8 = 1/5.034135e7 ergs.A    
    model_flux *= 5.034135e7*model_wave
    
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
    
    return flux


def get_kurucz_parameters(spfluxcalibfilename,starfibers) :
    print "reading kurucz parameters in",spfluxcalibfilename
    calibhdulist=pyfits.open(spfluxcalibfilename)
    table=calibhdulist[2].data
    columns=calibhdulist[2].columns
    modelparams={}
    for fiber in starfibers :
        row=numpy.where(table.field("FIBERID")==fiber+1)[0][0]
        #print fiber,row
        dico={}
        for k in columns.names :
            dico[k]=table.field(k)[row]
        modelparams[fiber]=dico
    return modelparams

import pyfits,sys,json,pylab,string,numpy,os,scipy,scipy.sparse,scipy.linalg
from scipy.sparse.linalg import spsolve
from math import *

if len(sys.argv)<6 :
    print sys.argv[0],"inspec.fits plPlugMapM.par kurucz_stds_raw_v5.fits spFluxCalib.fits outspec.fits"
    sys.exit(12);

infilename=sys.argv[1]
plplgmap=sys.argv[2]
modelfilename=sys.argv[3]
spfluxcalibfilename=sys.argv[4]
outfilename=sys.argv[5]

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


nstarfibers=len(starfibers)
nfibers=Rdata.shape[0]
d=Rdata.shape[1]/2
nwave=Rdata.shape[2]
offsets = range(d,-d-1,-1)

# get kurucz model parameters for each fiber
modelparams = get_kurucz_parameters(spfluxcalibfilename,starfibers)

models=specex_read_kurucz(modelfilename,modelparams,wave)

if os.path.isfile("models.fits") :
    os.unlink("models.fits")
pyfits.HDUList([pyfits.PrimaryHDU(models),pyfits.ImageHDU(models),pyfits.ImageHDU(wave)]).writeto("models.fits")

sys.exit(0)

# read kurucz models



# sys.exit(0)
