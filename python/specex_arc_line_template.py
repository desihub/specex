#!/usr/bin/env python

# use measured lines on testbench to create a DESI arc line template

import os,sys
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import desiutil.io

cameras=["b1","r1","z1"]
scales=[3.,1.,2.] # because of the VPH grating orientation errors
globalscale = 0.2 # to match current sims

waves=np.array([])
fluxes=np.array([])

for camera,scale in zip(cameras,scales) :
    filename="%s/20170118/psf-%s.fits"%(os.environ["DESI_CCD_CALIBRATION_DATA"],camera)
    
    h=pyfits.open(filename)
    t=h["SPOTS"].data
    wave=np.unique(t["WAVE"])
    flux=np.zeros(wave.size)
    for i,w in enumerate(wave) :
        flux[i]=np.mean(t["FLUX"][(t["WAVE"]==w)&(t["FIBER"]>=10)])
    flux *= globalscale*scale
    plt.plot(wave,flux,"o")
    waves=np.append(waves,wave)
    fluxes=np.append(fluxes,flux)

# we take the mean on overlapping cameras to make things simpler
wave=np.unique(waves)
flux=np.zeros(wave.size)
for i,w in enumerate(wave) :
    flux[i]=np.mean(fluxes[waves==w])
plt.plot(wave,flux,"+",c="r")

data={}
data["VACUUM_WAVE"]=wave
data["ELECTRONS"]=flux
data = desiutil.io.encode_table(data)
hdus = pyfits.HDUList([pyfits.PrimaryHDU()])
hdus.append(pyfits.convenience.table_to_hdu(data))
hdus[0].header.add_comment("generated using %s"%os.path.basename(sys.argv[0]))
hdus[0].header.add_comment("and arc lamp images from WINLIGHT in 2017/01/18")
ofilename="arc-lines-average-in-vacuum-from-winlight-20170118.fits"
hdus.writeto(ofilename,overwrite=True)
print("wrote",ofilename)

infile = os.getenv('DESI_ROOT')+'/spectro/templates/calib/v0.3/arc-lines-average-in-vacuum.fits'
h=pyfits.open(infile)
plt.plot(h[1].data["VACUUM_WAVE"],h[1].data["ELECTRONS"],"x",c="k")
plt.show()



