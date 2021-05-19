import numpy as np
import fitsio
from fitsio import FITS, FITSHDR

filename1 = 'fit-psf-b1-00068217-mst.fits'
filename2 = 'fit-psf-b1-00068217-dev.fits'

fits1 = FITS(filename1,'r')
fits2 = FITS(filename2,'r')

xtrace1 = fits1['xtrace'].read()
ytrace1 = fits1['ytrace'].read()
psf1    = fits1['psf'].read()['COEFF']

xtrace2 = fits2['xtrace'].read()
ytrace2 = fits2['ytrace'].read()
psf2    = fits2['psf'].read()['COEFF']

dxtrace = xtrace2 / xtrace1 - 1.0
dytrace = ytrace2 / ytrace1 - 1.0

epsilon=1e-16
dpsf = (psf2+epsilon) / (psf1+epsilon) - 1.0
num   = psf2 - psf1
denom = psf1
dpsf  = num / denom
dpsf[num==0] = 0

print('xtrace difference range',dxtrace.min(),dxtrace.max())
print('ytrace difference range',dytrace.min(),dytrace.max())
print('   psf difference range',   dpsf.min(),   dpsf.max())#dpsf.mean(),np.shape(dpsf),psf1.mean(),psf2.mean())

