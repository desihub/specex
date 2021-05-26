import numpy as np
import fitsio
from fitsio import FITS, FITSHDR
import sys
import math

filename1 = sys.argv[1]
filename2 = sys.argv[2]

fits1 = FITS(filename1,'r')
fits2 = FITS(filename2,'r')

xtrace1 = fits1['xtrace'].read()
ytrace1 = fits1['ytrace'].read()
psf1    = fits1['psf'].read()['COEFF']
pvals1  = fits1['psf'].read()['PARAM']
i1      = np.where((pvals1=='STATUS') | (pvals1=='STATUS  '))
stat1   = psf1[i1,:,0].astype(int)
sel1    = np.where(stat1==0)[0]
psf1    = psf1[:,sel1,:]

xtrace2 = fits2['xtrace'].read()
ytrace2 = fits2['ytrace'].read()
psf2    = fits2['psf'].read()['COEFF']
pvals2  = fits2['psf'].read()['PARAM']
i2      = np.where((pvals2=='STATUS') | (pvals2=='STATUS  '))
stat2   = psf2[i2,:,0].astype(int)
sel2    = np.where(stat2==0)[0]
psf2    = psf2[:,sel2,:]

dxtrace = xtrace2 / xtrace1 - 1.0
dytrace = ytrace2 / ytrace1 - 1.0

epsilon=1e-64
psf1[psf1<epsilon] = 0
psf2[psf2<epsilon] = 0

dpsf = (psf2+epsilon) / (psf1+epsilon) - 1.0
num   = psf2 - psf1
denom = 0.5 * (psf1 + psf2)
dpsf  = num / denom
dpsf[num==0] = 0

print('xtrace difference range',dxtrace.min(),dxtrace.max())
print('ytrace difference range',dytrace.min(),dytrace.max())
print('   psf difference range',   dpsf.min(),   dpsf.max())
print('   psf other           ',dpsf.mean(),psf1.mean(),psf2.mean())
