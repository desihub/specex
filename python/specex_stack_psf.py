#!/usr/bin/env python

import astropy.io.fits as pyfits
import numpy as np
import sys,pylab

if len(sys.argv)<3 :
    print sys.argv[0],"residuals.fits spots.list"
    sys.exit(0)


hdulist=pyfits.open(sys.argv[1])
hdulist.info()
data=hdulist[0].data
model=hdulist["MODEL"].data
ivar=(hdulist["PULL"].data/(hdulist["RESIDUAL"].data+(hdulist["RESIDUAL"].data==0)))**2 # (PULL/RES)**2 = IVAR
var=(hdulist["RESIDUAL"].data/(hdulist["PULL"].data+(hdulist["PULL"].data==0)))**2 # VAR

vals=np.loadtxt(sys.argv[2]).T
ww=vals[0]
xx=vals[3]
yy=vals[4]
ff=vals[5]

hx=7
hy=7

resampling=8
bins=np.linspace(-hx-1,hx+1,(2*hx+3)*resampling)
xbins=bins[:-1]+(bins[1]-bins[0])/2.

data_xprof=np.zeros((bins.size-1))
model_xprof=np.zeros((bins.size-1))

data_yprof=np.zeros((bins.size-1))
model_yprof=np.zeros((bins.size-1))

count=0
for w,x,y,f in zip(ww,xx,yy,ff) :
    print w,x,y
    
    count += 1
    if count %100 == 0 :
        print count,w,x,y,f
    xi=int(x)
    yi=int(y)
    
    if 0 : # stack on other dim
        data_stamp_xprof = np.sum(data[yi-2:yi+2+1,xi-hx:xi+hx+1],axis=0)
        var_stamp_xprof = np.sum(var[yi-2:yi+2+1,xi-hx:xi+hx+1],axis=0)
        model_stamp_xprof = np.sum(model[yi-2:yi+2+1,xi-hx:xi+hx+1],axis=0)
        data_stamp_yprof = np.sum(data[yi-hy:yi+hy+1,xi-2:xi+2+1],axis=1)
        var_stamp_yprof = np.sum(var[yi-2:yi+2+1,xi-hx:xi+hx+1],axis=1)        
        model_stamp_yprof = np.sum(model[yi-hy:yi+hy+1,xi-2:xi+2+1],axis=1)
    else : # middle slice
        data_stamp_xprof  = data[yi,xi-hx:xi+hx+1]
        model_stamp_xprof = model[yi,xi-hx:xi+hx+1]
        data_stamp_yprof  = data[yi-hy:yi+hy+1,xi]
        model_stamp_yprof = model[yi-hy:yi+hy+1,xi]
        
    if np.sum(model_stamp_xprof)<=0 :
        continue
    
    stamp_x = np.linspace(-hx-(x-xi+0.5),hx+1-(x-xi+0.5),2*hx+1)
    stamp_y = np.linspace(-hy-(y-yi+0.5),hy+1-(y-yi+0.5),2*hy+1)
    
    if stamp_y.size != data_stamp_yprof.size :
        continue

    
    if 0 : # interpolation version
        data_xprof += np.interp(xbins,stamp_x,data_stamp_xprof)
        model_xprof += np.interp(xbins,stamp_x,model_stamp_xprof)
        data_yprof += np.interp(xbins,stamp_y,data_stamp_yprof)
        model_yprof += np.interp(xbins,stamp_y,model_stamp_yprof)
    else : # stacking version
        
        
        txbins=int(np.floor((xi-x+1.)*resampling))+np.arange(stamp_x.size)*resampling
        tybins=int(np.floor((yi-y+1.)*resampling))+np.arange(stamp_y.size)*resampling
        for i in range(resampling) :
            data_xprof[txbins+i] += data_stamp_xprof
            model_xprof[txbins+i] += model_stamp_xprof
            data_yprof[tybins+i] += data_stamp_yprof
            model_yprof[tybins+i] += model_stamp_yprof


norme=np.max(data_xprof)
a=pylab.subplot(1,2,1)
a.plot(xbins,data_xprof/norme,"-",c="b",lw=2,label="stacked data")
a.plot(xbins,model_xprof/norme,"--",c="r",lw=2,label="stacked PSF model")
a.legend(loc="upper center")
a.set_ylim([-0.05,1.5])
a.set_xlabel("X CCD")

norme=np.max(data_yprof)
a=pylab.subplot(1,2,2)
a.plot(xbins,data_yprof/norme,"-",c="b",lw=2,label="stacked data")
a.plot(xbins,model_yprof/norme,"--",c="r",lw=2,label="stacked PSF model")
a.legend(loc="upper center")
a.set_ylim([-0.05,1.5])
a.set_xlabel("Y CCD")

pylab.show()
sys.exit(0)


pylab.plot(x,y,"o")
pylab.show()

