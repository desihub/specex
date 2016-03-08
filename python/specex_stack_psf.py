#!/usr/bin/env python

import astropy.io.fits as pyfits
import numpy as np
import sys,pylab

if len(sys.argv)<3 :
    print sys.argv[0],"residuals.fits spots.list"
    sys.exit(0)


hdulist=pyfits.open(sys.argv[1])
data=hdulist[0].data
model=hdulist[1].data

vals=np.loadtxt(sys.argv[2]).T
ww=vals[0]
xx=vals[3]
yy=vals[4]
ff=vals[5]

hx=7
hy=7

resampling=16
bins=np.linspace(-hx,hx,(2*hx+1)*resampling)
xbins=bins[:-1]+(bins[1]-bins[0])/2.

data_xprof=np.zeros((bins.size-1))
model_xprof=np.zeros((bins.size-1))

data_yprof=np.zeros((bins.size-1))
model_yprof=np.zeros((bins.size-1))

count=0
for w,x,y,f in zip(ww,xx,yy,ff) :
    
    count += 1
    if count %100 == 0 :
        print count,w,x,y,f
    xi=int(x)
    yi=int(y)
    
    data_stamp_xprof = np.sum(data[yi-2:yi+2+1,xi-hx:xi+hx+1],axis=0)
    model_stamp_xprof = np.sum(model[yi-2:yi+2+1,xi-hx:xi+hx+1],axis=0)
    stamp_x = np.linspace(-hx-(x-xi+0.5),hx+1-(x-xi+0.5),2*hx+1)
        
    data_stamp_yprof = np.sum(data[yi-hy:yi+hy+1,xi-2:xi+2+1],axis=1)
    model_stamp_yprof = np.sum(model[yi-hy:yi+hy+1,xi-2:xi+2+1],axis=1)
    stamp_y = np.linspace(-hy-(y-yi+0.5),hy+1-(y-yi+0.5),2*hy+1)
    
    if stamp_y.size != data_stamp_yprof.size :
        continue

    data_xprof += np.interp(xbins,stamp_x,data_stamp_xprof)
    model_xprof += np.interp(xbins,stamp_x,model_stamp_xprof)
    data_yprof += np.interp(xbins,stamp_y,data_stamp_yprof)
    model_yprof += np.interp(xbins,stamp_y,model_stamp_yprof)

norme=np.max(data_xprof)
a=pylab.subplot(1,2,1)
a.plot(xbins,data_xprof/norme,"-",c="b",lw=2,label="stacked data")
a.plot(xbins,model_xprof/norme,"--",c="r",lw=2,label="stacked PSF model")
a.legend(loc="upper center")
a.set_ylim([-0.05,1.2])
a.set_xlabel("X CCD")

norme=np.max(data_yprof)
a=pylab.subplot(1,2,2)
a.plot(xbins,data_yprof/norme,"-",c="b",lw=2,label="stacked data")
a.plot(xbins,model_yprof/norme,"--",c="r",lw=2,label="stacked PSF model")
a.legend(loc="upper center")
a.set_ylim([-0.05,1.2])
a.set_xlabel("Y CCD")

pylab.show()
sys.exit(0)


pylab.plot(x,y,"o")
pylab.show()

