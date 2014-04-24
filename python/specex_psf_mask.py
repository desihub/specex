#!/usr/bin/env python


import pyfits,sys,string,numpy

if len(sys.argv)!=4 :
    print sys.argv[0],"spX-input.fits psf.fits spX-output.fits"
    print "example: specex_psf_mask.py spXv-b1-00149384.fits /clusterfs/riemann/boss/nextgen/testdata/jguy/psf/PROD006/3647/56219/psf-b1-00149383-flexed-to-00149384.fits  spXvm-b1-00149384.fits"
    sys.exit(0)

nfibers_per_bundle=20 # hardcoded
chi2max=20 # beyond that it's very bad, we must mask out the extracted data

ifilename=sys.argv[1]
psffilename=sys.argv[2]
ofilename=sys.argv[3]

ihdulist=pyfits.open(ifilename)
nfibers=ihdulist[0].data.shape[0]
nbundles=nfibers/nfibers_per_bundle
#psffile=ihdulist[0].header["IN_PSF"]
#print psffile

print "nbundles=",nbundles

# read psf chi2
chi2=numpy.zeros((nbundles))
psfhdulist=pyfits.open(psffilename)
header=psfhdulist[1].header
#print header
for b in range(nbundles) :
    key="B%02dCCHI2"%b
    chi2[b]=header[key]
psfhdulist.close()
print chi2

mask=numpy.zeros((nfibers)).astype(int)

for b in range(nbundles) :
    if chi2[b]>chi2max :
        print "masking fibers",range(b*nfibers_per_bundle,(b+1)*nfibers_per_bundle)
        for fiber in range(b*nfibers_per_bundle,(b+1)*nfibers_per_bundle) :
            mask[fiber]=1
            #ihdulist[0].data[fiber]=0
            ihdulist[1].data[fiber]=0


ihdulist.append(pyfits.ImageHDU(mask,name="FMASK"))
ihdulist.writeto(ofilename)
