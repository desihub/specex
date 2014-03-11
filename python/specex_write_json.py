#!/usr/bin/env python


import pyfits,sys,json

if len(sys.argv)<5 :
    print sys.argv[0],"spec.fits image.fits psf.xml extract.json"
    sys.exit(12);

spec_filename=sys.argv[1]
image_filename=sys.argv[2]
psf_filename=sys.argv[3]
json_filename=sys.argv[4]

hdulist=pyfits.open(spec_filename)

(nspec,nwave) = hdulist[0].data.shape
wave =  hdulist[2].data

if len(wave) != nwave :
    print "error"
    sys.exit(12)

wavemin=wave[0]
wavemax=wave[-1]
wavebin=wave[1]-wave[0]
fibermin=0 # I assume this
fibermax=nspec-1
print nspec,nwave,wavemin,wavemax,wavebin


all={}

all["psf_type"]="specex"
psf={}
psf["type"]="xml"
psf["path"]=psf_filename
psf["wavebin"]=wavebin
psf["wavemin"]=wavemin
psf["wavemax"]=wavemax
psf["fibermin"]=fibermin
psf["fibermax"]=fibermax
all["psf"]=psf

all["image_type"]="fits"
image={}
image["path"]=image_filename
image["signal"]=1
image["noise"]=2
all["image"]=image

ofile=open(json_filename,"w")
ofile.write(json.dumps(all, sort_keys=True,indent=4, separators=(',', ': ')))
ofile.close()

print "wrote",json_filename
