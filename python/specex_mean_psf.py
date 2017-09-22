#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import argparse
import sys,os
from numpy.polynomial.legendre import legval,legfit

from desiutil.log import get_logger

def compatible(head1,head2) :
    log=get_logger()
    for k in ["PSFTYPE","NPIX_X","NPIX_Y","HSIZEX","HSIZEY","FIBERMAX","FIBERMIN","FIBERMAX","NPARAMS","LEGDEG","GHDEGX","GHDEGY"] :
        if (head1[k] != head2[k]) :
            log.warning("different %s : %s , %s"%(k,head1[k],head2[k]))
            return False
    return True

def main() :
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input', type = str, default = None, required=True, nargs='*', action='append', 
                    help = 'input specex psf (fits format)')
    parser.add_argument('-o','--output', type = str, default = None, required=True,
		                        help = 'output specex psf (fits format)')
    

    args = parser.parse_args()
    log=get_logger()
    
    inputs=[]
    for i in args.input :
        for j in i :
            inputs.append(j)
    npsf=len(inputs)
    log.info("Will compute the average of %d PSFs"%npsf)
    
    refhead=None
    tables=[]
    xtrace=[]
    ytrace=[]
    hdulist=None
    bundle_rchi2=[]
    nbundles=None
    nfibers_per_bundle=None
    for input in inputs :
        psf=pyfits.open(input)
        if refhead is None :
            hdulist=psf
            refhead=psf["PSF"].header            
            nfibers=(psf["PSF"].header["FIBERMAX"]-psf["PSF"].header["FIBERMIN"])+1
            PSFVER=int(refhead["PSFVER"])
            if(PSFVER<3) :
                print("ERROR NEED PSFVER>=3")
                sys.exit(1)
            
        else :
            if not compatible(psf["PSF"].header,refhead) :
                log.error("psfs %s and %s are not compatible"%(inputs[0],input))
                sys.exit(12)
        tables.append(psf["PSF"].data)
        if "XTRACE" in psf :
            xtrace.append(psf["XTRACE"].data)
        if "YTRACE" in psf :
            ytrace.append(psf["YTRACE"].data)

        rchi2=[]
        b=0
        while "B%02dRCHI2"%b in psf["PSF"].header :
            rchi2.append(psf["PSF"].header["B%02dRCHI2"%b])
            b += 1
        rchi2=np.array(rchi2)
        nbundles=rchi2.size
        bundle_rchi2.append(rchi2)
    
    bundle_rchi2=np.array(bundle_rchi2)
    log.info("bundle_rchi2= %s"%str(bundle_rchi2))
    median_bundle_rchi2 = np.median(bundle_rchi2)
    rchi2_threshold=median_bundle_rchi2+1.
    log.info("median chi2=%f threshold=%f"%(median_bundle_rchi2,rchi2_threshold))
    
    WAVEMIN=refhead["WAVEMIN"]
    WAVEMAX=refhead["WAVEMAX"]
    FIBERMIN=int(refhead["FIBERMIN"])
    FIBERMAX=int(refhead["FIBERMAX"])
    
    
    fibers_in_bundle={}
    
    if PSFVER>=3 :
        i=np.where(tables[0]["PARAM"]=="BUNDLE")[0][0]
        bundle_of_fibers=tables[0]["COEFF"][i][:,0].astype(int)
        bundles=np.unique(bundle_of_fibers)
        for b in bundles :
            fibers_in_bundle[b]=np.where(bundle_of_fibers==b)[0]
    else :
        bmin=refhead["BUNDLMIN"]
        bmax=refhead["BUNDLMAX"]
        nfiber_per_bundle=(bmax-bmin+1)/(FIBERMAX-FIBERMIN+1)
        for b in range(bmin,bmax+1) :
            fibers_in_bundle[b]=np.arange(nfiber_per_bundle*b,nfiber_per_bundle*(b+1))
    
    for b in bundles :
        print("%d : %s"%(b,fibers_in_bundle[b]))
        
    for entry in range(tables[0].size) :
        PARAM=tables[0][entry]["PARAM"]
        log.info("Averaging '%s' coefficients"%PARAM)        
        coeff=[tables[0][entry]["COEFF"]]
        npar=coeff[0][1].size
        for p in range(1,npsf) :
            coeff.append(tables[p][entry]["COEFF"])
        coeff=np.array(coeff)
        
        output_rchi2=np.zeros((bundle_rchi2.shape[1]))
        output_coeff=np.zeros(tables[0][entry]["COEFF"].shape)
        
        #log.info("input coeff.shape  = %d"%coeff.shape)
        #log.info("output coeff.shape = %d"%output_coeff.shape)
        
        # now merge, using rchi2 as selection score
        
        for bundle in fibers_in_bundle.keys() :
            
            ok=np.where(bundle_rchi2[:,bundle]<rchi2_threshold)[0]
            #ok=np.array([0,1]) # debug
            if entry==0 :
                log.info("for fiber bundle %d, %d valid PSFs"%(bundle,ok.size))
            
            
            if ok.size>=2 : # use median
                log.info("bundle #%d : use median"%bundle)
                for f in fibers_in_bundle[bundle]  :
                    output_coeff[f]=np.median(coeff[ok,f],axis=0)
                output_rchi2[bundle]=np.median(bundle_rchi2[ok,bundle])
            elif ok.size==1 : # copy
                log.info("bundle #%d : use only one psf "%bundle)
                for f in fibers_in_bundle[bundle]  :
                    output_coeff[f]=coeff[ok[0],f]
                output_rchi2[bundle]=bundle_rchi2[ok[0],bundle]
                    
            else : # we have a problem here, take the smallest rchi2
                log.info("bundle #%d : take smallest chi2 "%bundle)
                i=np.argmin(bundle_rchi2[:,bundle])
                for f in fibers_in_bundle[bundle]  :
                    output_coeff[f]=coeff[i,f]
                output_rchi2[bundle]=bundle_rchi2[i,bundle]
        
        # now copy this in output table
        hdulist["PSF"].data["COEFF"][entry]=output_coeff
        # change bundle chi2
        for bundle in range(output_rchi2.size) :
            hdulist["PSF"].header["B%02dRCHI2"%bundle]=output_rchi2[bundle]

        if len(xtrace)>0 :
            hdulist["xtrace"].data = np.median(np.array(xtrace),axis=0)
            hdulist["ytrace"].data = np.median(np.array(ytrace),axis=0)
            


        # alter other keys in header
        hdulist["PSF"].header["EXPID"]=0. # it's a mix , need to add the expids here
        
        


    # save output PSF
    hdulist.writeto(args.output,clobber=True)
    log.info("wrote %s"%args.output)
    
        


if __name__ == '__main__':
    main()
    
