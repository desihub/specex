#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import argparse
import sys,os
from numpy.polynomial.legendre import legval,legfit

from desispec.log import get_logger

def compatible(head1,head2) :
    log=get_logger()
    for k in ["PSFTYPE","NPIX_X","NPIX_Y","HSIZEX","HSIZEY","BUNDLMIN","BUNDLMAX","FIBERMAX","FIBERMIN","FIBERMAX","NPARAMS","LEGDEG","GHDEGX","GHDEGY"] :
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
    hdulist=None
    bundle_rchi2=[]
    nbundles=None
    nfibers_per_bundle=None
    for input in inputs :
        psf=pyfits.open(input)
        if refhead is None :
            hdulist=psf
            refhead=psf[1].header            
            nbundles=(psf[1].header["BUNDLMAX"]-psf[1].header["BUNDLMIN"])+1
            nfibers=(psf[1].header["FIBERMAX"]-psf[1].header["FIBERMIN"])+1
            nfibers_per_bundle=nfibers/nbundles
            print "nbundles=%d"%nbundles
            print "nfibers_per_bundle=%d"%nfibers_per_bundle
            
        else :
            if not compatible(psf[1].header,refhead) :
                log.error("psfs %s and %s are not compatible"%(inputs[0],input))
                sys.exit(12)
        tables.append(psf[1].data)

        
        
        
        
        rchi2=np.zeros((nbundles))
        for b in range(nbundles) :
            rchi2[b]=psf[1].header["B%02dRCHI2"%b]
        bundle_rchi2.append(rchi2)
    
    bundle_rchi2=np.array(bundle_rchi2)
    print "bundle_rchi2=",bundle_rchi2

    for entry in range(tables[0].size) :
        PARAM=tables[0][entry]["PARAM"]
        log.info("Averaging '%s' coefficients"%PARAM)
        # check WAVEMIN WAVEMAX compatibility
        WAVEMIN=tables[0][entry]["WAVEMIN"]
        WAVEMAX=tables[0][entry]["WAVEMAX"]
        """
        for p in range(1,npsf) :
            if tables[p][entry]["WAVEMIN"] != WAVEMIN :
                log.error("WAVEMIN not compatible for param %s : %f!=%f"%(PARAM,tables[p][entry]["WAVEMIN"],WAVEMIN)) 
                sys.exit(12)
            if tables[p][entry]["WAVEMAX"] != WAVEMAX :
                log.error("WAVEMAX not compatible for param %s : %f!=%f"%(PARAM,tables[p][entry]["WAVEMAX"],WAVEMAX))
                sys.exit(12)
        """
        # will need to readdress coefs ...         
        coeff=[tables[0][entry]["COEFF"]]
        npar=coeff[0][1].size
        for p in range(1,npsf) :
            if tables[p][entry]["WAVEMIN"] == WAVEMIN and tables[p][entry]["WAVEMAX"] == WAVEMAX :
                coeff.append(tables[p][entry]["COEFF"])
            else :
                icoeff=tables[p][entry]["COEFF"]
                ocoeff=np.zeros(icoeff.shape)
                # need to reshape legpol
                iu=np.linspace(-1,1,npar+3)
                iwavemin=tables[p][entry]["WAVEMIN"]
                iwavemax=tables[p][entry]["WAVEMAX"]
                wave=(iu+1.)/2.*(iwavemax-iwavemin)+iwavemin
                ou=(wave-WAVEMIN)/(WAVEMAX-WAVEMIN)*2.-1.                
                for f in range(icoeff.shape[0]) :
                    val=legval(iu,icoeff[f])
                    ocoeff[f]=legfit(ou,val,deg=npar-1)
                #print ""
                #print icoeff[2]
                #print ocoeff[2]
                coeff.append(ocoeff)
        coeff=np.array(coeff)
        
        output_rchi2=np.zeros((bundle_rchi2.shape[1]))
        output_coeff=np.zeros(tables[0][entry]["COEFF"].shape)
        
        #log.info("input coeff.shape  = %d"%coeff.shape)
        #log.info("output coeff.shape = %d"%output_coeff.shape)
        
        # now merge, using rchi2 as selection score
        rchi2_threshold=2.
        for bundle in range(bundle_rchi2.shape[1]) :
            
            ok=np.where(bundle_rchi2[:,bundle]<rchi2_threshold)[0]
            #ok=np.array([0,1]) # debug
            if entry==0 :
                log.info("for fiber bundle %d, %d valid PSFs"%(bundle,ok.size))
            
            fibers=np.arange(bundle*nfibers_per_bundle,(bundle+1)*nfibers_per_bundle)
            if ok.size>=2 : # use mean (median is dangerous)
                for f in fibers :
                    output_coeff[f]=np.median(coeff[ok,f],axis=0)
                output_rchi2[bundle]=np.median(bundle_rchi2[ok,bundle])
            elif ok.size==1 : # copy
                for f in fibers :
                    output_coeff[f]=coeff[ok[0],f]
                output_rchi2[bundle]=bundle_rchi2[ok[0],bundle]
                    
            else : # we have a problem here, take the smallest rchi2
                i=np.argmin(bundle_rchi2[:,bundle])
                for f in fibers :
                    output_coeff[f]=coeff[i,f]
                output_rchi2[bundle]=bundle_rchi2[i,bundle]
        
        # now copy this in output table
        hdulist[1].data["COEFF"][entry]=output_coeff
        # change bundle chi2
        for bundle in range(output_rchi2.size) :
            hdulist[1].header["B%02dRCHI2"%bundle]=output_rchi2[bundle]
        
        # alter other keys in header
        hdulist[1].header["EXPID"]=0. # it's a mix , need to add the expids here
        
        


    # save output PSF
    hdulist.writeto(args.output,clobber=True)
    log.info("wrote %s"%args.output)
    
        


if __name__ == '__main__':
    main()
    
