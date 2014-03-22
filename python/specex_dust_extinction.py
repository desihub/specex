import numpy


# from idlutils/pro/dust/ext_ccm.pro
def ccm_dust_extinction(wave,Rv=3.1) :
    
    xx = 10000./wave
    indices_LO  = numpy.where(xx>8.0)[0]                                                 # No data, lambda < 1250 Ang
    indices_FUV = numpy.intersect1d(numpy.where(xx>5.9)[0],numpy.where(xx<=8.0)[0])      # UV + FUV
    indices_NUV = numpy.intersect1d(numpy.where(xx>3.3)[0],numpy.where(xx<=5.9)[0])      # UV + FUV
    indices_OPT = numpy.intersect1d(numpy.where(xx>1.1)[0],numpy.where(xx<=3.3)[0])      # Optical/NIR
    indices_IR  = numpy.intersect1d(numpy.where(xx>0.3)[0],numpy.where(xx<=1.1)[0])      # IR
    indices_HI  = numpy.where(xx<=0.3)[0]                                                # No data, lambda > 33,333 Ang
    
    extinction = numpy.zeros(wave.shape)
    #tmp 
    yy    = numpy.zeros(wave.shape)
    afac  = numpy.zeros(wave.shape)
    bfac  = numpy.zeros(wave.shape)
    
    extinction[indices_LO]=5.0
    
    afac[indices_FUV] = 1.752 - 0.316*xx[indices_FUV] - 0.104 / ( (xx[indices_FUV]-4.67)**2 + 0.341 ) - 0.04473*(xx[indices_FUV]-5.9)**2 - 0.009779*(xx[indices_FUV]-5.9)**3
    bfac[indices_FUV] = -3.090 + 1.825*xx[indices_FUV] + 1.206 / ( (xx[indices_FUV]-4.62)**2 + 0.263 ) + 0.2130*(xx[indices_FUV]-5.9)**2 + 0.1207*(xx[indices_FUV]-5.9)**3
    
    afac[indices_NUV] = 1.752 - 0.316*xx[indices_NUV] - 0.104 / ( (xx[indices_NUV]-4.67)**2 + 0.341 ) 
    bfac[indices_NUV] = -3.090 + 1.825*xx[indices_NUV] + 1.206 / ( (xx[indices_NUV]-4.62)**2 + 0.263 ) 
    
    yy[indices_OPT] = xx[indices_OPT] - 1.82
    afac[indices_OPT] = 1.0 + 0.17699*yy[indices_OPT] \
        - 0.50447*yy[indices_OPT]**2 - 0.02427*yy[indices_OPT]**3 \
        + 0.72085*yy[indices_OPT]**4 + 0.01979*yy[indices_OPT]**5 \
        - 0.77530*yy[indices_OPT]**6 + 0.32999*yy[indices_OPT]**7
    bfac[indices_OPT] = 1.41338*yy[indices_OPT] \
        + 2.28305*yy[indices_OPT]**2 + 1.07233*yy[indices_OPT]**3 \
        - 5.38434*yy[indices_OPT]**4 - 0.62251*yy[indices_OPT]**5 \
        + 5.30260*yy[indices_OPT]**6 - 2.09002*yy[indices_OPT]**7
    
    yy[indices_IR] = xx[indices_IR]**1.61
    afac[indices_IR] = 0.574*yy[indices_IR]
    bfac[indices_IR] = -0.527*yy[indices_IR]
    
    yy[indices_HI] = xx[indices_HI]**1.61
    afac[indices_HI] = 0.574*yy[indices_HI]
    bfac[indices_HI] = -0.527*yy[indices_HI]
    
    extinction = afac + bfac / Rv
    
    return extinction

# from idlutils/pro/dust/ext_odonnell.pro
def odonnell_dust_extinction(wave,Rv=3.1) :
    xx = 10000./wave
    
    indices_optical = numpy.intersect1d(numpy.where(xx>=1.1)[0],numpy.where(xx<=3.3)[0])
    indices_other   = numpy.union1d(numpy.where(xx<1.1)[0],numpy.where(xx>3.3)[0])

    extinction = numpy.zeros(wave.shape)
    # tmp
    yy    = numpy.zeros(wave.shape)
    afac  = numpy.zeros(wave.shape)
    bfac  = numpy.zeros(wave.shape)
    
    yy[indices_optical] = xx[indices_optical] - 1.82

    afac[indices_optical] = 1.0 + 0.104*yy[indices_optical] \
        - 0.609*yy[indices_optical]**2 + 0.701*yy[indices_optical]**3 \
        + 1.137*yy[indices_optical]**4 - 1.718*yy[indices_optical]**5 \
        - 0.827*yy[indices_optical]**6 + 1.647*yy[indices_optical]**7 \
        - 0.505*yy[indices_optical]**8

    bfac[indices_optical] = 1.952*yy[indices_optical] \
        + 2.908*yy[indices_optical]**2 - 3.989*yy[indices_optical]**3 \
        - 7.985*yy[indices_optical]**4 + 11.102*yy[indices_optical]**5 \
        + 5.491*yy[indices_optical]**6 - 10.805*yy[indices_optical]**7 \
        + 3.347*yy[indices_optical]**8
    
    extinction[indices_optical] = afac[indices_optical] + bfac[indices_optical] / Rv
    
    if len(indices_other)>0 :
        extinction[indices_other] = ccm_dust_extinction(wave[indices_other],Rv)
    
    return extinction
    
    
