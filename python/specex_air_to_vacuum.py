import numpy

# convert air to vacuum, this is IDL routine airtovac for instance :
# http://idlastro.gsfc.nasa.gov/ftp/pro/astro/airtovac.pro
def convert_air_to_vacuum(air_wave) :
    # idl code :
    # for iter=0, 1 do begin
    # sigma2 = (1d4/double(wave_vac[g]) )^2.     ;Convert to wavenumber squared
    # ; Compute conversion factor
    # fact = 1.D +  5.792105D-2/(238.0185D0 - sigma2) + $
    #                        1.67917D-3/( 57.362D0 - sigma2)
    # wave_vac[g] = wave_air[g]*fact              ;Convert Wavelength
    # endfor
        
    sigma2 = (1e4/air_wave)**2
    fact = 1. +  5.792105e-2/(238.0185 - sigma2) +  1.67917e-3/( 57.362 - sigma2)
    vacuum_wave = air_wave*fact

    # comparison with http://www.sdss.org/dr7/products/spectra/vacwavelength.html
    # where : AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4)
    # air_wave=numpy.array([4861.363,4958.911,5006.843,6548.05,6562.801,6583.45,6716.44,6730.82])
    # expected_vacuum_wave=numpy.array([4862.721,4960.295,5008.239,6549.86,6564.614,6585.27,6718.29,6732.68])
    # test ok
    return vacuum_wave
