# SPECtral EXtraction

Package for PSF measurement in fiber-fed spectrograph (BOSS data, DESI simulations).
This package is intended to be used with [HARP](https://github.com/tskisner/HARP) and/or [specter](https://github.com/desihub/specter) extraction codes.
The installation procedure is described in the [INSTALL file](INSTALL.md), the code depends on HARP.


### Usage for DESI 

```
Allowed Options:
  -h [ --help ]                 display usage information
  -a [ --arc ] arg              arc pre-reduced fits image file name 
                                (mandatory), ex:  sdProc-b1-00108382.fits
  --flux-hdu arg (=1)            flux hdu in input arc fits
  --ivar-hdu arg (=2)            ivar hdu in input arc fits
  --mask-hdu arg (=3)            mask hdu in input arc fits
  --header-hdu arg (=1)          header hdu in input arc fits
  --xcoord-file arg             fits image file name with xcoord legendre 
                                polynomial of wavelength (mandatory)
  --xcoord-hdu arg (=1)         hdu of xcoord legendre polynomial of wavelength
  --ycoord-file arg             fits image file name with ycoord legendre 
                                polynomial of wavelength (mandatory)
  --ycoord-hdu arg (=1)         hdu of ycoord legendre polynomial of wavelength
  --first_bundle arg            first fiber bundle to fit
  --last_bundle arg             last fiber bundle to fit
  --first_fiber arg             first fiber (must be in bundle)
  --last_fiber arg              last fiber (must be in bundle)
  --half_size_x arg (=4)        half size of PSF stamp (full size is 
                                2*half_size+1)
  --half_size_y arg (=4)        half size of PSF stamp (full size is 
                                2*half_size+1)
  --psfmodel arg                PSF model, default is GAUSSHERMITE
  --positions                   fit positions of each spot individually after 
                                global fit for debugging
  -v [ --verbose ]              turn on verbose mode
  --lamplines arg               lamp lines ASCII file name (def. is 
                                $SPECEXDATA/opfiles/lamplines.par)
  --core                        dump core files when harp exception is thrown
  --gauss_hermite_deg arg (=3)  degree of Hermite polynomials (same for x and 
                                y, only if GAUSSHERMITE psf)
  --gauss_hermite_deg2 arg (=2) degree of Hermite polynomials (same for x and 
                                y, only if GAUSSHERMITE2 psf)
  --legendre_deg_wave arg (=4)  degree of Legendre polynomials along wavelength
                                (can be reduced if missing data)
  --legendre_deg_x arg (=1)     degree of Legendre polynomials along x_ccd (can
                                be reduced if missing data)
  --trace_deg_wave arg (=0)     degree of Legendre polynomials along wavelength
                                for fit of traces
  --trace_deg_x arg (=0)        degree of Legendre polynomials along x_ccd for 
                                fit of traces
  --psf_error arg (=0)          psf fractional uncertainty (default is 0.01, 
                                for weights in the fit)
  --psf_core_wscale arg         scale up the weight of pixels in 5x5 PSF core
  --fit_psf_tails               unable fit of psf tails
  --fit_continuum               unable fit of continuum
  --no_trace_fit                do not fit traces
  --out_xml arg                  output psf xml file name
  --out_fits arg                 output psf fits file name
  --out_spots arg                output spots file name
  --prior arg                    gaussian prior on a param : 'name' value error
  --tmp_results                  write tmp results
```

### Example for DESI

```
specex_desi_psf_fit -a $DIR/pix-b-arc.fits --xcoord-file psf-b-boot.fits --xcoord-hdu 1 --ycoord-file psf-b-boot.fits --ycoord-hdu 2 --out_xml psf-b0-1-bundle0-tracefit-notail.xml --out_spots psf-b0-1-bundle0-spots-tracefit-notail.xml --out_fits  psf-b0-1-bundle0-spots-tracefit-notail.fits --first_bundle 0 --last_bundle 0 --gauss_hermite_deg 6 --psfmodel GAUSSHERMITE --half_size_x 7 --half_size_y 4 -v --core  --legendre_deg_x 1 --legendre_deg_wave 4 --trace_deg_x 6 --trace_deg_wave 6
```

 * `pix-b-arc.fits` is a DESI CCD image for an arc lamp exposure, flux,ivar,mask are expected to be found in the first 3 HDUs.

 * `psf-b-boot.fits` contains the trace locations, it is the result of the code `desi_bootcalib.py` in [desispec](https://github.com/desihub/desispec)

 * Output PSF is both in xml and fits formats. The fits format is described in [gauss-hermite-psf-datamodel.rst](doc/gauss-hermite-psf-datamodel.rst). It is read by the [specter](https://github.com/desihub/specter) package.


### Parallel processing

specex is not MPI. It uses internally `openmp` to speed up the computation (note also specex runs significantly faster with the configuration option `--optimize=3`). The number of cpu used is given by the value environment variable `OMP_NUM_THREADS`. A convenient way to do parallel computing is to fit the PSF per block of fibers (or bundle, as called in the code), and then merge the results.

Example (in `bash`):

```
for BUNDLE in `seq 0 19`; do
  specex_desi_psf_fit ... --first_bundle $BUNDLE --last_bundle $BUNDLE --out_xml psf-of-bundle-${BUNDLE}.xml &
done
wait
specex_merge_psf psf-of-bundle-*.xml --out-fits psf.fits --out-xml  psf.xml
```





