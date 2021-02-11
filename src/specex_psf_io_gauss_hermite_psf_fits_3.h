#ifndef SPECEX_PSF_IO_GH3__H

void write_gauss_hermite_psf_fits_version_3(specex::GaussHermitePSF& psf, fitsfile* fp);
void read_gauss_hermite_psf_fits_version_3(specex::PSF_p& psf, fitsfile* fp, int hdu);

#endif
