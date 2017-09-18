#ifndef SPECEX_PSF_IO_GH2__H

void write_gauss_hermite_psf_fits_version_2(const specex::GaussHermitePSF& psf, fitsfile* fp);
void read_gauss_hermite_psf_fits_version_2(specex::PSF_p& psf, fitsfile* fp, int hdu, int first_bundle=-1, int last_bundle=-1);

#endif
