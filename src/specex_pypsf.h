#ifndef SPECEX_PYPSF__H
#define SPECEX_PYPSF__H

#include <vector>
#include <string>

#include <specex_unbls.h>

#include <specex_psf.h>
#include <specex_spot.h>

#include <specex_gauss_hermite_psf.h>

#include <specex_pyoptions.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace specex {

  class PyPSF : public std::enable_shared_from_this <PyPSF> {

  public :

    typedef std::shared_ptr <PyPSF> pshr;

    specex::PSF_p psf;
    vector<Spot_p> fitted_spots;

    PyPSF(){}

    int trace_ncoeff, table_nrows, nfibers;
    int FIBERMIN, FIBERMAX, LEGDEG, TRDEGW;
    double trace_WAVEMIN, trace_WAVEMAX;
    double table_WAVEMIN, table_WAVEMAX;
    long long mjd, plate_id, arc_exposure_id, NPIX_X, NPIX_Y,
      hSizeX, hSizeY, nparams_all, ncoeff, GHDEGX, GHDEGY;
    std::string camera_id;
    double psf_error, readout_noise, gain;
    std::vector<int>    bundle_id, bundle_ndata, bundle_nparams;
    std::vector<double> bundle_chi2pdf;

    image_data get_trace(std::string);
    void set_trace(py::array, int, int);
    void init_traces(specex::PyOptions);
    void synchronize_traces();
    void set_psf(
		 std::vector<std::string>&,
		 std::vector<double>&,
		 std::vector<int>&,
		 std::vector<int>&
		 );

    void SetParamsOfBundle();

  };

}

#endif
