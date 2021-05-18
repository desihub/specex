#include <string>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

#include <specex_pyfitting.h>
#include <specex_pyoptions.h>
#include <specex_pyio.h>
#include <specex_pyimage.h>
#include <specex_psf.h>
#include <specex_gauss_hermite_psf.h>
#include <specex_psf_io.h>
#include <specex_image_data.h>
#include <specex_pypsf.h>
#include <specex_fits.h>

namespace spx = specex;
namespace py  = pybind11;

py::array_t<double>  image_data2nparray(spx::image_data image_data){
  auto v = image_data.data;
  size_t size = v.size();
  double *foo = new double[size];
  for (unsigned i = 0; i < v.size(); i++) foo[i]=v[i];	
  
  // Create a Python object that will free the allocated
  // memory when destroyed:
  py::capsule free_when_done(foo, [](void *f) {
      double *foo = reinterpret_cast<double *>(f);
      delete[] foo;
    });
  
  return py::array_t<double>(
			     {size}, // shape
			     {8}, // C-style//
			     foo, // the data pointer
			     free_when_done); // array refs parent
  
}

template <class T>
std::vector<T> vecassign(T var){
  std::vector<T> varvec;
  T vararr[] = { var } ;
  varvec.assign(vararr,vararr+1);
  return varvec;
}

using string_double_map  = std::map<std::string,std::string>;
using string_vstring_map = std::map<std::string,std::vector<std::string>>;
using string_vint_map    = std::map<std::string,std::vector<int>>;
using string_vdouble_map = std::map<std::string,std::vector<double>>;

using    vec_string   = std::vector<std::string>;
using    vec_int      = std::vector<int>;
using    vec_double   = std::vector<double>;
using    vec_longlong = std::vector<long long>;
using vecvec_int      = std::vector<std::vector<int>>;
using vecvec_double   = std::vector<std::vector<double>>;

PYBIND11_MAKE_OPAQUE(string_double_map);
PYBIND11_MAKE_OPAQUE(string_vstring_map);
PYBIND11_MAKE_OPAQUE(string_vint_map);
PYBIND11_MAKE_OPAQUE(string_vdouble_map);

PYBIND11_MAKE_OPAQUE(vec_string);
PYBIND11_MAKE_OPAQUE(vec_int);
PYBIND11_MAKE_OPAQUE(vec_longlong);
PYBIND11_MAKE_OPAQUE(vec_double);
PYBIND11_MAKE_OPAQUE(vecvec_int);
PYBIND11_MAKE_OPAQUE(vecvec_double);

//PYBIND11_MAKE_OPAQUE(varstring);
//PYBIND11_MAKE_OPAQUE(varlonglong);
//PYBIND11_MAKE_OPAQUE(vardouble);

using ShapeContainer = py::detail::any_container<ssize_t>;

PYBIND11_MODULE(_libspecex, m) {
    m.doc() = R"(
    Internal wrapper around compiled specex code.

    The objects returned to python are wrapped in a shared_ptr, and wrapped
    functions are designed to take shared_ptrs as arguments.  If you are
    modifying this code, it is critical that shared_ptrs are used consistently
    everywhere.
    )";

    // stl bindings
    
    py::bind_map<std::map<std::string,             std::string>> (m, "MapStringString");
    py::bind_map<std::map<std::string, std::vector<std::string>>>(m, "MapStringVString");
    py::bind_map<std::map<std::string, std::vector<int        >>>(m, "MapStringVInt");
    py::bind_map<std::map<std::string, std::vector<double     >>>(m, "MapStringVDouble");
    py::bind_vector<std::vector<std::string                    >>(m, "VectorString");
    py::bind_vector<std::vector<int>>                            (m, "VectorInt");
    py::bind_vector<std::vector<long long>>                      (m, "VectorLongLong");
    py::bind_vector<std::vector<double>>                         (m, "VectorDouble");
    py::bind_vector<std::vector<std::vector<int               >>>(m, "VectorVectorInt");
    py::bind_vector<std::vector<std::vector<double>>>            (m, "VectorVectorDouble");
    
    // data interface functions

    m.def("tablewrite_init", [](spx::PyPSF&          pyps) {

	pyps.FIBERMIN = pyps.psf->pydata.FIBERMIN;
	pyps.FIBERMAX = pyps.psf->pydata.FIBERMAX;
	
	pyps.trace_WAVEMIN  = pyps.psf->pydata.trace_WAVEMIN;
	pyps.trace_WAVEMAX  = pyps.psf->pydata.trace_WAVEMAX;
	pyps.table_WAVEMIN  = pyps.psf->pydata.table_WAVEMIN;
	pyps.table_WAVEMAX  = pyps.psf->pydata.table_WAVEMAX;

	pyps.nfibers = pyps.FIBERMAX - pyps.FIBERMIN + 1;
	
	pyps.mjd = (long long)pyps.psf->mjd;
	pyps.plate_id = (long long)pyps.psf->plate_id;
	pyps.camera_id       = (std::string)pyps.psf->camera_id;
	pyps.arc_exposure_id = (long long)pyps.psf->arc_exposure_id;
	pyps.NPIX_X          = (long long)pyps.psf->ccd_image_n_cols;
	pyps.NPIX_Y          = (long long)pyps.psf->ccd_image_n_rows;
	pyps.hSizeX          = (long long)pyps.psf->hSizeX;
	pyps.hSizeY          = (long long)pyps.psf->hSizeY;
	pyps.nparams_all     = (long long)(pyps.psf->LocalNAllPar()+2);
	pyps.ncoeff          = (long long)pyps.psf->ncoeff;
	pyps.GHDEGX          = (long long)pyps.psf->Degree();
	pyps.GHDEGY          = (long long)pyps.psf->Degree();
	pyps.psf_error       = (double   )pyps.psf->psf_error;
	pyps.readout_noise   = (double   )pyps.psf->readout_noise;
	pyps.gain            = (double   )pyps.psf->gain;	
    });

    m.def("get_trace", [](spx::PyPSF& pyps,
			  std::string axis) {
	    
	specex::image_data trace = pyps.get_trace(axis);
	pyps.trace_ncoeff  = trace.n_cols();
	return image_data2nparray(trace);
	    
    });
    
    m.def("get_table", [](spx::PyPSF&  pyps,
			  std::vector<std::string> &table_col0,
			  std::vector<double>      &table_col1,
			  std::vector<int>         &table_col2,
			  std::vector<int>         &table_col3,
			  std::vector<int>         &bundle_id,
			  std::vector<int>         &bundle_ndata,
			  std::vector<int>         &bundle_nparams,
			  std::vector<double>      &bundle_chi2pdf
			  ) {

	pyps.SetParamsOfBundle();

	bundle_id      = pyps.bundle_id;
	bundle_ndata   = pyps.bundle_ndata;
	bundle_nparams = pyps.bundle_nparams;
	bundle_chi2pdf = pyps.bundle_chi2pdf;
	    
	spx::FitsTable table = pyps.psf->pydata.table;
	pyps.table_nrows = table.data.size();
	int nrows = pyps.table_nrows;
	
	std::map<std::string,spx::FitsColumnDescription>::const_iterator c1 =
	  table.columns.begin();

	for(int r=0;r<nrows;r++){

	  // row0, PARAM
	  const char *input_val = table.data[r][0].string_val.c_str();	  
	  table_col0.push_back(std::string(input_val));

	  // row2, COEFF	  
	  const spx::FitsColumnDescription& column = c1->second;
	  int nd = column.SizeOfVectorOfDouble();
	  for(int d=0;d<nd;d++) table_col1.push_back(table.data[r][1].double_vals[d]);

	  // row3, LEGDEGX
	  table_col2.push_back(table.data[r][2].int_vals(0));

	  // row4, LEGDEGW
	  table_col3.push_back(table.data[r][3].int_vals(0));
	  
	}
	
    });
    
    // classes

    py::class_ <spx::PSF, spx::PSF::pshr> (m, "PSF", R"(PSF base class)");
    
    py::class_ <spx::GaussHermitePSF, spx::GaussHermitePSF::pshr> (m,
					         "GaussHermitePSF", R"(
        Class for storing and processing PSF in specex.
        )")
        .def(py::init ())
        .def("Degree", &spx::GaussHermitePSF::Degree);
    
    py::class_ <spx::PyOptions, spx::PyOptions::pshr > (m, "PyOptions", R"(
        Class for storing and processing input options to specex.
        )")
        .def(py::init ())
        .def_readwrite("arc_image_filename",   &spx::PyOptions::arc_image_filename)
        .def_readwrite("input_psf_filename",   &spx::PyOptions::input_psf_filename)
        .def_readwrite("output_fits_filename", &spx::PyOptions::output_fits_filename)
        .def_readwrite("trace_deg_x",          &spx::PyOptions::trace_deg_x)
        .def_readwrite("trace_deg_wave",       &spx::PyOptions::trace_deg_wave)

        .def("parse", [](spx::PyOptions &self, std::vector<std::string>& args){
	    std::vector<char *> cstrs;
	    cstrs.reserve(args.size());
	    for (auto &s : args) cstrs.push_back(const_cast<char *>(s.c_str()));
	    return self.parse(cstrs.size(), cstrs.data());
	}
	);
    
    py::class_ <spx::PyImage, spx::PyImage::pshr > (m, "PyImage", R"(
        Class for specex image interchangeable with python implementations
        )")
        .def(py::init ())
        .def(py::init<py::array,py::array,py::array,py::array,
	     std::map<std::string,std::string>>())
        .def("get_header", &spx::PyImage::get_header)
        .def("get_data", [](spx::PyImage &self, std::string tag){
	  auto v = self.get_data(tag);
	  return py::array(v.size(),v.data());
	}
	);

    py::class_ <spx::PyPSF, spx::PyPSF::pshr > (m, "PyPSF", R"(
        Class for specex PSF interchangeable with python implementations
        )")
        .def(py::init ())

        .def("set_trace",          &spx::PyPSF::set_trace)
        .def("synchronize_traces", &spx::PyPSF::synchronize_traces)

        .def_readwrite("psf", &spx::PyPSF::psf)

        .def_readwrite("nfibers",          &spx::PyPSF::nfibers)    

        .def_readwrite("table_nrows",      &spx::PyPSF::table_nrows)    
        .def_readwrite("trace_ncoeff",     &spx::PyPSF::trace_ncoeff)
    
        .def_readwrite("FIBERMIN",         &spx::PyPSF::FIBERMIN)
        .def_readwrite("FIBERMAX",         &spx::PyPSF::FIBERMAX)
        .def_readwrite("trace_WAVEMIN",    &spx::PyPSF::trace_WAVEMIN)
        .def_readwrite("trace_WAVEMAX",    &spx::PyPSF::trace_WAVEMAX)
        .def_readwrite("table_WAVEMIN",    &spx::PyPSF::table_WAVEMIN)
        .def_readwrite("table_WAVEMAX",    &spx::PyPSF::table_WAVEMAX)
      
        .def_readwrite("mjd",              &spx::PyPSF::mjd)
        .def_readwrite("plate_id",         &spx::PyPSF::plate_id)
        .def_readwrite("arc_exposure_id",  &spx::PyPSF::arc_exposure_id)
        .def_readwrite("NPIX_X",           &spx::PyPSF::NPIX_X)
        .def_readwrite("NPIX_Y",           &spx::PyPSF::NPIX_Y)
        .def_readwrite("hSizeX",           &spx::PyPSF::hSizeX)
        .def_readwrite("hSizeY",           &spx::PyPSF::hSizeY)
        .def_readwrite("nparams_all",      &spx::PyPSF::nparams_all)
        .def_readwrite("ncoeff",           &spx::PyPSF::ncoeff)
        .def_readwrite("GHDEGX",           &spx::PyPSF::GHDEGX)
        .def_readwrite("GHDEGY",           &spx::PyPSF::GHDEGY)
        .def_readwrite("camera_id",        &spx::PyPSF::camera_id)
        .def_readwrite("psf_error",        &spx::PyPSF::psf_error)
        .def_readwrite("readout_noise",    &spx::PyPSF::readout_noise)
        .def_readwrite("gain",             &spx::PyPSF::gain);

    py::class_ <spx::PyIO, spx::PyIO::pshr > (m, "PyIO", R"(
        Class for specex IO interchangeable with python implementations
        )")
        .def(py::init ())
        .def("set_inputpsf", [](spx::PyIO &self, spx::PyOptions opts){
	    return self.set_inputpsf(opts);
	}
	)
        .def("read_preproc", [](spx::PyIO &self, spx::PyOptions opts,
				   spx::PyImage& pyimg){
	    return self.read_preproc(opts,pyimg);
	}
	)
        .def("read_psf", [](spx::PyIO &self, spx::PyOptions opts,
				 spx::PyPSF& pypsf){
	    return self.read_psf(opts,pypsf);
	}
	)
        .def("read_traces", [](spx::PyIO &self, spx::PyOptions opts,
				 spx::PyPSF& pypsf){
	    return self.read_traces(opts,pypsf);
	}
	)
        .def("init_traces", [](spx::PyIO &self, spx::PyOptions opts,
				 spx::PyPSF& pypsf){
	    return self.init_traces(opts,pypsf);
	}
	)
        .def("load_psf",[](spx::PyIO &self, spx::PyOptions opts,
				 spx::PyPSF& pypsf){
	    return self.load_psf(opts,pypsf);
	}
	)
        .def("write_spots",   [](spx::PyIO &self, spx::PyOptions opts,
				 spx::PyPSF& pypsf){
	    return self.write_spots(opts,pypsf);
	}
	);

    py::class_ <spx::PyPrior, spx::PyPrior::pshr > (m, "PyPrior", R"(
        Class for specex priors interchangeable with python implementations
        )")
        .def(py::init ())
        .def("set_priors", [](spx::PyPrior &self, spx::PyOptions opts){
	  return self.set_priors(opts);
	}
	);

    py::class_ <spx::PyFitting, spx::PyFitting::pshr > (m, "PyFitting", R"(
        Class for fitting PSFs in specex.
        )")
        .def(py::init ())
        .def("fit_psf", [](spx::PyFitting &self,
			   spx::PyOptions opts,
			   spx::PyIO      pyio,
			   spx::PyPrior   pypr,
			   spx::PyImage   pymg,
			   spx::PyPSF     pyps){
	       return self.fit_psf(opts,pyio,pypr,pymg,pyps);
	}
	);
    

}
