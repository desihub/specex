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

namespace spx = specex;
namespace py  = pybind11;

py::array_t<double>  image_data2nparray(spx::image_data image_data){
  auto v = image_data.data;
  cout << "size in get_pyps_image: " << v.size() << endl;
  size_t size = v.size();
  double *foo = new double[size];
  for (unsigned i = 0; i < v.size(); i++) foo[i]=v[i];	
  
  // Create a Python object that will free the allocated
  // memory when destroyed:
  py::capsule free_when_done(foo, [](void *f) {
      double *foo = reinterpret_cast<double *>(f);
      std::cerr << "Element [0] = " << foo[0] << "\n";
      std::cerr << "freeing memory @ " << f << "\n";
      delete[] foo;
    });
  
  return py::array_t<double>(
			     {size}, // shape
			     {8}, // C-style//
			     foo, // the data pointer
			     free_when_done); // array refs parent
  
}

using string_double_map  = std::map<std::string,std::string>;
using string_vstring_map = std::map<std::string,std::vector<std::string>>;
using string_vint_map    = std::map<std::string,std::vector<int>>;
using string_vdouble_map = std::map<std::string,std::vector<double>>;

using vec_double = std::vector<double>;
using vec_int    = std::vector<int>;

PYBIND11_MAKE_OPAQUE(string_double_map);
PYBIND11_MAKE_OPAQUE(string_vstring_map);
PYBIND11_MAKE_OPAQUE(string_vint_map);
PYBIND11_MAKE_OPAQUE(string_vdouble_map);

PYBIND11_MAKE_OPAQUE(vec_double);
PYBIND11_MAKE_OPAQUE(vec_int);

using ShapeContainer = py::detail::any_container<ssize_t>;

PYBIND11_MODULE(specex, m) {
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
    py::bind_vector<std::vector<double>>                         (m, "VectorDouble");
    py::bind_vector<std::vector<int>>                            (m, "VectorInt");

    // data interface functions

    m.def("get_tracekeys", [](spx::PyPSF pyps,
			      std::vector<int>    &fiberkeys,
			      std::vector<double> &wavekeys) {

	int   fiberarr[] = { pyps.psf->pydata.FIBERMIN, pyps.psf->pydata.FIBERMAX };
	double wavearr[] = { pyps.psf->pydata.WAVEMIN,  pyps.psf->pydata.WAVEMAX  };
	
	fiberkeys.assign(fiberarr,fiberarr+2);
	wavekeys.assign(  wavearr, wavearr+2);
	
    });
    
    m.def("get_trace", [](spx::PyPSF  pyps,
			  std::string axis) {

	 specex::image_data trace = pyps.get_trace(axis);
	 return image_data2nparray(trace);
	    
    });
    
    m.def("get_tablekeys", [](spx::PyPSF pyps,
			      std::map<std::string, std::string> &tablekeys_comment,
			      std::map<std::string, std::string> &tablekeys_string,
			      std::map<std::string, int>         &tablekeys_int,
			      std::map<std::string, double>      &tablekeys_double) {

	    tablekeys_comment = pyps.psf->pydata.tablekeys_comment;
	    tablekeys_string  = pyps.psf->pydata.tablekeys_string;
	    tablekeys_int     = pyps.psf->pydata.tablekeys_int;
	    tablekeys_double  = pyps.psf->pydata.tablekeys_double;
	    
    });

    m.def("get_table", [](spx::PyPSF  pyps,
	                  std::vector<std::string>         &table_string,
			  std::vector<std::vector<int>>    &table_int,
			  std::vector<std::vector<double>> &table_double) {

	 table_string = pyps.get_tablestring();
	 table_int    = pyps.get_tableint();
	 table_double = pyps.get_tabledouble();
	    
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
        .def_readwrite("arc_image_filename", &spx::PyOptions::arc_image_filename)
        .def_readwrite("input_psf_filename", &spx::PyOptions::input_psf_filename)
        .def_readwrite("output_fits_filename", &spx::PyOptions::output_fits_filename)
        .def("parse", [](spx::PyOptions &self, std::vector<std::string> args){
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
        .def_readwrite("psf", &spx::PyPSF::psf);

    py::class_ <spx::PyIO, spx::PyIO::pshr > (m, "PyIO", R"(
        Class for specex IO interchangeable with python implementations
        )")
        .def(py::init ())
        .def("check_input_psf", [](spx::PyIO &self, spx::PyOptions opts){
	    return self.check_input_psf(opts);
	}
	)
        .def("read_img_data", [](spx::PyIO &self, spx::PyOptions opts,
				   spx::PyImage& pyimg){
	    return self.read_img_data(opts,pyimg);
	}
	)
        .def("read_psf_data", [](spx::PyIO &self, spx::PyOptions opts,
				 spx::PyPSF& pypsf){
	    return self.read_psf_data(opts,pypsf);
	}
	)
        .def("write_psf_data",[](spx::PyIO &self, spx::PyOptions opts,
				 spx::PyPSF& pypsf){
	    return self.write_psf_data(opts,pypsf);
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
        .def("deal_with_priors", [](spx::PyPrior &self, spx::PyOptions opts){
	  return self.deal_with_priors(opts);
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
