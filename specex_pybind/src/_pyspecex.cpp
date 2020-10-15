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

namespace spx = specex;
namespace py  = pybind11;

PYBIND11_MAKE_OPAQUE(std::map<std::string,std::string>);

using ShapeContainer = py::detail::any_container<ssize_t>;

PYBIND11_MODULE(_internal, m) {
    m.doc() = R"(
    Internal wrapper around compiled specex code.

    The objects returned to python are wrapped in a shared_ptr, and wrapped
    functions are designed to take shared_ptrs as arguments.  If you are
    modifying this code, it is critical that shared_ptrs are used consistently
    everywhere.
    )";
    
    // stl bindings
   
    py::bind_map<std::map<std::string, std::string>>(m, "MapStringString");

    // classes

    py::class_ <spx::PyOptions, spx::PyOptions::pshr > (m, "PyOptions", R"(
        Class for storing and processing input options to specex.
        )")
        .def(py::init ())
        .def_readwrite("arc_image_filename", &spx::PyOptions::arc_image_filename)
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
        .def(py::init ());

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
