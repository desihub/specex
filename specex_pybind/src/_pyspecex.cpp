#include <string>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <pybind11/stl_bind.h>

#include <mymath.h>
#include <myargs.h>
#include <ext.h>
#include <specex_desi.h>

#include <specex_pyoptions.h>
#include <specex_pyio.h>
#include <specex_pyimage.h>

namespace spx = specex;
namespace py  = pybind11;

using ShapeContainer = py::detail::any_container<ssize_t>;

PYBIND11_MODULE(_internal, m) {
    m.doc() = R"(
    Internal wrapper around compiled specex code.

    The objects returned to python are wrapped in a shared_ptr, and wrapped
    functions are designed to take shared_ptrs as arguments.  If you are
    modifying this code, it is critical that shared_ptrs are used consistently
    everywhere.
    )";

    m.def("specex_desi_psf_fit_main", [](spx::PyOptions opts,
					 spx::PyIO      pyio,
					 spx::PyPrior   pypr,
					 spx::PyImage   pymg,
					 spx::PyPSF     pyps){
	    return specex_desi_psf_fit_main(opts,pyio,pypr,pymg,pyps);
    }
    );
    
    // wrap the valid static definitions

    m.attr("MYMATH_INT") = py::int_(MYMATH_INT);
    m.attr("MYARGS_INT") = py::int_(MYARGS_INT);

    // classes
    
    py::class_ <spx::PyOptions, spx::PyOptions::pshr > (m, "PyOptions", R"(
        Class for storing and processing input options to specex.
        )")
        .def(py::init ())
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
        .def(py::init ());

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
				 spx::PyPSF& pypsf, spx::PyIO pyio){
	    return self.read_psf_data(opts,pypsf,pyio);
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

    py::class_ <spx::MyMath, spx::MyMath::pshr > (m, "MyMath", R"(
        Simple mymath class.

        This class is just a mymath to do math type tests.
        )")
        .def(py::init <int> ())
        .def("add", &spx::MyMath::add);
    
    py::class_ <spx::MyArgs, spx::MyArgs::pshr > (m, "MyArgs", R"(
        Simple myargs class.

        This class is just a myargs to do args type tests.
        )")
        .def(py::init <int> ())
        .def("print_args", &spx::MyArgs::print_args)
        .def("get_args", [](spx::MyArgs &self, std::vector<std::string> args){
	    std::vector<char *> cstrs;
	    cstrs.reserve(args.size());
	    for (auto &s : args) cstrs.push_back(const_cast<char *>(s.c_str()));
	    return self.get_args(cstrs.size(), cstrs.data());
	}
	);

}
