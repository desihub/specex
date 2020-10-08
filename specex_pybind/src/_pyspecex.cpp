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

#include <specex_options.h>
#include <specex_pyio.h>

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

    //m.def("specex_desi_psf_fit_main",&specex_desi_psf_fit_main);
    m.def("specex_desi_psf_fit_main", [](spx::Options opts){
	return specex_desi_psf_fit_main(opts);
    }
    );
    /*
    m.def("specex_desi_psf_fit_main", [](std::vector<std::string> args){
	    std::vector<char *> cstrs;
	    cstrs.reserve(args.size());
	    for (auto &s : args) cstrs.push_back(const_cast<char *>(s.c_str()));
	    return specex_desi_psf_fit_main(cstrs.size(), cstrs.data());
    }
    );
    */
    
    // wrap the valid static definitions

    m.attr("MYMATH_INT") = py::int_(MYMATH_INT);
    m.attr("MYARGS_INT") = py::int_(MYARGS_INT);

    // classes
    
    py::class_ <spx::Options, spx::Options::pshr > (m, "Options", R"(
        Class for storing and processing input options to specex.
        )")
        .def(py::init ())
        .def("parse", [](spx::Options &self, std::vector<std::string> args){
	    std::vector<char *> cstrs;
	    cstrs.reserve(args.size());
	    for (auto &s : args) cstrs.push_back(const_cast<char *>(s.c_str()));
	    return self.parse(cstrs.size(), cstrs.data());
	}
	);
    
    py::class_ <spx::PyIO, spx::PyIO::pshr > (m, "PyIO", R"(
        Class for specex IO interchangeable with python implementations
        )")
        .def(py::init ())
        .def("check_input_psf", [](spx::PyIO &self, spx::Options opts){
	    return self.check_input_psf(opts);
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
