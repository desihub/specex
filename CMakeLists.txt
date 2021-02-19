cmake_minimum_required(VERSION 2.8.12)
project(specex)
set(Python_ROOT_DIR $ENV{DESICONDA})
find_package (Python 3.8 EXACT REQUIRED)
find_package (PythonInterp 3.8 EXACT REQUIRED)
find_package (PythonLibs 3.8 EXACT REQUIRED)
add_subdirectory(src/pybind11)
add_compile_options(-pthread -Wl,--sysroot=/ -Wsign-compare -DNDEBUG -g -fwrapv -O3 -Wall)
pybind11_add_module(_libspecex NO_EXTRAS 
	    src/_libspecex.cpp
	    src/specex_psfpy.cc
	    src/specex_pyimage.cc
            src/specex_pyio.cc
	    src/specex_pyfitting.cc
	    src/specex_pypsf.cc
	    src/specex_pyprior.cc
	    src/specex_tokens.cc
	    src/specex_pyoptions.cc
	    src/mymath.cpp
	    src/myargs.cpp
	    src/ext.cpp
	    src/specex_gauss_hermite_psf.cc
	    src/specex_psf.cc
	    src/specex_psf_fitter.cc
	    src/specex_message.cc
	    src/specex_legendre.cc
	    src/specex_linalg.cc
	    src/specex_psf_io.cc
	    src/specex_model_image.cc
	    src/specex_fits.cc
	    src/specex_image_data.cc
	    src/specex_hermite.cc
	    src/specex_trace.cc
	    src/specex_spot_array.cc
	    src/specex_desi_io.cc
	    src/specex_psf_io_gauss_hermite_psf_fits_1.cc
	    src/specex_psf_io_gauss_hermite_psf_fits_2.cc
	    src/specex_psf_io_gauss_hermite_psf_fits_3.cc
	    src/specex_spot.cc
	    src/specex_vector_utils.cc
	    src/specex_mask.cc
	    src/specex_lamp_lines_utils.cc
	    src/specex_brent.cc)
target_include_directories(_libspecex
	PUBLIC
	src
	src/deboosted
	)

target_link_directories(_libspecex
	PUBLIC)
target_link_libraries(_libspecex PUBLIC "-lharp -lboost_regex-mt -lboost_program_options-mt -lboost_serialization-mt")
set(INTERPROCEDURAL_OPTIMIZATION FALSE)
set(INTERPROCEDURAL_OPTIMIZATION_<CONFIG> FALSE)
set(CMAKE_C_COMPILE_OPTIONS_IPO "")
target_compile_options(_libspecex PRIVATE -w)
