cmake_minimum_required(VERSION 2.8.12)
project(specex)
set(Python_ROOT_DIR "/global/common/software/desi/cori/desiconda/20200801-1.4.0-spec/conda")
find_package (Python 3.8 EXACT REQUIRED)
find_package (PythonInterp 3.8 EXACT REQUIRED)
find_package (PythonLibs 3.8 EXACT REQUIRED)
add_subdirectory(pybind11)
add_compile_options(-std=c++17 -std=c++14 -std=c++11 -pthread -B /global/common/cori_cle7/software/python/3.7-anaconda-2019.10/compiler_compat -Wl,--sysroot=/ -Wsign-compare -DNDEBUG -g -fwrapv -O3 -Wall -fPIC)
find_library(BOOSTHARP
	NAMES harp boost_regex-mt boost_program_options-mt
	      boost_serialization-mt
	PATHS /global/common/software/desi/cori/desiconda/20200801-1.4.0-spec/aux/lib
)
pybind11_add_module(specex
	    src/_pyspecex.cpp
	    src/specex_pyimage.cc
            src/specex_pyio.cc
	    src/specex_pyfitting.cc
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
target_include_directories(specex
	PUBLIC
	src
	/global/common/software/desi/cori/desiconda/20200801-1.4.0-spec/aux/lib)
target_link_libraries(specex PUBLIC "${BOOSTHARP}")

