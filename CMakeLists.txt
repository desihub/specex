cmake_minimum_required(VERSION 2.8.12)

project(specex)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(Python_ROOT_DIR $ENV{DESICONDA})

find_package (Python 3.8 REQUIRED)
find_package (PythonInterp 3.8 REQUIRED)
find_package (PythonLibs 3.8 REQUIRED)

add_subdirectory(src/pybind11)
pybind11_add_module(_libspecex NO_EXTRAS
	    src/_libspecex.cpp
	    src/specex_psfpy.cc
	    src/specex_pyimage.cc
            src/specex_pyio.cc
	    src/specex_pyfitting.cc
	    src/specex_pypsf.cc
	    src/specex_pyprior.cc
	    src/specex_pyoptions.cc
	    src/specex_gauss_hermite_psf.cc
	    src/specex_psf.cc
	    src/specex_psf_fitter.cc
	    src/specex_message.cc
	    src/specex_legendre.cc
	    src/specex_linalg.cc
	    src/specex_psf_proc.cc
	    src/specex_model_image.cc
	    src/specex_image_data.cc
	    src/specex_hermite.cc
	    src/specex_trace.cc
	    src/specex_spot_array.cc
	    src/specex_spot.cc
	    src/specex_mask.cc
	    src/specex_lamp_lines_utils.cc
	    src/specex_brent.cc
	    src/specex_tokens.cc
	    src/specex_vector_utils.cc
	    src/specex_blas.c
	    src/specex_lapack.c
	    src/specex_image_data_base.cc
    	    src/specex_fits.cc
	    src/specex_unbst.cc
	    src/specex_unbls.cc
	    )

# look for BLAS/LAPACK with MKL first

set(BLA_VENDOR Intel10_64lp)
find_package (BLAS)
find_package (LAPACK)
if(BLAS_FOUND)

	# this also means MKL with BLAS libraries were found
	# now see if MKL include file can be found as well
	if(EXISTS "$ENV{DESICONDA}/include/mkl.h")
	    set(MKL_FOUND TRUE)
	    set(MKL_INCLUDE "$ENV{DESICONDA}/include")
	    target_compile_definitions(_libspecex PUBLIC USE_MKL)
	elseif(EXISTS "/opt/intel/mkl/include/mkl.h")
	    set(MKL_FOUND TRUE)
	    set(MKL_INCLUDE "/opt/intel/mkl/include")
	    target_compile_definitions(_libspecex PUBLIC USE_MKL)
	endif()

endif()

if(NOT MKL_FOUND)

	# if MKL BLAS and mkl.h not both found, look for OpenBLAS/LAPACK
	set(BLAS_FOUND FALSE)
	set(BLA_VENDOR OpenBLAS)
	find_package (BLAS)
	find_package (LAPACK)

endif()

if(NOT BLAS_FOUND)
       message (FATAL_ERROR "BLAS not found, CMake will exit.")
endif()

message ("MKL found: ${MKL_FOUND}")
message ("MKL include: ${MKL_INCLUDE}")
message ("BLAS libraries: ${BLAS_LIBRARIES}")
message ("LAPACK libraries: ${LAPACK_LIBRARIES}")

target_include_directories(_libspecex
	PUBLIC
	src
	${MKL_INCLUDE}
	)

target_link_directories(_libspecex
	PUBLIC)

target_link_libraries(_libspecex PUBLIC ${BLAS_LIBRARIES})

set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS})

set(INTERPROCEDURAL_OPTIMIZATION FALSE)
set(INTERPROCEDURAL_OPTIMIZATION_<CONFIG> FALSE)
set(CMAKE_C_COMPILE_OPTIONS_IPO "")
target_compile_options(_libspecex PRIVATE -w)
