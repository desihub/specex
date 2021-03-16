cmake_minimum_required(VERSION 2.8.12)

project(specex)

set(Python_ROOT_DIR $ENV{DESICONDA})

find_package (Python 3.8 EXACT REQUIRED)
find_package (PythonInterp 3.8 EXACT REQUIRED)
find_package (PythonLibs 3.8 EXACT REQUIRED)

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
	    src/specex_psf_io.cc
	    src/specex_model_image.cc
	    src/specex_image_data.cc
	    src/specex_hermite.cc
	    src/specex_trace.cc
	    src/specex_spot_array.cc
	    src/specex_psf_io_gauss_hermite_psf_fits_1.cc
	    src/specex_psf_io_gauss_hermite_psf_fits_2.cc
	    src/specex_psf_io_gauss_hermite_psf_fits_3.cc
	    src/specex_spot.cc
	    src/specex_mask.cc
	    src/specex_lamp_lines_utils.cc
	    src/specex_fits.cc
	    src/deboosted/specex_brent.cc
	    src/deboosted/specex_desi_io.cc
	    src/deboosted/specex_tokens.cc
	    src/deboosted/specex_vector_utils.cc
	    src/specex_blas.c)
# look for BLAS/LAPACK with MKL first

set(BLA_VENDOR Intel10_64lp)
#set(BLA_VENDOR Intel10_64_dyn)
find_package (BLAS)
find_package (LAPACK)

if(BLAS_FOUND)

	# this also means MKL was found
	# currently only "standard" MKL location supported
	set(MKL_FOUND TRUE)
	set(MKL_INCLUDE /opt/intel/mkl/include)
	target_compile_definitions(_libspecex PUBLIC USE_MKL)
	
else()

	# otherwise look for generic BLAS/LAPACK
	
	set(BLA_VENDOR "Generic")
	set(BLA_VENDOR OpenBLAS)
	find_package (BLAS)
	find_package (LAPACK)	

endif()

#	message ("BLAS found: ${BLAS_FOUND}")
	message ("BLAS libraries: ${BLAS_LIBRARIES}")
#	message ("BLAS linker flags: ${BLAS_LINKER_FLAGS}")

#	message ("LAPACK found: ${LAPACK_FOUND}")
	message ("LAPACK libraries: ${LAPACK_LIBRARIES}")
#	message ("LAPACK linker flags: ${LAPACK_LINKER_FLAGS}")

if(NOT BLAS_FOUND)
       message (FATAL_ERROR "BLAS not found, CMake will exit.")
endif()

#set(MKLROOT /opt/intel/mkl)

#set(MKLLINK "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl")

#set(MKLCOMPILE -DMKL_ILP64  -I${MKLROOT}/include)

#add_compile_options(-pthread -Wl,--sysroot=/ -Wsign-compare -DNDEBUG -g -fwrapv -O3 -Wall ${MKLCOMPILE})

target_include_directories(_libspecex
	PUBLIC
	src
	src/deboosted
	${MKL_INCLUDE}
	)

target_link_directories(_libspecex
	PUBLIC)
target_link_libraries(_libspecex PUBLIC "-lharp -lboost_regex-mt -lboost_program_options-mt -lboost_serialization-mt" ${BLAS_LIBRARIES})
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS})
#set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MKLLINK})
set(INTERPROCEDURAL_OPTIMIZATION FALSE)
set(INTERPROCEDURAL_OPTIMIZATION_<CONFIG> FALSE)
set(CMAKE_C_COMPILE_OPTIONS_IPO "")
target_compile_options(_libspecex PRIVATE -w)
