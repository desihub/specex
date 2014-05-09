
HEADERS = \
specex_boss_io.h \
specex_brent.h \
specex_fits.h \
specex_gauss_hermite_psf.h \
specex_gauss_hermite_two_psf.h \
specex_hat_hermite_psf.h \
specex_hat_moffat_psf.h \
specex_hermite.h \
specex_image_data.h \
specex_lamp_lines_utils.h \
specex_legendre.h \
specex_linalg.h \
specex_mask.h \
specex_message.h \
specex_model_image.h \
specex_psf.h \
specex_psf_fitter.h \
specex_psf_io.h \
specex_serialization.h \
specex_spectrograph.h \
specex_spot.h \
specex_spot_array.h \
specex_stamp.h \
specex_trace.h \
specex_vector_utils.h

OBJS = \
specex_boss_io.o \
specex_brent.o \
specex_fits.o \
specex_gauss_hermite_psf.o \
specex_gauss_hermite_two_psf.o \
specex_hat_hermite_psf.o \
specex_hat_moffat_psf.o \
specex_hermite.o \
specex_image_data.o \
specex_lamp_lines_utils.o \
specex_legendre.o \
specex_linalg.o \
specex_mask.o \
specex_message.o \
specex_model_image.o \
specex_psf.o \
specex_psf_fitter.o \
specex_psf_io.o \
specex_spectrograph.o \
specex_spot.o \
specex_spot_array.o \
specex_stamp.o \
specex_trace.o \
specex_vector_utils.o

LIBS = libspecex.$(PLUG_EXT)

all : $(LIBS)

libspecex.$(PLUG_EXT) : $(OBJS)
	$(CXX) $(PLUG_LINK) -o $@ $(OBJS)

install : all
	@mkdir -p "$(SPECEX_PREFIX)/lib"; \
	cp $(LIBS) "$(SPECEX_PREFIX)/lib/"

uninstall :
	@rm -f "$(SPECEX_PREFIX)/lib/libspecex*"

clean :
	@rm -f $(LIBS) *.o *~


%.o : %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<


