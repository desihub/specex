
PROGS = \
specex_cross_profile \
specex_extract \
specex_fiber_profile \
specex_flexure_fit_for_flat \
specex_flexure_fit_for_science \
specex_merge_psf \
specex_merge_spot \
specex_psf_fit \
specex_psf_residuals \
specex_spot_list


all : $(PROGS)

install : all
	@mkdir -p "$(SPECEX_PREFIX)/bin"; \
	cp $(PROGS) "$(SPECEX_PREFIX)/bin/"

uninstall :
	@for p in $(PROGS); do \
	rm -f "$(SPECEX_PREFIX)/bin/$${p}"; \
	done

clean :
	@rm -f $(PROGS) *.o *~


specex_extract : specex_extract.cc
	$(CXX) $(CXXFLAGS) $(PLUG_FLAGS) -I../library -I../plugin -o $@ $< ../plugin/harp_plugin_specex.$(PLUG_EXT) ../library/libspecex.$(PLUG_EXT) $(PLUG_LINK) $(LINK)

% : %.cc
	$(CXX) $(CXXFLAGS) -I../library -I../plugin -o $@ $< ../library/libspecex.$(PLUG_EXT) $(LINK)


