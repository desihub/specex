
PROGS = \
specex_psf_image \
specex_test_boost \
specex_test_derivatives \
specex_test_derivatives2 \
specex_test_psf_specter_io \
specex_test_psf_xml_io


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


% : %.cc
	$(CXX) $(CXXFLAGS) -I../library -o $@ $< ../library/libspecex.a $(LINK)

