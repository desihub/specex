
HEADERS = harp_plugin_specex.h

PLUGS = harp_plugin_specex.$(PLUG_EXT)
ifneq "x$(MPICXX)" "x"
  PLUGS += harp_mpi_plugin_specex.$(PLUG_EXT)
endif 


all : $(PLUGS)

install : all
	@mkdir -p "$(SPECEX_PREFIX)/lib"; \
	cp $(PLUGS) "$(SPECEX_PREFIX)/lib/"

uninstall :
	@for p in $(PLUGS); do \
	rm -f "$(SPECEX_PREFIX)/lib/$${p}"; \
	done

clean :
	@rm -f $(PLUGS) *.o *~


harp_plugin_specex.$(PLUG_EXT) : harp_plugin_specex.o ../library/libspecex.$(PLUG_EXT)
	$(CXX) $(PLUG_FLAGS) $(PLUG_LINK) -o $@ $< ../library/libspecex.$(PLUG_EXT) $(LINK)


harp_mpi_plugin_specex.$(PLUG_EXT) : harp_mpi_plugin_specex.o ../library/libspecex_mpi.$(PLUG_EXT)
	$(MPICXX) $(PLUG_FLAGS) $(PLUG_LINK) -o $@ $< ../library/libspecex_mpi.$(PLUG_EXT) $(MPILINK)


harp_plugin_specex.o : harp_plugin_specex.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) $(PLUG_FLAGS) -I../library -o $@ -c $<


harp_mpi_plugin_specex.o : harp_plugin_specex.cc $(HEADERS)
	$(MPICXX) $(CXXFLAGS) $(PLUG_FLAGS) -DUSE_MPI -I../library -o $@ -c $<



