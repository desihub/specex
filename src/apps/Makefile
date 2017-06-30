
PROGS := $(patsubst %.cc,%,$(wildcard *.cc))

.PHONY : all install uninstall clean

all : $(PROGS)

install : all
	mkdir -p $(SPECEX_PREFIX)/bin
	cp $(PROGS) $(SPECEX_PREFIX)/bin

uninstall :
	@ for p in $(PROGS); do $(RM) "$(SPECEX_PREFIX)/bin/$${p}"; done

clean :
	- $(RM) -f $(PROGS) *.o *~


specex_extract : specex_extract.cc ../plugin/harp_plugin_specex_lib.a ../library/libspecex.a
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I. -I../library -I../plugin -o $@ $< ../plugin/harp_plugin_specex_lib.a ../library/libspecex.a $(LINK)


% : %.cc ../library/libspecex.a
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I. -I../library -I../plugin -o $@ $< ../library/libspecex.a $(LINK)

