
OS := $(shell uname)

HEADERS = harp_plugin_specex.h

PLUGS = harp_plugin_specex.$(PLUG_EXT)
LIBS = harp_plugin_specex_lib.a

.PHONY : all install uninstall clean

all : $(PLUGS) $(LIBS)

install : all
	mkdir -p $(SPECEX_PREFIX)/lib
	cp $(PLUGS) $(SPECEX_PREFIX)/lib

uninstall :
	@ for p in $(PLUGS); do $(RM) "$(SPECEX_PREFIX)/lib/$${p}"; done

clean :
	- $(RM) $(PLUGS) *.o *~ *.$(PLUG_EXT) *.a


# We build two versions here, one is an internal library used by applications
# and the other is a loadable module.

harp_plugin_specex.$(PLUG_EXT) : harp_plugin_specex.o harp_plugin_specex_reg.o ../library/libspecex.a
	$(CXX) $(PLUG_FLAGS) $(PLUG_LINK) -o $@ $< harp_plugin_specex_reg.o ../library/libspecex.a $(LINK)

harp_plugin_specex_lib.a : harp_plugin_specex.o
	ar cr $@ harp_plugin_specex.o
	ranlib $@


%.o : %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -I. -I../library -fPIC -o $@ -c $<