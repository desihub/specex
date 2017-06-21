
OS := $(shell uname)

HEADERS = harp_plugin_specex.h

PLUGS = harp_plugin_specex.$(PLUG_EXT)

.PHONY : all install uninstall clean

all : $(PLUGS)

install : all
	mkdir -p $(SPECEX_PREFIX)/lib
	cp $(PLUGS) $(SPECEX_PREFIX)/lib

uninstall :
	@ for p in $(PLUGS); do $(RM) "$(SPECEX_PREFIX)/lib/$${p}"; done

clean :
	- $(RM) $(PLUGS) *.o *~

ifeq ($(OS), Darwin)

harp_plugin_specex.$(PLUG_EXT) : harp_plugin_specex.o ../library/libspecex.a
	touch harp_plugin_specex.$(PLUG_EXT)

else

harp_plugin_specex.$(PLUG_EXT) : harp_plugin_specex.o ../library/libspecex.a
	$(CXX) $(PLUG_FLAGS) $(PLUG_LINK) -o $@ $< ../library/libspecex.a $(LINK)

endif

harp_plugin_specex.o : harp_plugin_specex.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) $(PLUG_FLAGS) -I../library -o $@ -c $<

