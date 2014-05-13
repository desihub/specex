
HEADERS = harp_plugin_specex.h

PLUGS = harp_plugin_specex.$(PLUG_EXT)


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


harp_plugin_specex.$(PLUG_EXT) : harp_plugin_specex.o ../library/libspecex.a
	$(CXX) $(PLUG_FLAGS) $(PLUG_LINK) -o $@ $< ../library/libspecex.a $(LINK)

harp_plugin_specex.o : harp_plugin_specex.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) $(PLUG_FLAGS) -I../library -o $@ -c $<

