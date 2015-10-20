
PROGS := $(patsubst %.cc,%,$(wildcard *.cc))


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

