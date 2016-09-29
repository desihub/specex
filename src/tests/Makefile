
PROGS := $(patsubst %.cc,%,$(wildcard *.cc))

.PHONY : all install uninstall clean

all : $(PROGS)

install : all
	mkdir -p $(SPECEX_PREFIX)/bin
	cp $(PROGS) $(SPECEX_PREFIX)/bin

uninstall :
	@ for p in $(PROGS); do $(RM) "$(SPECEX_PREFIX)/bin/$${p}"; done

clean :
	- $(RM) $(PROGS) *.o *~


% : %.cc
	$(CXX) $(CXXFLAGS) -I../library -o $@ $< ../library/libspecex.a $(LINK)
