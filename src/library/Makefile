
HEADERS := $(wildcard *.h)

OBJS := $(patsubst %.cc,%.o,$(wildcard *.cc))

LIBS = libspecex.a libspecex.$(PLUG_EXT)

.PHONY : all install uninstall clean

all : $(LIBS)

libspecex.a : $(OBJS)
	ar cr $@ $(OBJS)
	ranlib $@

libspecex.$(PLUG_EXT) : $(OBJS)
	$(CXX) -Wl,--no-undefined -fPIC -shared -o $@ $(OBJS) $(LINK)

install : all
	mkdir -p $(SPECEX_PREFIX)/lib
	cp $(LIBS) $(SPECEX_PREFIX)/lib

uninstall :
	- $(RM) "$(SPECEX_PREFIX)/lib/libspecex*"

clean :
	- $(RM) $(LIBS) *.o *~


%.o : %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -fPIC -o $@ -c $<
