
HEADERS := $(wildcard *.h)

OBJS := $(patsubst %.cc,%.o,$(wildcard *.cc))

LIBS = libspecex.a libspecex.$(PLUG_EXT)

OS := $(shell uname)

ifeq ($(OS), Darwin)
  NOUNDEF = -Wl,-undefined,error
else
  NOUNDEF = -Wl,--no-undefined
endif


.PHONY : all install uninstall clean

all : $(LIBS)

libspecex.a : $(OBJS)
	ar cr $@ $(OBJS)
	ranlib $@

libspecex.$(PLUG_EXT) : $(OBJS)
	$(CXX) $(NOUNDEF) -fPIC -shared -o $@ $(OBJS) $(LINK)

install : all
	mkdir -p $(SPECEX_PREFIX)/lib
	cp $(LIBS) $(SPECEX_PREFIX)/lib

uninstall :
	- $(RM) "$(SPECEX_PREFIX)/lib/libspecex*"

clean :
	- $(RM) $(LIBS) *.o *~


%.o : %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -I. -fPIC -o $@ -c $<
