
HEADERS := $(wildcard *.h)

OBJS := $(patsubst %.cc,%.o,$(wildcard *.cc))

LIBS = libspecex.a libspecex.$(PLUG_EXT)

OS := $(shell uname)
ifeq ($(OS), Darwin)
  NOUNDEF = -Wl,-undefined,error
else
  NOUNDEF = -Wl,--no-undefined
endif


.PHONY : all install uninstall clean version

all : $(LIBS)

libspecex.a : $(OBJS)
	ar cr $@ $(OBJS)
	ranlib $@

libspecex.$(PLUG_EXT) : $(OBJS)
	$(CXX) $(NOUNDEF) $(PLUG_FLAGS) $(PLUG_LINK) -o $@ $(OBJS) $(LINK)

install : all
	mkdir -p $(SPECEX_PREFIX)/lib
	cp $(LIBS) $(SPECEX_PREFIX)/lib

uninstall :
	- $(RM) "$(SPECEX_PREFIX)/lib/libspecex*"

clean :
	- $(RM) $(LIBS) *.o *~

#
# Enable 'make version' to update the version string.
# Do make TAG=0.1.2 version to set the tag explicitly.
#
LASTTAG := $(shell git describe --tags --dirty --always | cut -d- -f1)
COUNT := $(shell git rev-list --count HEAD)
version :
	- $(RM) version.h
	@ if test -n "$(TAG)"; then v=$(TAG); else v=$(LASTTAG).dev$(COUNT); fi; \
		echo "#define VERSION_STRING \"$$v\"" > version.h

%.o : %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -I. -fPIC -o $@ -c $<

