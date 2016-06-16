
HEADERS := $(wildcard *.h)

OBJS := $(patsubst %.cc,%.o,$(wildcard *.cc))

DYEXT := so

LIBS = libspecex.a libspecex.$(DYEXT)

all : $(LIBS)

libspecex.a : $(OBJS)
	ar cr $@ $(OBJS)
	ranlib $@

libspecex.$(DYEXT) : $(OBJS)
	$(CXX) -Wl,--no-undefined -fPIC -shared -o $@ $(OBJS) $(LINK)

install : all
	@mkdir -p "$(SPECEX_PREFIX)/lib"; \
	cp $(LIBS) "$(SPECEX_PREFIX)/lib/"

uninstall :
	@rm -f "$(SPECEX_PREFIX)/lib/libspecex*"

clean :
	@rm -f $(LIBS) *.o *~


%.o : %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -fPIC -o $@ -c $<


