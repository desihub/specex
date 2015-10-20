
HEADERS := $(wildcard *.h)

OBJS := $(patsubst %.cc,%.o,$(wildcard *.cc))


LIBS = libspecex.a

all : $(LIBS)

libspecex.a : $(OBJS)
	ar cru $@ $(OBJS)
	ranlib $@

install : all
	@mkdir -p "$(SPECEX_PREFIX)/lib"; \
	cp $(LIBS) "$(SPECEX_PREFIX)/lib/"

uninstall :
	@rm -f "$(SPECEX_PREFIX)/lib/libspecex*"

clean :
	@rm -f $(LIBS) *.o *~


%.o : %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<


