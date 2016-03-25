ifndef SPECEX_PREFIX
SPECEX_PREFIX = $(error SPECEX_PREFIX undefined)UNDEFINED
endif

# get the variables from harpconfig

export CXX := $(shell harpconfig --cxx)

export CXXFLAGS := $(shell harpconfig --cxxflags --cppflags) -I. -Wuninitialized -Wunused-value -Wunused-variable

export PLUG_FLAGS := $(shell harpconfig --plugflags)
export PLUG_LINK := $(shell harpconfig --pluglink)
export PLUG_EXT := $(shell harpconfig --plugext)

export LINK := $(shell harpconfig --link)

PYSCRIPTS = specex_mean_psf.py

# descend to src directory

.PHONY : all clean install uninstall

all : 
	@cd src/ ; $(MAKE)

install : all
	@cd src/ ; $(MAKE) install; \
	cd ../python; cp $(PYSCRIPTS) $(SPECEX_PREFIX)

uninstall :
	@cd src/ ; $(MAKE) uninstall

clean :
	@cd src/ ; $(MAKE) clean

