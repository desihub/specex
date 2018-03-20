#
# Template Makefile for use with desiInstall.  You can assume that
# desiInstall will set these environment variables:
#
# WORKING_DIR   : The directory containing the svn export
# INSTALL_DIR   : The directory the installed product will live in.
# (PRODUCT)     : Where (PRODUCT) is replaced with the name of the
#                 product in upper case, e.g. DESITEMPLATE.  This should
#                 be the same as WORKING_DIR for typical installs.
#
# Use this shell to interpret shell commands, & pass its value to sub-make
#
export SHELL = /bin/sh
#
# This is like doing 'make -w' on the command line.  This tells make to
# print the directory it is in.
#
MAKEFLAGS = w
#
# This is a list of subdirectories that make should descend into.  Makefiles
# in these subdirectories should also understand 'make all' & 'make clean'.
# This list can be empty, but should still be defined.
#
SUBDIRS = src
#
# Set SPECEX_PREFIX, using INSTALL_DIR if available.
#
ifdef INSTALL_DIR
export SPECEX_PREFIX = $(INSTALL_DIR)
endif
ifndef SPECEX_PREFIX
NO_SPECEX_PREFIX = $(error SPECEX_PREFIX is undefined!)
endif
#
# Get the variables from harpconfig
#
export CXX := $(shell harpconfig --cxx)
ifndef CXX
NO_CXX = $(error CXX is undefined! Is harpconfig in your PATH?)
endif
export CXXFLAGS := $(shell harpconfig --cxxflags --cppflags) -I. -Wuninitialized -Wunused-value -Wunused-variable
export PLUG_FLAGS := $(shell harpconfig --plugflags)
export PLUG_LINK := $(shell harpconfig --pluglink)
export PLUG_EXT := $(shell harpconfig --plugext)
export LINK := $(shell harpconfig --link)
export CPPFLAGS := $(shell harpconfig --cppflags)
#
# Copy these scripts to the bin directory.
#
PYSCRIPTS = specex_mean_psf.py
#
# This is a message to make that these targets are 'actions' not files.
#
.PHONY : all clean install uninstall version
#
# This should compile all code prior to it being installed.
#
all : ; $(NO_CXX)
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f all ; done

install : all ; $(NO_SPECEX_PREFIX)
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f install; done
	@ for f in $(PYSCRIPTS); do cp python/$$f $(SPECEX_PREFIX)/bin; done
	@ cp -a data $(SPECEX_PREFIX)/
	@ chmod +x $(SPECEX_PREFIX)/bin/specex*

uninstall : ; $(NO_SPECEX_PREFIX)
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f uninstall; done
#
# GNU make pre-defines $(RM).  The - in front of $(RM) causes make to
# ignore any errors produced by $(RM).
#
clean :
	- $(RM) *~ core
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f clean ; done
#
# Enable 'make version' to update the version string.
# Do make TAG=0.1.2 version to set the tag explicitly.
#
version :
	$(MAKE) -C src version
