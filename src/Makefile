
DIRS = library plugin tests apps

.PHONY : all install clean uninstall version $(DIRS)

all : $(DIRS)

apps : plugin library
	$(MAKE) -C apps

tests : library
	$(MAKE) -C tests

plugin : library
	$(MAKE) -C plugin

library :
	$(MAKE) -C library

install : all
	@ for d in $(DIRS); do $(MAKE) -C $$d install; done

uninstall :
	@ for d in $(DIRS); do $(MAKE) -C $$d uninstall; done

clean :
	@ for d in $(DIRS); do $(MAKE) -C $$d clean; done

version :
	$(MAKE) -C library version
