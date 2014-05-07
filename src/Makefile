
DIRS = library plugin tests apps

.PHONY : all install clean uninstall $(DIRS)

all : $(DIRS)

apps : plugin library
	@cd apps; $(MAKE); cd ..

tests : library
	@cd tests; $(MAKE); cd ..

plugin : library
	@cd plugin; $(MAKE); cd ..

library :
	@cd library; $(MAKE); cd ..

install : all
	@for dir in $(DIRS); do \
	cd $${dir}; $(MAKE) install; \
	cd ..; \
	done

uninstall :
	@for dir in $(DIRS); do \
	cd $${dir}; $(MAKE) uninstall; \
	cd ..; \
	done

clean :
	@for dir in $(DIRS); do \
	cd $${dir}; $(MAKE) clean; \
	cd ..; \
	done

