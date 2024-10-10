include config.mk

VERSION_MAJOR = 1
VERSION_MINOR = 0
VERSION_PATCH = 0

VERSION = $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

MAN3    = $(wildcard man/man3/*)

all: $(ANAME)

$(ANAME): fvqhull.o
	$(AR) -rc $@ $<
	$(RANLIB) $@

fvqhull.o: fvqhull.f90
	$(FC) -c $< -o $@

install: all
	mkdir -p "$(DESTDIR)$(LIBPREFIX)"
	cp -f $(ANAME) "$(DESTDIR)$(LIBPREFIX)"
	cp -f $(MODNAME) "$(DESTDIR)$(INCPREFIX)"

uninstall:
	$(RM) -f $(ANAME) "$(DESTDIR)$(LIBPREFIX)"
	$(RM) -f $(MODNAME) "$(DESTDIR)$(INCPREFIX)"

clean:
	$(RM) fvqhull.o $(ANAME) $(MODNAME)
