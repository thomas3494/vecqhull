include config.mk

VERSION_MAJOR = 1
VERSION_MINOR = 0
VERSION_PATCH = 0
MAN_DATE      = 2024-09-23

VERSION = $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

MAN3    = $(wildcard man/man3/*)

.PHONY: all clean install uninstall

all: $(ANAME) $(SONAME)

$(ANAME): vqhull.o
	$(AR) -rc $@ $<
	$(RANLIB) $@

$(SONAME): vqhull.o
	$(CXX) $< $(SOFLAGS) $(LDFLAGS) -o $@

vqhull.o: src/vqhull.cpp include/vqhull.h
	$(CXX) -Iinclude $(FLAGS) $(SHFLAGS) -c $< -o $@ $(LFLAGS)

install: all
	mkdir -p "$(DESTDIR)$(LIBPREFIX)"
	mkdir -p "$(DESTDIR)$(INCPREFIX)"
	mkdir -p "$(DESTDIR)$(MANPREFIX)/man3"
	cp -f $(MAN3) "$(DESTDIR)$(MANPREFIX)/man3"
	cp -f $(ANAME) "$(DESTDIR)$(LIBPREFIX)"
	cp -f $(SONAME) "$(DESTDIR)$(LIBPREFIX)/$(SONAME)"
	cp -f include/vqhull.h "$(DESTDIR)$(INCPREFIX)"

clean:
	$(RM) vqhull.o $(ANAME) $(SONAME)

uninstall:
	for m in $(MAN3); do $(RM) "$(DESTDIR)$(MANPREFIX)/man3/`basename $$m`"; done
	$(RM) -f "$(DESTDIR)$(LIBPREFIX)/$(ANAME)"
	$(RM) -f "$(DESTDIR)$(LIBPREFIX)/$(SONAME)"
	$(RM) -f "$(DESTDIR)$(INCPREFIX)/vqhull.h"
