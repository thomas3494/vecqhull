# Customize below to fit your system

# Paths
PREFIX    ?= $(HOME)/.local
INCPREFIX  = $(PREFIX)/include
LIBPREFIX  = $(PREFIX)/lib
MANPREFIX  = $(PREFIX)/share/man

# names
ANAME     = libvqhull.a
SONAME    = libvqhull.so.$(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

# Flags. For debugging, call BUILD=DEBUG make. We rely on IEEE-754 here, so
# do not use -ffast-math or -Ofast
BUILD      ?= RELEASE
FLAGS       = -O3 -march=native -mtune=native -Wall -Wextra -ffunction-sections -fdata-sections -fopenmp -nodefaultlibs -lgcc -fno-exceptions
DEBUG_FLAGS = -fsanitize=address -fsanitize=undefined -ggdb
RELEASE_FLAGS = -D NDEBUG
FLAGS += $(${BUILD}_FLAGS)
LDFLAGS  = -lm -fopenmp -lhwy -Wl,--gc-sections
LDFLAGS += $(${BUILD}_FLAGS)

SHFLAGS   = -fPIC
SOFLAGS   = -shared -nostdlib -Wl,--soname=libvqhull.so.$(VERSION_MAJOR).$(VERSION_MINOR)
SOSYMLINK = true

CXX      = g++
AR       = ar
RANLIB   = ranlib
