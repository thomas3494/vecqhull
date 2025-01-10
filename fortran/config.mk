# Customize below to fit your system

# paths
PREFIX    ?= $(HOME)/.local
INCPREFIX  = $(PREFIX)/include
LIBPREFIX  = $(PREFIX)/lib
MANPREFIX  = $(PREFIX)/share/man

# names
ANAME     = libfvqhull.a
MODNAME   = fvqhull.mod
#SONAME    = libfvqhull.so.$(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

FC      = g++
AR       = ar
RANLIB   = ranlib
