# Prerequisites

1. The C-library `vqhull` of the parent directory must be installed.

2. A Fortran compiler (e.g. gfortran)

3. A Modula-2 compiler (e.g. gm2)

# Installation

1. Set the compiler and installation prefix in `config.mk`

2. Run `[sudo] make install`

# Compiling

In `examples` contains a Makefile for a GNU system with some example programs.
There are some pitfalls.

1. The linking order matters, `-lfvqhull` must come before `-lvqhull`.

2. The search path to the `.mod` file must be specified explicitly with the `-J` option.
