# Installation

## Prerequisites

For building the library, you need [Highway](https://github.com/google/highway) and a C++ compiler.

For using the library you need a C compiler, standard library, and an OpenMP
implementation.

## Unix-like

For Unix and derivatives like Linux, macOS, \*BSD.

1. First install 

2. Set the compiler and installation prefix in `config.mk`

3. Run `[sudo] make install`

## Windows

Pull requests welcome. For now, you can install WSL
and follow the steps above.

# Usage

You can compile programs by linking against `-lvqhull -lm -fopenmp`.

The provided functions are documented in the man pages. The directory
`examples` contains some examples.
