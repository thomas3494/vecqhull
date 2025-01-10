# Installation

## Prerequisites

For building the library, you need [Highway](https://github.com/google/highway) and a C++ compiler.

For using the library you need a C compiler, standard library, and an OpenMP
implementation.

## Unix-like

For Unix and derivatives like Linux, macOS, \*BSD.

1. First install [Highway](https://github.com/google/highway)

2. Run `sudo make install` on your own system, or something like
`PREFIX=$HOME/.local make install` on a remote system you have no root
access for. See `config.mk` for more options.

## Windows

Pull requests welcome. For now, you can install WSL and follow the steps above.

# Usage

You can compile programs by linking against `-lvqhull -lm -fopenmp`.

The provided functions are documented in the man pages. The directory
`examples` contains some examples.
