# Installation

## Prerequisites

[Highway](https://github.com/google/highway) and a C++ compiler.


## Unix-like

For Unix and derivatives like Linux, macOS, \*BSD.

1. First install 

2. Set the compiler and installation prefix in `config.mk`

3. Run `[sudo] make install`

## Windows

Pull requests welcome. For now, you can install WSL
and follow the steps above.

# Usage

You can compile programs by linking against `-lvqhull -lm -fopenmp -libstd++`.
TODO: can we kill additional dependencies? I don't think we use anything from
the C++ standard library.

The provided functions are documented in the man pages. The directory
`examples` contains some examples.
