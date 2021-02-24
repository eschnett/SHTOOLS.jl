# SHTOOLS.jl

A Julia package wrapping
[SHTOOLS](https://github.com/SHTOOLS/SHTOOLS), the Spherical Harmonic
Tools.

[![Documenter](https://img.shields.io/badge/docs-dev-blue.svg)](https://eschnett.github.io/SHTOOLS.jl/dev)
* [![GitHub
  CI](https://github.com/eschnett/SHTOOLS.jl/workflows/CI/badge.svg)](https://github.com/eschnett/SHTOOLS.jl/actions)
* [![Codecov](https://codecov.io/gh/eschnett/SHTOOLS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/eschnett/SHTOOLS.jl)

Most functions come in two versions, one that mutates its arguments
and one that allocates its output.

Note that the function arguments are not always the same as in Fortran
or C. For example, output arguments for mutating functions come first
in the argument list, and are omitted in non-mutating functions.

The documentation lists the implemented functions as well as their
Julia signatures.
