# SHTOOLS.jl

A Julia package wrapping
[SHTOOLS](https://github.com/SHTOOLS/SHTOOLS), the Spherical Harmonic
Tools.

* [![GitHub
  CI](https://github.com/eschnett/SHTOOLS.jl/workflows/CI/badge.svg)](https://github.com/eschnett/SHTOOLS.jl/actions)
* [![Codecov](https://codecov.io/gh/eschnett/SHTOOLS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/eschnett/SHTOOLS.jl)

Note that the function arguments are not always the same. For example,
output arguments to mutating functions come first in the argument
list, and are omitted in non-mutating functions.

## Functions implemented so far:

### Legendre Functions:

```Julia
PlmBar
PlmBar_d1
PlBar
PlBar_d1

PlmON
PlmON_d1
PlON
PlON_d1

PlmSchmidt
PlmSchmidt_d1
PlSchmidt
PlSchmidt_d1

PLegendreA
PLegendreA_d1
PLegendre
PLegendre_d1

PlmIndex
```

### Spherical harmonic transforms

```Julia
SHExpandDH
MakeGridDH
SHExpandDHC
MakeGridDHC
MakeGradientDH

SHGLQ
SHExpandGLQ
MakeGridGLQ
SHExpandGLQC
MakeGridGLQC
GLQGridCoord

SHExpandLSQ
MakeGrid2d
MakeGridPoint
MakeGridPointC
SHMultiply
```
