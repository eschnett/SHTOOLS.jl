# SHTOOLS.jl

The [`SHTOOLS.jl`](https://github.com/eschnett/SHTOOLS.jl) package
wraps [SHTOOLS](https://shtools.github.io/SHTOOLS/index-fortran.html),
the Spherical Harmonic Tools.

*SHTOOLS is an archive of Fortran 95 software that can be used to
perform spherical harmonic transforms, multitaper spectral analyses,
expansions of functions into Slepian bases, and standard operations on
global gravitational and magnetic field data.*

## Legendre functions

### 4π normalized

```@docs
PlmBar!
PlmBar
PlmBar_d1!
PlmBar_d1
PlBar!
PlBar
PlBar_d1!
PlBar_d1
```

### Orthonormalized

```@docs
PlmON!
PlmON
PlmON_d1!
PlmON_d1
PlON!
PlON
PlON_d1!
PlON_d1
```

### Schmidt semi-normalized

```@docs
PlmSchmidt!
PlmSchmidt
PlmSchmidt_d1!
PlmSchmidt_d1
PlSchmidt!
PlSchmidt
PlSchmidt_d1!
PlSchmidt_d1
```

### Unnormalized

```@docs
PLegendreA!
PLegendreA
PLegendreA_d1!
PLegendreA_d1
PLegendre!
PLegendre
PLegendre_d1!
PLegendre_d1
```

### Utilities

```@docs
PlmIndex
```

## Spherical harmonic transforms

### Equally sampled (N×N) and equally spaced (N×2N) grids

```@docs
SHExpandDH
MakeGridDH
SHExpandDHC
MakeGridDHC
MakeGradientDH
```

### Gauss-Legendre quadrature grids 

```@docs
SHGLQ
SHExpandGLQ
MakeGridGLQ
SHExpandGLQC
MakeGridGLQC
GLQGridCoord
```

### Other routines

```@docs
SHExpandLSQ
MakeGrid2d
MakeGridPoint
MakeGridPointC
SHMultiply
```

## Spherical harmonic I/O, storage, and conversions

### Spherical harmonic storage

```@docs
SHCilmToVector!
SHCilmToVector
SHVectorToCilm!
SHVectorToCilm
```

### Spherical harmonic conversions
