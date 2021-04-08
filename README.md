# SHTOOLS.jl

A Julia package wrapping
[SHTOOLS](https://github.com/SHTOOLS/SHTOOLS), the Spherical Harmonic
Tools.

* [![Documenter](https://img.shields.io/badge/docs-dev-blue.svg)](https://eschnett.github.io/SHTOOLS.jl/dev)
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

## Example: Calculate gradient of a scalar field

```Julia
julia> using SHTOOLS

julia> n = 4;

julia> nθ = n;

julia> nϕ = 2n;

julia> Θ = [π*(i-1)/nθ for i in 1:nθ];

julia> Φ = [2π*(j-1)/nϕ for j in 1:nϕ];

julia> # z + 2x
       F = [cos(θ) + 2*sin(θ)*cos(ϕ) for θ in Θ, ϕ in Φ];

julia> chop(x) = abs2(x) < 10eps(x) ? zero(x) : x;

julia> C,lmax = SHExpandDH(F, n; sampling=2);

julia> chop.(C[1,:,:])

julia> chop.(C[2,:,:])

julia> ∂θF, ∂ϕF, _ = MakeGradientDH(C, lmax; sampling=2);

julia> chop.(∂θF)

julia> chop.(∂ϕF)
```
These are wrong???
