export Wigner3j!
@doc raw"""
    w3j, jmin, jmax = Wigner3j!(w3j::AbstractVector{Cdouble}, j2::Integer, j3::Integer,
                                m1::Integer, m2::Integer, m3::Integer,
                                exitstatus::Union{Nothing,Ref{<:Integer}} = nothing)

Compute the Wigner 3j symbol

```math
\left(\begin{array}{ccc}
j & j_{2} & j_{3}\\
m_{1} & m_{2} & m_{3}
\end{array}\right)
```

for all integral values of `j` that satisfy `max(abs(m1), abs(j2 - j3)) <= j <= j2 + j3`.
The values of `jmin` and `jmax` are set to `max(abs(m1), abs(j2 - j3))` and `j2 + j3` respectively.
The vector `w3j` is filled with the computed 3j symbols.
Its first index corresponds to `j = jmin`,
and it has its first `jmax - jmin + 1` indices populated.

!!! note
    `w3j` must be sufficiently long to contain all the elements.
    Pre-existing values in `w3j` may be overwritten.

See also: [`Wigner3j`](@ref)

SHTOOLS page: [`Wigner3j`](https://shtools.oca.eu/shtools/public/wigner3j.html)
"""
function Wigner3j!(w3j::AbstractVector{Cdouble},
    j2::Integer, j3::Integer,
    m1::Integer, m2::Integer, m3::Integer,
    exitstatus::Optional{Ref{<:Integer}}=nothing)

    @assert j2 ≥ 0 "j2 must be >= 0"
    @assert j3 ≥ 0 "j2 must be >= 0"
    @assert abs(m2) ≤ j2 "m2 must satisfy $(-j2) <= m2 <= $j2 for j2 = $j2"
    @assert abs(m3) ≤ j3 "m3 must satisfy $(-j3) <= m3 <= $j3 for j3 = $j3"
    @assert abs(m1) ≤ j2 + j3 "m1 must satisfy $(-(j2 + j3)) <= m1 <= $(j2 + j3)"
    w3j_minlength = j2 + j3 - max(abs(m1), abs(j2 - j3)) + 1
    @assert length(w3j) ≥ w3j_minlength "length of w3j must be >= $w3j_minlength"

    exitstatus′ = Ref{Cint}()
    jmin = Ref{Cint}()
    jmax = Ref{Cint}()
    ccall((:Wigner3j, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Ref{Cint}, Ref{Cint},
            Cint, Cint, Cint, Cint, Cint, Ref{Cint}),
          w3j, jmin, jmax, j2, j3, m1, m2, m3, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 &&
            error("Wigner3j!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end

    return w3j, jmin[], jmax[]
end

export Wigner3j
@doc raw"""
    w3j, jmin, jmax = Wigner3j(j2::Integer, j3::Integer,
                                m1::Integer, m2::Integer, m3::Integer,
                                exitstatus::Union{Nothing,Ref{<:Integer}} = nothing)

Compute the Wigner 3j symbol

```math
\left(\begin{array}{ccc}
j & j_{2} & j_{3}\\
m_{1} & m_{2} & m_{3}
\end{array}\right)
```

for all integral values of `j` that satisfy `max(abs(m1), abs(j2 - j3)) <= j <= j2 + j3`.
The values of `jmin` and `jmax` are set to `max(abs(m1), abs(j2 - j3))` and `j2 + j3` respectively.
The vector `w3j` contains the computed 3j symbols. Its first index corresponds to `j = jmin`,
and last index corresponds to `j = jmax`.

See also: [`Wigner3j!`](@ref)

SHTOOLS page: [`Wigner3j`](https://shtools.oca.eu/shtools/public/wigner3j.html)
"""
function Wigner3j(j2::Integer, j3::Integer,
    m1::Integer, m2::Integer, m3::Integer, exitstatus::Optional{Ref{<:Integer}}=nothing)

    @assert j2 ≥ 0 "j2 must be >= 0"
    @assert j3 ≥ 0 "j2 must be >= 0"
    @assert abs(m1) ≤ j2 + j3 "m1 must satisfy $(-(j2 + j3)) <= m1 <= $(j2 + j3)"

    w3j_minlength = j2 + j3 - max(abs(m1), abs(j2 - j3)) + 1
    w3j = Vector{Cdouble}(undef, w3j_minlength)
    w3j, jmin, jmax = Wigner3j!(w3j, j2, j3, m1, m2, m3, exitstatus)
    return w3j, jmin, jmax
end
