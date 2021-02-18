module SHTOOLS

using SHTOOLS_jll

const Optional{T} = Union{Nothing,T}
# For optional arguments
optional(x::Optional, x0) = x !== nothing ? x : x0

################################################################################

# Legendre functions

# 4π normalized
# Orthonormalized
# Schmidt semi-normalized
# Unnormalized
# Utilities

for op in [:PlmBar, :PlmON, :PlmSchmidt]
    op! = Symbol(op, "!")
    op_d1 = Symbol(op, "_d1")
    op_d1! = Symbol(op, "_d1!")
    @eval begin
        #
        export $op!
        function $op!(p::AbstractVector{Cdouble}, lmax::Integer,
                      z::Union{AbstractFloat,Integer}; csphase::Integer=1,
                      cnorm::Integer=0,
                      exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            @assert csphase ∈ (1, -1)
            @assert cnorm ∈ (0, 1)
            @assert length(p) ≥ (lmax + 1) * (lmax + 2) ÷ 2
            exitstatus′ = Ref{Cint}()
            ccall(($(QuoteNode(op)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Cint, Cdouble, Ref{Cint}, Ref{Cint},
                   Ref{Cint}), p, lmax, z, csphase, cnorm, exitstatus′)
            if exitstatus === nothing
                exitstatus′[] ≠ 0 &&
                    error("$($op!): Error code $(exitstatus′[])")
            else
                exitstatus[] = exitstatus′[]
            end
            return p
        end
        #
        export $op
        function $op(lmax::Integer, z::Union{AbstractFloat,Integer};
                     csphase::Integer=1, cnorm::Integer=0,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            p = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            $op!(p, lmax, z; csphase=csphase, cnorm=cnorm,
                 exitstatus=exitstatus)
            return p
        end
        #
        export $op_d1!
        function $op_d1!(p::AbstractVector{Cdouble},
                         dp::AbstractVector{Cdouble}, lmax::Integer,
                         z::Union{AbstractFloat,Integer}; csphase::Integer=1,
                         cnorm::Integer=0,
                         exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            @assert csphase ∈ (1, -1)
            @assert cnorm ∈ (0, 1)
            @assert length(p) ≥ lmax + 1
            @assert length(dp) ≥ lmax + 1
            exitstatus′ = Ref{Cint}()
            ccall(($(QuoteNode(op_d1)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Ref{Cint},
                   Ref{Cint}, Ref{Cint}), p, dp, lmax, z, csphase, cnorm,
                  exitstatus′)
            if exitstatus === nothing
                exitstatus′[] ≠ 0 &&
                    error("$($op_d1!): Error code $(exitstatus′[])")
            else
                exitstatus[] = exitstatus′[]
            end
            return p, dp
        end
        #
        export $op_d1
        function $op_d1(lmax::Integer, z::Union{AbstractFloat,Integer};
                        csphase::Integer=1, cnorm::Integer=0,
                        exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            p = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            dp = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            $op_d1!(p, dp, lmax, z; csphase=csphase, cnorm=cnorm,
                    exitstatus=exitstatus)
            return p, dp
        end
        #
    end
end

for op in [:PLegendreA]
    op! = Symbol(op, "!")
    op_d1 = Symbol(op, "_d1")
    op_d1! = Symbol(op, "_d1!")
    @eval begin
        #
        export $op!
        function $op!(p::AbstractVector{Cdouble}, lmax::Integer,
                      z::Union{AbstractFloat,Integer}; csphase::Integer=1,
                      exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            @assert csphase ∈ (1, -1)
            @assert length(p) ≥ (lmax + 1) * (lmax + 2) ÷ 2
            exitstatus′ = Ref{Cint}()
            ccall(($(QuoteNode(op)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Cint, Cdouble, Ref{Cint}, Ref{Cint}), p, lmax,
                  z, csphase, exitstatus′)
            if exitstatus === nothing
                exitstatus′[] ≠ 0 &&
                    error("$($op!): Error code $(exitstatus′[])")
            else
                exitstatus[] = exitstatus′[]
            end
            return p
        end
        #
        export $op
        function $op(lmax::Integer, z::Union{AbstractFloat,Integer};
                     csphase::Integer=1,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            p = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            $op!(p, lmax, z; csphase=csphase, exitstatus=exitstatus)
            return p
        end
        #
        export $op_d1!
        function $op_d1!(p::AbstractVector{Cdouble},
                         dp::AbstractVector{Cdouble}, lmax::Integer,
                         z::Union{AbstractFloat,Integer}; csphase::Integer=1,
                         exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            @assert csphase ∈ (1, -1)
            @assert length(p) ≥ lmax + 1
            @assert length(dp) ≥ lmax + 1
            exitstatus′ = Ref{Cint}()
            ccall(($(QuoteNode(op_d1)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Ref{Cint},
                   Ref{Cint}), p, dp, lmax, z, csphase, exitstatus′)
            if exitstatus === nothing
                exitstatus′[] ≠ 0 &&
                    error("$($op_d1!): Error code $(exitstatus′[])")
            else
                exitstatus[] = exitstatus′[]
            end
            return p, dp
        end
        #
        export $op_d1
        function $op_d1(lmax::Integer, z::Union{AbstractFloat,Integer};
                        csphase::Integer=1,
                        exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            p = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            dp = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            $op_d1!(p, dp, lmax, z; csphase=csphase, exitstatus=exitstatus)
            return p, dp
        end
        #
    end
end

for op in [:PlBar, :PlON, :PlSchmidt, :PLegendre]
    op! = Symbol(op, "!")
    op_d1 = Symbol(op, "_d1")
    op_d1! = Symbol(op, "_d1!")
    @eval begin
        #
        export $op!
        function $op!(p::AbstractVector{Cdouble}, lmax::Integer,
                      z::Union{AbstractFloat,Integer};
                      exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            @assert length(p) ≥ lmax + 1
            exitstatus′ = Ref{Cint}()
            ccall(($(QuoteNode(op)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Cint, Cdouble, Ref{Cint}), p, lmax, z,
                  exitstatus′)
            if exitstatus === nothing
                exitstatus′[] ≠ 0 &&
                    error("$($op!): Error code $(exitstatus′[])")
            else
                exitstatus[] = exitstatus′[]
            end
            return p
        end
        #
        export $op
        function $op(lmax::Integer, z::Union{AbstractFloat,Integer};
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            p = Array{Cdouble}(undef, lmax + 1)
            $op!(p, lmax, z; exitstatus=exitstatus)
            return p
        end
        #
        export $op_d1!
        function $op_d1!(p::AbstractVector{Cdouble},
                         dp::AbstractVector{Cdouble}, lmax::Integer,
                         z::Union{AbstractFloat,Integer};
                         exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            @assert length(p) ≥ lmax + 1
            @assert length(dp) ≥ lmax + 1
            exitstatus′ = Ref{Cint}()
            ccall(($(QuoteNode(op_d1)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Ref{Cint}), p, dp,
                  lmax, z, exitstatus′)
            if exitstatus === nothing
                exitstatus′[] ≠ 0 &&
                    error("$($op_d1!): Error code $(exitstatus′[])")
            else
                exitstatus[] = exitstatus′[]
            end
            return p, dp
        end
        #
        export $op_d1
        function $op_d1(lmax::Integer, z::Union{AbstractFloat,Integer};
                        exitstatus::Optional{Ref{<:Integer}}=nothing)
            @assert lmax ≥ 0
            p = Array{Cdouble}(undef, lmax + 1)
            dp = Array{Cdouble}(undef, lmax + 1)
            $op_d1!(p, dp, lmax, z; exitstatus=exitstatus)
            return p, dp
        end
        #
    end
end

export PlmIndex
function PlmIndex(l::Integer, m::Integer)
    0 ≤ m ≤ l ||
        error("m must be greater than or equal to zero and less than or equal to l: m=$m, l=$l")
    index = ccall((:PlmIndex, libSHTOOLS), Cint, (Cint, Cint), l, m)
    return Int(index)
end

################################################################################

# Spherical harmonic transforms

# Equally sampled (N×N) and equally spaced (N×2N) grids

export SHExpandDH!
function SHExpandDH!(cilm::AbstractArray{Cdouble,3},
                     griddh::AbstractArray{Cdouble,2}, n::Integer;
                     norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                     lmax_calc::Optional{Integer}=nothing,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert n > 0 && iseven(n)
    @assert norm ∈ (1, 2, 3, 4)
    @assert sampling ∈ (1, 2)
    @assert csphase ∈ (1, -1)
    lmax_calc′ = optional(lmax_calc, n ÷ 2 - 1)
    @assert lmax_calc′ ≥ 0
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax_calc′ + 1
    @assert size(cilm, 3) == size(cilm, 2)
    @assert size(griddh, 1) ≥ n
    @assert size(griddh, 2) ≥ sampling * n
    lmax′ = Ref{Cint}()
    exitstatus′ = Ref{Cint}()
    ccall((:SHExpandDH, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}, Cint, Ref{Cint},
           Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}), griddh,
          size(griddh, 1), size(griddh, 2), n, cilm, size(cilm, 2), lmax′, norm,
          sampling, csphase, lmax_calc′, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("SHExpandDH!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return cilm, Int(lmax′[])
end

export SHExpandDH
function SHExpandDH(griddh::AbstractArray{Cdouble,2}, n::Integer;
                    norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                    lmax_calc::Optional{Integer}=nothing,
                    exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert n > 0 && iseven(n)
    lmax_calc′ = optional(lmax_calc, n ÷ 2 - 1)
    cilm = Array{Cdouble}(undef, 2, lmax_calc′ + 1, lmax_calc′ + 1)
    _, lmax = SHExpandDH!(cilm, griddh, n; norm=norm, sampling=sampling,
                          csphase=csphase, lmax_calc=lmax_calc,
                          exitstatus=exitstatus)
    return cilm, lmax
end

export MakeGridDH!
function MakeGridDH!(griddh::AbstractArray{Cdouble,2},
                     cilm::AbstractArray{Cdouble,3}, lmax::Integer;
                     norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                     lmax_calc::Optional{Integer}=nothing, extend::Integer=0,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert norm ∈ (1, 2, 3, 4)
    @assert sampling ∈ (1, 2)
    @assert csphase ∈ (1, -1)
    @assert extend ∈ (0, 1)
    n′ = Int(2 * lmax + 2)
    @assert size(griddh, 1) ≥ n′ + extend
    @assert size(griddh, 2) ≥ sampling * n′ + extend
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax + 1
    @assert size(cilm, 3) == size(cilm, 2)
    lmax_calc′ = optional(lmax_calc, lmax)
    exitstatus′ = Ref{Cint}()
    ccall((:MakeGridDH, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Cint, Ref{Cint}, Ptr{Cdouble}, Cint, Cint,
           Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}),
          griddh, size(griddh, 1), size(griddh, 2), n′, cilm, size(cilm, 2),
          lmax, norm, sampling, csphase, lmax_calc′, extend, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("MakeGridDH!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return griddh, n′
end

export MakeGridDH
function MakeGridDH(cilm::AbstractArray{Cdouble,3}, lmax::Integer;
                    norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                    lmax_calc::Optional{Integer}=nothing, extend::Integer=0,
                    exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert sampling ∈ (1, 2)
    @assert extend ∈ (0, 1)
    n′ = 2 * lmax + 2
    griddh = Array{Cdouble}(undef, n′ + extend, sampling * n′ + extend)
    _, n = MakeGridDH!(griddh, cilm, lmax; norm=norm, sampling=sampling,
                       csphase=csphase, lmax_calc=lmax_calc, extend=extend,
                       exitstatus=exitstatus)
    return griddh, n
end

function SHExpandDH!(cilm::AbstractArray{Complex{Cdouble},3},
                     griddh::AbstractArray{Complex{Cdouble},2}, n::Integer;
                     norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                     lmax_calc::Optional{Integer}=nothing,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert n > 0 && iseven(n)
    @assert norm ∈ (1, 2, 3, 4)
    @assert sampling ∈ (1, 2)
    @assert csphase ∈ (1, -1)
    lmax_calc′ = optional(lmax_calc, n ÷ 2 - 1)
    @assert lmax_calc′ ≥ 0
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax_calc′ + 1
    @assert size(cilm, 3) == size(cilm, 2)
    @assert size(griddh, 1) ≥ n
    @assert size(griddh, 2) ≥ sampling * n
    lmax′ = Ref{Cint}()
    exitstatus′ = Ref{Cint}()
    ccall((:SHExpandDHC, libSHTOOLS), Cvoid,
          (Ptr{Complex{Cdouble}}, Cint, Cint, Cint, Ptr{Complex{Cdouble}}, Cint,
           Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}),
          griddh, size(griddh, 1), size(griddh, 2), n, cilm, size(cilm, 2),
          lmax′, norm, sampling, csphase, lmax_calc′, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("SHExpandDHC!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return cilm, Int(lmax′[])
end

function SHExpandDH(griddh::AbstractArray{Complex{Cdouble},2}, n::Integer;
                    norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                    lmax_calc::Optional{Integer}=nothing,
                    exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert n > 0 && iseven(n)
    lmax_calc′ = optional(lmax_calc, n ÷ 2 - 1)
    cilm = Array{Complex{Cdouble}}(undef, 2, lmax_calc′ + 1, lmax_calc′ + 1)
    _, lmax = SHExpandDH!(cilm, griddh, n; norm=norm, sampling=sampling,
                          csphase=csphase, lmax_calc=lmax_calc,
                          exitstatus=exitstatus)
    return cilm, lmax
end

export SHExpandDHC!
const SHExpandDHC! = SHExpandDH!
export SHExpandDHC
const SHExpandDHC = SHExpandDH

function MakeGridDH!(griddh::AbstractArray{Complex{Cdouble},2},
                     cilm::AbstractArray{Complex{Cdouble},3}, lmax::Integer;
                     norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                     lmax_calc::Optional{Integer}=nothing, extend::Integer=0,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert norm ∈ (1, 2, 3, 4)
    @assert sampling ∈ (1, 2)
    @assert csphase ∈ (1, -1)
    @assert extend ∈ (0, 1)
    n′ = Int(2 * lmax + 2)
    @assert size(griddh, 1) ≥ n′ + extend
    @assert size(griddh, 2) ≥ sampling * n′ + extend
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax + 1
    @assert size(cilm, 3) == size(cilm, 2)
    lmax_calc′ = optional(lmax_calc, lmax)
    exitstatus′ = Ref{Cint}()
    ccall((:MakeGridDHC, libSHTOOLS), Cvoid,
          (Ptr{Complex{Cdouble}}, Cint, Cint, Ref{Cint}, Ptr{Complex{Cdouble}},
           Cint, Cint, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint},
           Ref{Cint}), griddh, size(griddh, 1), size(griddh, 2), n′, cilm,
          size(cilm, 2), lmax, norm, sampling, csphase, lmax_calc′, extend,
          exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("MakeGridDHC!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return griddh, n′
end

function MakeGridDH(cilm::AbstractArray{Complex{Cdouble},3}, lmax::Integer;
                    norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                    lmax_calc::Optional{Integer}=nothing, extend::Integer=0,
                    exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert sampling ∈ (1, 2)
    @assert extend ∈ (0, 1)
    n′ = 2 * lmax + 2
    griddh = Array{Complex{Cdouble}}(undef, n′ + extend, sampling * n′ + extend)
    _, n = MakeGridDH!(griddh, cilm, lmax; norm=norm, sampling=sampling,
                       csphase=csphase, lmax_calc=lmax_calc, extend=extend,
                       exitstatus=exitstatus)
    return griddh, n
end

export MakeGridDHC!
const MakeGridDHC! = MakeGridDH!
export MakeGridDHC
const MakeGridDHC = MakeGridDH

# TODO: Add MakeGradientDH once there is a C wrapper

export SHGLQ!
function SHGLQ!(zero::AbstractVector{Cdouble}, w::AbstractVector{Cdouble},
                plx::Optional{AbstractArray{Cdouble,2}}, lmax::Integer;
                norm::Integer=1, csphase::Integer=1, cnorm::Integer=0,
                exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert length(zero) ≥ lmax + 1
    @assert length(w) == length(zero)
    @assert lmax ≥ 0
    if plx !== nothing
        @assert size(plx) == (lmax + 1, (lmax + 1) * (lmax + 2) ÷ 2)
    end
    @assert norm ∈ (1, 2, 3, 4)
    @assert csphase ∈ (1, -1)
    @assert cnorm ∈ (0, 1)
    exitstatus′ = Ref{Cint}()
    ccall((:SHGLQ, libSHTOOLS), Cvoid,
          (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}, Ref{Cint},
           Ref{Cint}, Ref{Cint}), lmax, zero, w, optional(plx, Ptr{Cdouble}()),
          norm, csphase, cnorm, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("SHGLQ!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return zero, w
end

# Gauss-Legendre quadrature grids

export SHGLQ
function SHGLQ(plx::Optional{AbstractArray{Cdouble,2}}, lmax::Integer;
               norm::Integer=1, csphase::Integer=1, cnorm::Integer=0,
               exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    zero = Array{Cdouble}(undef, lmax + 1)
    w = Array{Cdouble}(undef, lmax + 1)
    if plx !== nothing
        @assert size(plx) == (lmax + 1, (lmax + 1) * (lmax + 2) ÷ 2)
    end
    SHGLQ!(zero, w, plx, lmax; norm=norm, csphase=csphase, cnorm=cnorm,
           exitstatus=exitstatus)
    return zero, w
end

export SHExpandGLQ!
function SHExpandGLQ!(cilm::AbstractArray{Cdouble,3}, lmax::Integer,
                      gridglq::AbstractArray{Cdouble,2},
                      w::AbstractVector{Cdouble},
                      plx::Optional{AbstractArray{Cdouble,2}},
                      zero::Optional{AbstractVector{Cdouble}}; norm::Integer=1,
                      csphase::Integer=1, lmax_calc::Optional{Integer}=nothing,
                      exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    lmax_calc′ = optional(lmax_calc, lmax)
    @assert lmax_calc′ ≥ 0
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax_calc′ + 1
    @assert size(cilm, 3) == size(cilm, 2)
    @assert size(gridglq) == (lmax + 1, 2 * lmax + 1)
    @assert length(w) == lmax + 1
    @assert (plx !== nothing) + (zero !== nothing) == 1
    if plx !== nothing
        @assert size(plx) == (lmax + 1, (lmax + 1) * (lmax + 2) ÷ 2)
    end
    if zero !== nothing
        @assert length(zero) == lmax + 1
    end
    @assert norm ∈ (1, 2, 3, 4)
    @assert csphase ∈ (1, -1)
    exitstatus′ = Ref{Cint}()
    ccall((:SHExpandGLQ, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}), cilm,
          size(cilm, 2), lmax, gridglq, w, optional(plx, Ptr{Cdouble}()),
          optional(zero, Ptr{Cdouble}()), norm, csphase, lmax_calc′,
          exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("SHExpandGLQ!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return cilm
end

export SHExpandGLQ
function SHExpandGLQ(lmax::Integer, gridglq::AbstractArray{Cdouble,2},
                     w::AbstractVector{Cdouble},
                     plx::Optional{AbstractArray{Cdouble,2}},
                     zero::Optional{AbstractVector{Cdouble}}; norm::Integer=1,
                     csphase::Integer=1, lmax_calc::Optional{Integer}=nothing,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    lmax_calc′ = optional(lmax_calc, lmax)
    @assert lmax_calc′ ≥ 0
    cilm = Array{Cdouble}(undef, 2, lmax_calc′ + 1, lmax_calc′ + 1)
    SHExpandGLQ!(cilm, lmax, gridglq, w, plx, zero; norm=norm, csphase=csphase,
                 lmax_calc=lmax_calc, exitstatus=exitstatus)
    return cilm
end

export MakeGridGLQ!
function MakeGridGLQ!(gridglq::AbstractArray{Cdouble,2},
                      cilm::AbstractArray{Cdouble,3}, lmax::Integer,
                      plx::Optional{AbstractArray{Cdouble,2}},
                      zero::Optional{AbstractVector{Cdouble}}; norm::Integer=1,
                      csphase::Integer=1, lmax_calc::Optional{Integer}=nothing,
                      extend::Integer=0,
                      exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert extend ∈ (0, 1)
    @assert size(gridglq) == (lmax + 1, 2 * lmax + 1 + extend)
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ 1
    @assert size(cilm, 3) == size(cilm, 2)
    @assert (plx !== nothing) + (zero !== nothing) == 1
    if plx !== nothing
        @assert size(plx) == (lmax + 1, (lmax + 1) * (lmax + 2) ÷ 2)
    end
    if zero !== nothing
        @assert length(zero) == lmax + 1
    end
    @assert norm ∈ (1, 2, 3, 4)
    @assert csphase ∈ (1, -1)
    lmax_calc′ = optional(lmax_calc, lmax)
    exitstatus′ = Ref{Cint}()
    ccall((:MakeGridGLQ, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble},
           Ptr{Cdouble}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}),
          gridglq, size(gridglq, 1), size(gridglq, 2), cilm, size(cilm, 2),
          lmax, optional(plx, Ptr{Cdouble}()), optional(zero, Ptr{Cdouble}()),
          norm, csphase, lmax_calc′, extend, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("MakeGridGLQ!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return gridglq
end

export MakeGridGLQ
function MakeGridGLQ(cilm::AbstractArray{Cdouble,3}, lmax::Integer,
                     plx::Optional{AbstractArray{Cdouble,2}},
                     zero::Optional{AbstractVector{Cdouble}}; norm::Integer=1,
                     csphase::Integer=1, lmax_calc::Optional{Integer}=nothing,
                     extend::Integer=0,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert extend ∈ (0, 1)
    gridglq = Array{Cdouble}(undef, lmax + 1, 2 * lmax + 1 + extend)
    MakeGridGLQ!(gridglq, cilm, lmax, plx, zero; norm=norm, csphase=csphase,
                 lmax_calc=lmax_calc, extend=extend, exitstatus=exitstatus)
    return gridglq
end

function SHExpandGLQ!(cilm::AbstractArray{Complex{Cdouble},3}, lmax::Integer,
                      gridglq::AbstractArray{Complex{Cdouble},2},
                      w::AbstractVector{Cdouble},
                      plx::Optional{AbstractArray{Cdouble,2}},
                      zero::Optional{AbstractVector{Cdouble}}; norm::Integer=1,
                      csphase::Integer=1, lmax_calc::Optional{Integer}=nothing,
                      exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    lmax_calc′ = optional(lmax_calc, lmax)
    @assert lmax_calc′ ≥ 0
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax_calc′ + 1
    @assert size(cilm, 3) == size(cilm, 2)
    @assert size(gridglq) == (lmax + 1, 2 * lmax + 1)
    @assert length(w) == lmax + 1
    @assert (plx !== nothing) + (zero !== nothing) == 1
    if plx !== nothing
        @assert size(plx) == (lmax + 1, (lmax + 1) * (lmax + 2) ÷ 2)
    end
    if zero !== nothing
        @assert length(zero) == lmax + 1
    end
    @assert norm ∈ (1, 2, 3, 4)
    @assert csphase ∈ (1, -1)
    exitstatus′ = Ref{Cint}()
    ccall((:SHExpandGLQC, libSHTOOLS), Cvoid,
          (Ptr{Complex{Cdouble}}, Cint, Cint, Ptr{Complex{Cdouble}},
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}, Ref{Cint},
           Ref{Cint}, Ref{Cint}), cilm, size(cilm, 2), lmax, gridglq, w,
          optional(plx, Ptr{Cdouble}()), optional(zero, Ptr{Cdouble}()), norm,
          csphase, lmax_calc′, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("SHExpandGLQC!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return cilm
end

function SHExpandGLQ(lmax::Integer, gridglq::AbstractArray{Complex{Cdouble},2},
                     w::AbstractVector{Cdouble},
                     plx::Optional{AbstractArray{Cdouble,2}},
                     zero::Optional{AbstractVector{Cdouble}}; norm::Integer=1,
                     csphase::Integer=1, lmax_calc::Optional{Integer}=nothing,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    lmax_calc′ = optional(lmax_calc, lmax)
    @assert lmax_calc′ ≥ 0
    cilm = Array{Complex{Cdouble}}(undef, 2, lmax_calc′ + 1, lmax_calc′ + 1)
    SHExpandGLQ!(cilm, lmax, gridglq, w, plx, zero; norm=norm, csphase=csphase,
                 lmax_calc=lmax_calc, exitstatus=exitstatus)
    return cilm
end

export SHExpandGLQC!
const SHExpandGLQC! = SHExpandGLQ!
export SHExpandGLQC
const SHExpandGLQC = SHExpandGLQ

function MakeGridGLQ!(gridglq::AbstractArray{Complex{Cdouble},2},
                      cilm::AbstractArray{Complex{Cdouble},3}, lmax::Integer,
                      plx::Optional{AbstractArray{Cdouble,2}},
                      zero::Optional{AbstractVector{Cdouble}}; norm::Integer=1,
                      csphase::Integer=1, lmax_calc::Optional{Integer}=nothing,
                      extend::Integer=0,
                      exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert extend ∈ (0, 1)
    @assert size(gridglq) == (lmax + 1, 2 * lmax + 1 + extend)
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ 1
    @assert size(cilm, 3) == size(cilm, 2)
    @assert (plx !== nothing) + (zero !== nothing) == 1
    if plx !== nothing
        @assert size(plx) == (lmax + 1, (lmax + 1) * (lmax + 2) ÷ 2)
    end
    if zero !== nothing
        @assert length(zero) == lmax + 1
    end
    @assert norm ∈ (1, 2, 3, 4)
    @assert csphase ∈ (1, -1)
    lmax_calc′ = optional(lmax_calc, lmax)
    exitstatus′ = Ref{Cint}()
    ccall((:MakeGridGLQC, libSHTOOLS), Cvoid,
          (Ptr{Complex{Cdouble}}, Cint, Cint, Ptr{Complex{Cdouble}}, Cint, Cint,
           Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}, Ref{Cint}, Ref{Cint},
           Ref{Cint}, Ref{Cint}), gridglq, size(gridglq, 1), size(gridglq, 2),
          cilm, size(cilm, 2), lmax, optional(plx, Ptr{Cdouble}()),
          optional(zero, Ptr{Cdouble}()), norm, csphase, lmax_calc′, extend,
          exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("MakeGridGLQC!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return gridglq
end

function MakeGridGLQ(cilm::AbstractArray{Complex{Cdouble},3}, lmax::Integer,
                     plx::Optional{AbstractArray{Cdouble,2}},
                     zero::Optional{AbstractVector{Cdouble}}; norm::Integer=1,
                     csphase::Integer=1, lmax_calc::Optional{Integer}=nothing,
                     extend::Integer=0,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert extend ∈ (0, 1)
    gridglq = Array{Complex{Cdouble}}(undef, lmax + 1, 2 * lmax + 1 + extend)
    MakeGridGLQ!(gridglq, cilm, lmax, plx, zero; norm=norm, csphase=csphase,
                 lmax_calc=lmax_calc, extend=extend, exitstatus=exitstatus)
    return gridglq
end

export MakeGridGLQC!
const MakeGridGLQC! = MakeGridGLQ!
export MakeGridGLQC
const MakeGridGLQC = MakeGridGLQ

export GLQGridCoord!
function GLQGridCoord!(latglq::AbstractVector{Cdouble},
                       longlq::AbstractVector{Cdouble}, lmax::Integer;
                       extend::Integer=0,
                       exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert extend ∈ (0, 1)
    @assert length(latglq) ≥ lmax + 1
    @assert length(longlq) ≥ 2 * lmax + 1 + extend
    nlat′ = Ref{Cint}()
    nlong′ = Ref{Cint}()
    exitstatus′ = Ref{Cint}()
    ccall((:GLQGridCoord, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint, Ref{Cint}, Ref{Cint},
           Ref{Cint}, Ref{Cint}), latglq, length(latglq), longlq,
          length(longlq), lmax, nlat′, nlong′, extend, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("GLQGridCoord!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return latglq, longlq, Int(nlat′[]), Int(nlong′[])
end

export GLQGridCoord
function GLQGridCoord(lmax::Integer; extend::Integer=0,
                      exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert extend ∈ (0, 1)
    latglq = Array{Cdouble}(undef, lmax + 1)
    longlq = Array{Cdouble}(undef, 2 * lmax + 1 + extend)
    GLQGridCoord!(latglq, longlq, lmax; extend=extend, exitstatus=exitstatus)
    return latglq, longlq
end

# Other routines

export SHExpandLSQ!
function SHExpandLSQ!(cilm::AbstractArray{Cdouble,3},
                      d::AbstractVector{Cdouble}, lat::AbstractVector{Cdouble},
                      lon::AbstractVector{Cdouble}, nmax::Integer,
                      lmax::Integer; norm::Integer=1, csphase::Integer=1,
                      weights::Vector{Cdouble},
                      exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax + 1
    @assert size(cilm, 3) == size(cilm, 2)
    @assert nmax ≥ 0
    @assert nmax ≥ (lmax + 1)^2
    @assert length(d) ≥ nmax
    @assert length(lat) ≥ nmax
    @assert length(lon) ≥ nmax
    @assert norm ∈ (1, 2, 3, 4)
    @assert length(weights) == nmax
    chi2′ = Ref{Cdouble}()
    exitstatus′ = Ref{Cint}()
    ccall((:SHExpandLSQ, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint,
           Cint, Ref{Cint}, Ref{Cdouble}, Ref{Cint}, Ptr{Cdouble}, Ref{Cint}),
          cilm, size(cilm, 2), d, lat, lon, nmax, lmax, norm, chi2′, csphase,
          weights, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("SHExpandLSQ!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return cilm, Cdouble(chi2′[])
end

export SHExpandLSQ
function SHExpandLSQ(d::AbstractVector{Cdouble}, lat::AbstractVector{Cdouble},
                     lon::AbstractVector{Cdouble}, nmax::Integer, lmax::Integer;
                     norm::Integer=1, csphase::Integer=1,
                     weights::Vector{Cdouble},
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    cilm = Array{Cdouble}(undef, 2, lmax + 1, lmax + 1)
    _, chi2 = SHExpandLSQ!(cilm, d, lat, lon, nmax, lmax; norm=norm,
                           csphase=csphase, weights=weights,
                           exitstatus=exitstatus)
    return cilm, chi2
end

export MakeGrid2d!
function MakeGrid2d!(grid::AbstractArray{Cdouble,2},
                     cilm::AbstractArray{Cdouble,3}, lmax::Integer,
                     interval::Union{AbstractFloat,Integer}; norm::Integer=1,
                     csphase::Integer=1,
                     f::Optional{Union{AbstractFloat,Integer}}=nothing,
                     a::Optional{Union{AbstractFloat,Integer}}=nothing,
                     north::Union{AbstractFloat,Integer}=90.0,
                     south::Union{AbstractFloat,Integer}=-90.0,
                     east::Union{AbstractFloat,Integer}=360.0,
                     west::Union{AbstractFloat,Integer}=0.0,
                     dealloc::Bool=false,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert interval > 0
    @assert size(grid, 1) ≥ floor(Int, 180 / interval + 1)
    @assert size(grid, 2) ≥ floor(Int, 360 / interval + 1)
    @assert size(cilm, 1) == 2
    @assert lmax ≥ 0
    @assert size(cilm, 2) ≥ lmax + 1
    @assert size(cilm, 3) == size(cilm, 2)
    @assert (f !== nothing) == (a !== nothing)
    nlat = Ref{Cint}()
    nlong = Ref{Cint}()
    exitstatus′ = Ref{Cint}()
    ccall((:MakeGrid2d, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Cint, Cint, Cdouble,
           Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cdouble},
           Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble},
           Ref{Cint}, Ref{Cint}), grid, size(grid, 1), size(grid, 2), cilm,
          size(cilm, 2), lmax, interval, nlat, nlong, norm, csphase,
          optional(f, Ptr{Cdouble}()), optional(a, Ptr{Cdouble}()), north,
          south, east, west, dealloc, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("MakeGrid2d!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return grid, Int(nlat[]), Int(nlong[])
end

export MakeGrid2d
function MakeGrid2d(cilm::AbstractArray{Cdouble,3}, lmax::Integer,
                    interval::Union{AbstractFloat,Integer}; norm::Integer=1,
                    csphase::Integer=1,
                    f::Optional{Union{AbstractFloat,Integer}}=nothing,
                    a::Optional{Union{AbstractFloat,Integer}}=nothing,
                    north::Union{AbstractFloat,Integer}=90.0,
                    south::Union{AbstractFloat,Integer}=-90.0,
                    east::Union{AbstractFloat,Integer}=360.0,
                    west::Union{AbstractFloat,Integer}=0.0, dealloc::Bool=false,
                    exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert interval > 0
    grid = Array{Cdouble}(undef, floor(Int, 180 / interval + 1),
                          floor(Int, 360 / interval + 1))
    MakeGrid2d!(grid, cilm, lmax, interval; norm=norm, csphase=csphase, f=f,
                a=a, north=north, south=south, east=east, west=west,
                dealloc=dealloc, exitstatus=exitstatus)
    return grid
end

export MakeGridPoint
function MakeGridPoint(cilm::AbstractArray{Cdouble,3}, lmax::Integer,
                       lat::Union{AbstractFloat,Integer},
                       lon::Union{AbstractFloat,Integer}; norm::Integer=1,
                       csphase::Integer=1, dealloc::Bool=false)
    @assert lmax ≥ 0
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax + 1
    @assert size(cilm, 3) ≥ size(cilm, 2)
    @assert norm ∈ (1, 2, 3, 4)
    @assert csphase ∈ (1, -1)
    value = ccall((:MakeGridPoint, libSHTOOLS), Cdouble,
                  (Ptr{Cdouble}, Cint, Cint, Cdouble, Cdouble, Ref{Cint},
                   Ref{Cint}, Ref{Cint}), cilm, size(cilm, 2), lmax, lat, lon,
                  norm, csphase, dealloc)
    return Float64(value)
end

function MakeGridPoint(cilm::AbstractArray{Complex{Cdouble},3}, lmax::Integer,
                       lat::Union{AbstractFloat,Integer},
                       lon::Union{AbstractFloat,Integer}; norm::Integer=1,
                       csphase::Integer=1, dealloc::Bool=false)
    @assert lmax ≥ 0
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax + 1
    @assert size(cilm, 3) ≥ size(cilm, 2)
    @assert norm ∈ (1, 2, 3, 4)
    @assert csphase ∈ (1, -1)
    value = ccall((:MakeGridPointC, libSHTOOLS), Complex{Cdouble},
                  (Ptr{Complex{Cdouble}}, Cint, Cint, Cdouble, Cdouble,
                   Ref{Cint}, Ref{Cint}, Ref{Cint}), cilm, size(cilm, 2), lmax,
                  lat, lon, norm, csphase, dealloc)
    return Complex{Float64}(value)
end

export MakeGridPointC
const MakeGridPointC = MakeGridPoint

export SHMultiply!
function SHMultiply!(cilmout::AbstractArray{Cdouble,3},
                     cilm1::AbstractArray{Cdouble,3}, lmax1::Integer,
                     cilm2::AbstractArray{Cdouble,3}, lmax2::Integer;
                     precomp::Bool=false, norm::Integer=1, csphase::Integer=1,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax1 ≥ 0
    @assert lmax2 ≥ 0
    @assert size(cilmout, 1) == 2
    @assert size(cilmout, 2) ≥ lmax1 + lmax2 + 1
    @assert size(cilmout, 3) == size(cilmout, 2)
    @assert size(cilm1, 1) == 2
    @assert size(cilm1, 2) ≥ lmax1 + 1
    @assert size(cilm1, 3) == size(cilm1, 2)
    @assert size(cilm2, 1) == 2
    @assert size(cilm2, 2) ≥ lmax2 + 1
    @assert size(cilm2, 3) == size(cilm2, 2)
    @assert norm ∈ (1, 2, 3, 4)
    @assert csphase ∈ (1, -1)
    exitstatus′ = Ref{Cint}()
    ccall((:SHMultiply, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Cint,
           Cint, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}), cilmout,
          size(cilmout, 2), cilm1, size(cilm1, 2), lmax1, cilm2, size(cilm2, 2),
          lmax2, precomp, norm, csphase, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 && error("SHMultiply!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return cilmout
end

export SHMultiply
function SHMultiply(cilm1::AbstractArray{Cdouble,3}, lmax1::Integer,
                    cilm2::AbstractArray{Cdouble,3}, lmax2::Integer;
                    precomp::Bool=false, norm::Integer=1, csphase::Integer=1,
                    exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax1 ≥ 0
    @assert lmax2 ≥ 0
    cilmout = Array{Cdouble}(undef, 2, lmax1 + lmax2 + 1, lmax1 + lmax2 + 1)
    SHMultiply!(cilmout, cilm1, lmax1, cilm2, lmax2; precomp=precomp, norm=norm,
                csphase=csphase, exitstatus=exitstatus)
    return cilmout
end

end
