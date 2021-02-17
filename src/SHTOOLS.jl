module SHTOOLS

using SHTOOLS_jll

const Optional{T} = Union{Nothing,T}
# For optional arguments
optional(x::Optional, x0) = x === nothing ? x0 : x
function optional(::Type{T}, x::Optional{T}) where {T}
    return x === nothing ? Ptr{T}() : Ref{T}(x)
end

################################################################################

# Legendre functions

for op in [:PlmBar, :PlmON, :PlmSchmidt]
    op! = Symbol(op, "!")
    op_d1 = Symbol(op, "_d1")
    op_d1! = Symbol(op, "_d1!")
    @eval begin
        #
        export $op!
        function $op!(p::AbstractVector{Cdouble}, lmax::Integer,
                      z::Union{AbstractFloat,Integer}; csphase::Integer=1,
                      cnorm::Integer=1,
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
                     csphase::Integer=1, cnorm::Integer=1,
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
                         z::Union{AbstractVector,Integer}; csphase::Integer=1,
                         cnorm::Integer=1,
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
                        csphase::Integer=1, cnorm::Integer=1,
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
                         z::Union{AbstractVector,Integer}; csphase::Integer=1,
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

export SHExpandDH!
function SHExpandDH!(cilm::Array{Cdouble,3}, lmax::Optional{Ref{<:Integer}},
                     griddh::Array{Cdouble,2}, n::Integer; norm::Integer=1,
                     sampling::Integer=1, csphase::Integer=1,
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
    if lmax !== nothing
        lmax[] = lmax′[]
    end
    return cilm, Int(lmax′[])
end

export SHExpandDH
function SHExpandDH(griddh::Array{Cdouble,2}, n::Integer; norm::Integer=1,
                    sampling::Integer=1, csphase::Integer=1,
                    lmax_calc::Optional{Integer}=nothing,
                    exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert n > 0 && iseven(n)
    lmax_calc′ = optional(lmax_calc, n ÷ 2 - 1)
    cilm = Array{Cdouble}(undef, 2, lmax_calc′ + 1, lmax_calc′ + 1)
    lmax = Ref{Int}()
    SHExpandDH!(cilm, lmax, griddh, n; norm=norm, sampling=sampling,
                csphase=csphase, lmax_calc=lmax_calc, exitstatus=exitstatus)
    return cilm, lmax[]
end

export MakeGridDH!
function MakeGridDH!(griddh::Array{Cdouble,2}, n::Optional{Ref{<:Integer}},
                     cilm::Array{Cdouble,3}, lmax::Integer; norm::Integer=1,
                     sampling::Integer=1, csphase::Integer=1,
                     lmax_calc::Optional{Integer}=nothing, extend::Integer=0,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert norm ∈ (1, 2, 3, 4)
    @assert sampling ∈ (1, 2)
    @assert csphase ∈ (1, -1)
    @assert extend ∈ (0, 1)
    n′ = 2 * lmax + 2
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
    if n !== nothing
        n[] = n′
    end
    return griddh, n′
end

export MakeGridDH
function MakeGridDH(cilm::Array{Cdouble,3}, lmax::Integer; norm::Integer=1,
                    sampling::Integer=1, csphase::Integer=1,
                    lmax_calc::Optional{Integer}=nothing, extend::Integer=0,
                    exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert sampling ∈ (1, 2)
    @assert extend ∈ (0, 1)
    n′ = 2 * lmax + 2
    griddh = Array{Cdouble}(undef, n′ + extend, sampling * n′ + extend)
    _, n = MakeGridDH!(griddh, nothing, cilm, lmax; norm=norm,
                       sampling=sampling, csphase=csphase, lmax_calc=lmax_calc,
                       extend=extend, exitstatus=exitstatus)
    return griddh, n
end

export SHExpandDHC!
function SHExpandDHC!(cilm::Array{Complex{Cdouble},3},
                      lmax::Optional{Ref{<:Integer}},
                      griddh::Array{Complex{Cdouble},2}, n::Integer;
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
    if lmax !== nothing
        lmax[] = lmax′[]
    end
    return cilm, Int(lmax′[])
end

export SHExpandDHC
function SHExpandDHC(griddh::Array{Complex{Cdouble},2}, n::Integer;
                     norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                     lmax_calc::Optional{Integer}=nothing,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert n > 0 && iseven(n)
    lmax_calc′ = optional(lmax_calc, n ÷ 2 - 1)
    cilm = Array{Complex{Cdouble}}(undef, 2, lmax_calc′ + 1, lmax_calc′ + 1)
    lmax = Ref{Int}()
    SHExpandDHC!(cilm, lmax, griddh, n; norm=norm, sampling=sampling,
                 csphase=csphase, lmax_calc=lmax_calc, exitstatus=exitstatus)
    return cilm, lmax[]
end

export MakeGridDHC!
function MakeGridDHC!(griddh::Array{Complex{Cdouble},2},
                      n::Optional{Ref{<:Integer}},
                      cilm::Array{Complex{Cdouble},3}, lmax::Integer;
                      norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                      lmax_calc::Optional{Integer}=nothing, extend::Integer=0,
                      exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert norm ∈ (1, 2, 3, 4)
    @assert sampling ∈ (1, 2)
    @assert csphase ∈ (1, -1)
    @assert extend ∈ (0, 1)
    n′ = 2 * lmax + 2
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
    if n !== nothing
        n[] = n′
    end
    return griddh, n′
end

export MakeGridDHC
function MakeGridDHC(cilm::Array{Complex{Cdouble},3}, lmax::Integer;
                     norm::Integer=1, sampling::Integer=1, csphase::Integer=1,
                     lmax_calc::Optional{Integer}=nothing, extend::Integer=0,
                     exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert sampling ∈ (1, 2)
    @assert extend ∈ (0, 1)
    n′ = 2 * lmax + 2
    griddh = Array{Complex{Cdouble}}(undef, n′ + extend, sampling * n′ + extend)
    _, n = MakeGridDHC!(griddh, nothing, cilm, lmax; norm=norm,
                        sampling=sampling, csphase=csphase, lmax_calc=lmax_calc,
                        extend=extend, exitstatus=exitstatus)
    return griddh, n
end

# TODO: Add MakeGradientDH once there is a C wrapper

end
