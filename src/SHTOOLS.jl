module SHTOOLS

using SHTOOLS_jll

const Maybe{T} = Union{Nothing,T}
cmaybe(::Type{T}, x::Maybe) where {T} = x === nothing ? Ptr{T}() : Ref{T}(x)
function cmayberef(::Type{T}, x::Maybe{Ref{T}}) where {T}
    return x === nothing ? Ptr{T}() : (x::Ref{T})
end

################################################################################

for op in [:PlmBar, :PlmON, :PlmSchmidt, :PLegendreA]
    op_d1 = Symbol(op, "_d1")
    @eval begin
        #
        export $op
        function $op(lmax::Union{Cint,Int}, z::Union{Cdouble,Float64,Int};
                     csphase::Maybe{Int}=nothing, cnorm::Maybe{Int}=nothing,
                     exitstatus::Maybe{Ref{Int}}=nothing)
            p = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            ccall(($(QuoteNode(op)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Cint, Cdouble, Ref{Cint}, Ref{Cint},
                   Ref{Cint}), p, lmax, z, cmaybe(Cint, csphase),
                  cmaybe(Cint, cnorm), cmayberef(Cint, exitstatus))
            return p
        end
        #
        export $op_d1
        function $op_d1(lmax::Union{Cint,Int}, z::Union{Cdouble,Float64,Int};
                        csphase::Maybe{Int}=nothing, cnorm::Maybe{Int}=nothing,
                        exitstatus::Maybe{Ref{Int}}=nothing)
            p = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            dp = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            ccall(($(QuoteNode(op_d1)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Ref{Cint},
                   Ref{Cint}, Ref{Cint}), p, dp, lmax, z, cmaybe(Cint, csphase),
                  cmaybe(Cint, cnorm), cmayberef(Cint, exitstatus))
            return p, dp
        end
        #
    end
end

for op in [:PlBar, :PlON, :PlSchmidt, :PLegendre]
    op_d1 = Symbol(op, "_d1")
    @eval begin
        #
        export $op
        function $op(lmax::Union{Cint,Int}, z::Union{Cdouble,Float64,Int};
                     exitstatus::Maybe{Ref{Int}}=nothing)
            p = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            ccall(($(QuoteNode(op)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Cint, Cdouble, Ref{Cint}), p, lmax, z,
                  cmayberef(Cint, exitstatus))
            return p
        end
        #
        export $op_d1
        function $op_d1(lmax::Union{Cint,Int}, z::Union{Cdouble,Float64,Int};
                        exitstatus::Maybe{Ref{Int}}=nothing)
            p = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            dp = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            ccall(($(QuoteNode(op_d1)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Ref{Cint}), p, dp,
                  lmax, z, cmayberef(Cint, exitstatus))
            return p, dp
        end
        #
    end
end

export PlmIndex
function PlmIndex(l::Union{Cint,Int}, m::Union{Cint,Int})::Int
    0 ≤ m ≤ l ||
        error("m must be greater than or equal to zero and less than or equal to l: m=$m, l=$l")
    index = ccall((:PlmIndex, libSHTOOLS), Cint, (Cint, Cint), l, m)
    return index
end

end
