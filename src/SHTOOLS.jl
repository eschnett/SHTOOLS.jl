module SHTOOLS

using SHTOOLS_jll

const Maybe{T} = Union{Nothing,T}
# For optional arguments
cmaybe(::Type{T}, x::Maybe{T}) where {T} = x === nothing ? Ptr{T}() : Ref{T}(x)

################################################################################

for op in [:PlmBar, :PlmON, :PlmSchmidt, :PLegendreA]
    op! = Symbol(op, "!")
    op_d1 = Symbol(op, "_d1")
    op_d1! = Symbol(op, "_d1!")
    @eval begin
        #
        export $op!
        function $op!(p::AbstractVector{Cdouble}, lmax::Integer,
                      z::Union{AbstractFloat,Integer};
                      csphase::Maybe{Integer}=nothing,
                      cnorm::Maybe{Integer}=nothing,
                      exitstatus::Maybe{Ref{<:Integer}}=nothing)
            @assert length(p) ≥ (lmax + 1) * (lmax + 2) ÷ 2
            exitstatus′ = 0
            ccall(($(QuoteNode(op)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Cint, Cdouble, Ref{Cint}, Ref{Cint},
                   Ref{Cint}), p, lmax, z, cmaybe(Cint, csphase),
                  cmaybe(Cint, cnorm), exitstatus′)
            if exitstatus === nothing
                exitstatus′ ≠ 0 && error("$($op!): Error code $exitstatus′")
            else
                exitstatus[] = exitstatus′
            end
            return p
        end
        #
        export $op
        function $op(lmax::Integer, z::Union{AbstractFloat,Integer};
                     csphase::Maybe{Integer}=nothing,
                     cnorm::Maybe{Integer}=nothing,
                     exitstatus::Maybe{Ref{<:Integer}}=nothing)
            p = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            $op!(p, lmax, z; csphase=csphase, cnorm=cnorm,
                 exitstatus=exitstatus)
            return p
        end
        #
        export $op_d1!
        function $op_d1!(p::AbstractVector{Cdouble},
                         dp::AbstractVector{Cdouble}, lmax::Integer,
                         z::Union{AbstractVector,Integer};
                         csphase::Maybe{Integer}=nothing,
                         cnorm::Maybe{Integer}=nothing,
                         exitstatus::Maybe{Ref{<:Integer}}=nothing)
            @assert length(p) ≥ lmax + 1
            @assert length(dp) ≥ lmax + 1
            exitstatus′ = 0
            ccall(($(QuoteNode(op_d1)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Ref{Cint},
                   Ref{Cint}, Ref{Cint}), p, dp, lmax, z, cmaybe(Cint, csphase),
                  cmaybe(Cint, cnorm), exitstatus′)
            if exitstatus === nothing
                exitstatus′ ≠ 0 && error("$($op_d1!): Error code $exitstatus′")
            else
                exitstatus[] = exitstatus′
            end
            return p, dp
        end
        #
        export $op_d1
        function $op_d1(lmax::Integer, z::Union{AbstractFloat,Integer};
                        csphase::Maybe{Integer}=nothing,
                        cnorm::Maybe{Integer}=nothing,
                        exitstatus::Maybe{Ref{<:Integer}}=nothing)
            p = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            dp = Array{Cdouble}(undef, (lmax + 1) * (lmax + 2) ÷ 2)
            $op_d1!(p, dp, lmax, z; csphase=csphase, cnorm=cnorm,
                    exitstatus=exitstatus)
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
                      exitstatus::Maybe{Ref{<:Integer}}=nothing)
            @assert length(p) ≥ lmax + 1
            exitstatus′ = 0
            ccall(($(QuoteNode(op)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Cint, Cdouble, Ref{Cint}), p, lmax, z,
                  exitstatus′)
            if exitstatus === nothing
                exitstatus′ ≠ 0 && error("$($op!): Error code $exitstatus′")
            else
                exitstatus[] = exitstatus′
            end
            return p
        end
        #
        export $op
        function $op(lmax::Integer, z::Union{AbstractFloat,Integer};
                     exitstatus::Maybe{Ref{<:Integer}}=nothing)
            p = Array{Cdouble}(undef, lmax + 1)
            $op!(p, lmax, z; exitstatus=exitstatus)
            return p
        end
        #
        export $op_d1!
        function $op_d1!(p::AbstractVector{Cdouble},
                         dp::AbstractVector{Cdouble}, lmax::Integer,
                         z::Union{AbstractFloat,Integer};
                         exitstatus::Maybe{Ref{<:Integer}}=nothing)
            @assert length(p) ≥ lmax + 1
            @assert length(dp) ≥ lmax + 1
            exitstatus′ = 0
            ccall(($(QuoteNode(op_d1)), libSHTOOLS), Cvoid,
                  (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Ref{Cint}), p, dp,
                  lmax, z, exitstatus′)
            if exitstatus === nothing
                exitstatus′ ≠ 0 && error("$($op_d1!): Error code $exitstatus′")
            else
                exitstatus[] = exitstatus′
            end
            return p, dp
        end
        #
        export $op_d1
        function $op_d1(lmax::Integer, z::Union{AbstractFloat,Integer};
                        exitstatus::Maybe{Ref{<:Integer}}=nothing)
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

end
