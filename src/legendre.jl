# Legendre functions

# 4π normalized
# Orthonormalized
# Schmidt semi-normalized
# Unnormalized
# Utilities

"""
    PlmBar!(p::AbstractVector{Cdouble},
            lmax::Integer,
            z::Union{AbstractFloat,Integer};
            csphase::Integer=1,
            cnorm::Integer=0,
            exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the 4π-normalized associated Legendre functions.

See also: [`PlmBar`](@ref), [`PlmBar_d1!`](@ref), [`PlBar!`](@ref),
[`PlmON!`](@ref), [`PlmSchmidt!`](@ref), [`PLegendreA!`](@ref),
[`PlmIndex`](@ref)
"""
function PlmBar! end

"""
    p = PlmBar(lmax::Integer,
               z::Union{AbstractFloat,Integer};
               csphase::Integer=1,
               cnorm::Integer=0,
               exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}

Compute all the 4π-normalized associated Legendre functions.

See also: [`PlmBar!`](@ref), [`PlmBar_d1`](@ref), [`PlBar`](@ref),
[`PlmON`](@ref), [`PlmSchmidt`](@ref), [`PLegendreA`](@ref),
[`PlmIndex`](@ref)
"""
function PlmBar end

"""
    PlmBar_d1!(p::AbstractVector{Cdouble},
               dp::AbstractVector{Cdouble}, 
               lmax::Integer,
               z::Union{AbstractFloat,Integer};
               csphase::Integer=1,
               cnorm::Integer=0,
               exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the 4π-normalized associated Legendre functions and first
derivatives.

See also: [`PlmBar_d1`](@ref), [`PlmBar_d1`](@ref),
[`PlBar_d1!`](@ref), [`PlmON_d1!`](@ref), [`PlmSchmidt_d1!`](@ref),
[`PLegendreA_d1!`](@ref), [`PlmIndex`](@ref)
"""
function PlmBar_d1! end

"""
    p, dp = PlmBar_d1(lmax::Integer,
                      z::Union{AbstractFloat,Integer};
                      csphase::Integer=1,
                      cnorm::Integer=0,
                      exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}
    dp::Vector{Cdouble}

Compute all the 4π-normalized associated Legendre functions and first
derivatives.

See also: [`PlmBar`](@ref), [`PlmBar_d1!`](@ref), [`PlBar_d1`](@ref),
[`PlmON_d1`](@ref), [`PlmSchmidt_d1`](@ref),
[`PLegendreA_d1`](@ref), [`PlmIndex`](@ref)
"""
function PlmBar_d1 end

"""
    PlmON!(p::AbstractVector{Cdouble},
           lmax::Integer,
           z::Union{AbstractFloat,Integer};
           csphase::Integer=1,
           cnorm::Integer=0,
           exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the orthonormalized associated Legendre functions.

See also: [`PlmON`](@ref), [`PlmON_d1!`](@ref), [`PlON!`](@ref),
[`PlmBar!`](@ref), [`PlmSchmidt!`](@ref), [`PLegendreA!`](@ref),
[`PlmIndex`](@ref)
"""
function PlmON! end

"""
    p = PlmON(lmax::Integer,
              z::Union{AbstractFloat,Integer};
              csphase::Integer=1,
              cnorm::Integer=0,
              exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}

Compute all the orthonormalized associated Legendre functions.

See also: [`PlmON!`](@ref), [`PlmON_d1`](@ref), [`PlON`](@ref),
[`PlmBar`](@ref), [`PlmSchmidt`](@ref), [`PLegendreA`](@ref),
[`PlmIndex`](@ref)
"""
function PlmON end

"""
    PlmON_d1!(p::AbstractVector{Cdouble},
              dp::AbstractVector{Cdouble}, 
              lmax::Integer,
              z::Union{AbstractFloat,Integer};
              csphase::Integer=1,
              cnorm::Integer=0,
              exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the orthonormalized associated Legendre functions and
first derivatives.

See also: [`PlmON_d1`](@ref), [`PlmON_d1`](@ref), [`PlON_d1!`](@ref),
[`PlmBar_d1!`](@ref), [`PlmSchmidt_d1!`](@ref),
[`PLegendreA_d1!`](@ref), [`PlmIndex`](@ref)
"""
function PlmON_d1! end

"""
    p, dp = PlmON_d1(lmax::Integer,
                     z::Union{AbstractFloat,Integer};
                     csphase::Integer=1,
                     cnorm::Integer=0,
                     exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}
    dp::Vector{Cdouble}

Compute all the orthonormalized associated Legendre functions and first
derivatives.

See also: [`PlmON`](@ref), [`PlmON_d1!`](@ref), [`PlON_d1`](@ref),
[`PlmBar_d1`](@ref), [`PlmSchmidt_d1`](@ref), [`PLegendreA_d1`](@ref),
[`PlmIndex`](@ref)
"""
function PlmON_d1 end

"""
    PlmSchmidt!(p::AbstractVector{Cdouble},
                lmax::Integer,
                z::Union{AbstractFloat,Integer};
                csphase::Integer=1,
                cnorm::Integer=0,
                exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the Schmidt-normalized associated Legendre functions.

See also: [`PlmSchmidt`](@ref), [`PlmSchmidt_d1!`](@ref),
[`PlSchmidt!`](@ref), [`PlmBar!`](@ref), [`PlmON!`](@ref),
[`PLegendreA!`](@ref), [`PlmIndex`](@ref)
"""
function PlmSchmidt! end

"""
    p = PlmSchmidt(lmax::Integer,
                   z::Union{AbstractFloat,Integer};
                   csphase::Integer=1,
                   cnorm::Integer=0,
                   exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}

Compute all the Schmidt-normalized associated Legendre functions.

See also: [`PlmSchmidt!`](@ref), [`PlmSchmidt_d1`](@ref),
[`PlSchmidt`](@ref), [`PlmBar`](@ref), [`PlmON`](@ref),
[`PLegendreA`](@ref), [`PlmIndex`](@ref)
"""
function PlmSchmidt end

"""
    PlmSchmidt_d1!(p::AbstractVector{Cdouble},
                   dp::AbstractVector{Cdouble}, 
                   lmax::Integer,
                   z::Union{AbstractFloat,Integer};
                   csphase::Integer=1,
                   cnorm::Integer=0,
                   exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the Schmidt-normalized associated Legendre functions and
first derivatives.

See also: [`PlmSchmidt_d1`](@ref), [`PlmSchmidt_d1`](@ref),
[`PlSchmidt_d1!`](@ref), [`PlmBar_d1!`](@ref), [`PlmON_d1!`](@ref),
[`PLegendreA_d1!`](@ref), [`PlmIndex`](@ref)
"""
function PlmSchmidt_d1! end

"""
    p, dp = PlmSchmidt_d1(lmax::Integer,
                          z::Union{AbstractFloat,Integer};
                          csphase::Integer=1,
                          cnorm::Integer=0,
                          exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}
    dp::Vector{Cdouble}

Compute all the Schmidt-normalized associated Legendre functions and
first derivatives.

See also: [`PlmSchmidt`](@ref), [`PlmSchmidt_d1!`](@ref),
[`PlSchmidt_d1`](@ref), [`PlmBar_d1`](@ref), [`PlmON_d1`](@ref),
[`PLegendreA_d1`](@ref), [`PlmIndex`](@ref)
"""
function PlmSchmidt_d1 end

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

"""
    PLegendreA!(p::AbstractVector{Cdouble},
                lmax::Integer,
                z::Union{AbstractFloat,Integer};
                csphase::Integer=1,
                exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the unnormalized associated Legendre functions.

See also: [`PLegendreA`](@ref), [`PLegendreA_d1!`](@ref),
[`PLegendre!`](@ref), [`PlmBar!`](@ref), [`PlmON!`](@ref),
[`PlmSchmidt!`](@ref), [`PlmIndex`](@ref)
"""
function PLegendreA! end

"""
    p = PLegendreA(lmax::Integer,
                   z::Union{AbstractFloat,Integer};
                   csphase::Integer=1,
                   exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}

Compute all the unnormalized associated Legendre functions.

See also: [`PLegendreA!`](@ref), [`PLegendreA_d1`](@ref),
[`PLegendre`](@ref), [`PlmBar`](@ref), [`PlmON`](@ref),
[`PlmSchmidt`](@ref), [`PlmIndex`](@ref)
"""
function PLegendreA end

"""
    PLegendreA_d1!(p::AbstractVector{Cdouble},
                   dp::AbstractVector{Cdouble}, 
                   lmax::Integer,
                   z::Union{AbstractFloat,Integer};
                   csphase::Integer=1,
                   exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the unnormalized associated Legendre functions and first
derivatives.

See also: [`PLegendreA_d1`](@ref), [`PLegendreA_d1`](@ref),
[`PLegendre_d1!`](@ref), [`PlmBar_d1!`](@ref), [`PlmON_d1!`](@ref),
[`PlmSchmidt_d1!`](@ref), [`PlmIndex`](@ref)
"""
function PLegendreA_d1! end

"""
    p, dp = PLegendreA_d1(lmax::Integer,
                            z::Union{AbstractFloat,Integer};
                            csphase::Integer=1,
                            exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}
    dp::Vector{Cdouble}

Compute all the unnormalized associated Legendre functions and first
derivatives.

See also: [`PLegendreA`](@ref), [`PLegendreA_d1!`](@ref),
[`PLegendre_d1`](@ref), [`PlmBar_d1`](@ref), [`PlmON_d1`](@ref),
[`PlmSchmidt_d1`](@ref), [`PlmIndex`](@ref)
"""
function PLegendreA_d1 end

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

"""
    PlBar!(p::AbstractVector{Cdouble},
           lmax::Integer,
           z::Union{AbstractFloat,Integer};
           exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the 4π-normalized Legendre polynomials.

See also: [`PlBar`](@ref), [`PlBar_d1!`](@ref), [`PlmBar!`](@ref),
[`PlON!`](@ref), [`PlSchmidt!`](@ref), [`PLegendre!`](@ref)
"""
function PlBar! end

"""
    p = PlBar(lmax::Integer,
               z::Union{AbstractFloat,Integer};
               exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}

Compute all the 4π-normalized associated Legendre functions.

See also: [`PlBar!`](@ref), [`PlBar_d1`](@ref), [`PlmBar`](@ref),
[`PlON`](@ref), [`PlSchmidt`](@ref), [`PLegendre`](@ref)
"""
function PlBar end

"""
    PlBar_d1!(p::AbstractVector{Cdouble},
               dp::AbstractVector{Cdouble}, 
               lmax::Integer,
               z::Union{AbstractFloat,Integer};
               exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the 4π-normalized associated Legendre functions and first
derivatives.

See also: [`PlBar_d1`](@ref), [`PlBar_d1`](@ref),
[`PlmBar_d1!`](@ref), [`PlON_d1!`](@ref), [`PlSchmidt_d1!`](@ref),
[`PLegendre_d1!`](@ref)
"""
function PlBar_d1! end

"""
    p, dp = PlBar_d1(lmax::Integer,
                      z::Union{AbstractFloat,Integer};
                      exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}
    dp::Vector{Cdouble}

Compute all the 4π-normalized associated Legendre functions and first
derivatives.

See also: [`PlBar`](@ref), [`PlBar_d1!`](@ref), [`PlmBar_d1`](@ref),
[`PlON_d1`](@ref), [`PlSchmidt_d1`](@ref), [`PLegendre_d1`](@ref)
"""
function PlBar_d1 end

"""
    PlON!(p::AbstractVector{Cdouble},
           lmax::Integer,
           z::Union{AbstractFloat,Integer};
           exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the orthonormalized associated Legendre functions.

See also: [`PlON`](@ref), [`PlON_d1!`](@ref), [`PlmON!`](@ref),
[`PlBar!`](@ref), [`PlSchmidt!`](@ref), [`PLegendre!`](@ref)
"""
function PlON! end

"""
    p = PlON(lmax::Integer,
              z::Union{AbstractFloat,Integer};
              exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}

Compute all the orthonormalized associated Legendre functions.

See also: [`PlON!`](@ref), [`PlON_d1`](@ref), [`PlmON`](@ref),
[`PlBar`](@ref), [`PlSchmidt`](@ref), [`PLegendre`](@ref)
"""
function PlON end

"""
    PlON_d1!(p::AbstractVector{Cdouble},
              dp::AbstractVector{Cdouble}, 
              lmax::Integer,
              z::Union{AbstractFloat,Integer};
              exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the orthonormalized associated Legendre functions and
first derivatives.

See also: [`PlON_d1`](@ref), [`PlON_d1`](@ref), [`PlmON_d1!`](@ref),
[`PlBar_d1!`](@ref), [`PlSchmidt_d1!`](@ref), [`PLegendre_d1!`](@ref)
"""
function PlON_d1! end

"""
    p, dp = PlON_d1(lmax::Integer,
                     z::Union{AbstractFloat,Integer};
                     exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}
    dp::Vector{Cdouble}

Compute all the orthonormalized associated Legendre functions and first
derivatives.

See also: [`PlON`](@ref), [`PlON_d1!`](@ref), [`PlmON_d1`](@ref),
[`PlBar_d1`](@ref), [`PlSchmidt_d1`](@ref), [`PLegendre_d1`](@ref)
"""
function PlON_d1 end

"""
    PlSchmidt!(p::AbstractVector{Cdouble},
                lmax::Integer,
                z::Union{AbstractFloat,Integer};
                exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the Schmidt-normalized associated Legendre functions.

See also: [`PlSchmidt`](@ref), [`PlSchmidt_d1!`](@ref),
[`PlmSchmidt!`](@ref), [`PlBar!`](@ref), [`PlON!`](@ref),
[`PLegendre!`](@ref)
"""
function PlSchmidt! end

"""
    p = PlSchmidt(lmax::Integer,
                   z::Union{AbstractFloat,Integer};
                   exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}

Compute all the Schmidt-normalized associated Legendre functions.

See also: [`PlSchmidt!`](@ref), [`PlSchmidt_d1`](@ref),
[`PlmSchmidt`](@ref), [`PlBar`](@ref), [`PlON`](@ref),
[`PLegendre`](@ref)
"""
function PlSchmidt end

"""
    PlSchmidt_d1!(p::AbstractVector{Cdouble},
                   dp::AbstractVector{Cdouble}, 
                   lmax::Integer,
                   z::Union{AbstractFloat,Integer};
                   exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the Schmidt-normalized associated Legendre functions and
first derivatives.

See also: [`PlSchmidt_d1`](@ref), [`PlSchmidt_d1`](@ref),
[`PlmSchmidt_d1!`](@ref), [`PlBar_d1!`](@ref), [`PlON_d1!`](@ref),
[`PLegendre_d1!`](@ref)
"""
function PlSchmidt_d1! end

"""
    p, dp = PlSchmidt_d1(lmax::Integer,
                          z::Union{AbstractFloat,Integer};
                          exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}
    dp::Vector{Cdouble}

Compute all the Schmidt-normalized associated Legendre functions and
first derivatives.

See also: [`PlSchmidt`](@ref), [`PlSchmidt_d1!`](@ref),
[`PlmSchmidt_d1`](@ref), [`PlBar_d1`](@ref), [`PlON_d1`](@ref),
[`PLegendre_d1`](@ref)
"""
function PlSchmidt_d1 end

"""
    PLegendre!(p::AbstractVector{Cdouble},
                lmax::Integer,
                z::Union{AbstractFloat,Integer};
                exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the unnormalized associated Legendre functions.

See also: [`PLegendre`](@ref), [`PLegendre_d1!`](@ref),
[`PLegendreA!`](@ref), [`PlBar!`](@ref), [`PlON!`](@ref),
[`PlSchmidt!`](@ref)
"""
function PLegendre! end

"""
    p = PLegendre(lmax::Integer,
                   z::Union{AbstractFloat,Integer};
                   exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}

Compute all the unnormalized associated Legendre functions.

See also: [`PLegendre!`](@ref), [`PLegendre_d1`](@ref),
[`PLegendreA`](@ref), [`PlBar`](@ref), [`PlON`](@ref),
[`PlSchmidt`](@ref)
"""
function PLegendre end

"""
    PLegendre_d1!(p::AbstractVector{Cdouble},
                   dp::AbstractVector{Cdouble}, 
                   lmax::Integer,
                   z::Union{AbstractFloat,Integer};
                   exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Compute all the unnormalized associated Legendre functions and first
derivatives.

See also: [`PLegendre_d1`](@ref), [`PLegendre_d1`](@ref),
[`PLegendreA_d1!`](@ref), [`PlBar_d1!`](@ref), [`PlON_d1!`](@ref),
[`PlSchmidt_d1!`](@ref)
"""
function PLegendre_d1! end

"""
    p, dp = PLegendre_d1(lmax::Integer,
                            z::Union{AbstractFloat,Integer};
                            exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    p::Vector{Cdouble}
    dp::Vector{Cdouble}

Compute all the unnormalized associated Legendre functions and first
derivatives.

See also: [`PLegendre`](@ref), [`PLegendre_d1!`](@ref),
[`PLegendre_d1!`](@ref), [`PlBar_d1`](@ref), [`PlON_d1`](@ref),
[`PlSchmidt_d1`](@ref)
"""
function PLegendre_d1 end

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
"""
    index = PlmIndex(l::Integer,
                     m::Integer)
    index::Int

Compute the index of an array of Legendre functions corresponding to
degree `l` and angular order `m`.
"""
function PlmIndex(l::Integer, m::Integer)
    0 ≤ m ≤ l ||
        error("m must be greater than or equal to zero and less than or equal to l: m=$m, l=$l")
    index = ccall((:PlmIndex, libSHTOOLS), Cint, (Cint, Cint), l, m)
    return Int(index)
end
