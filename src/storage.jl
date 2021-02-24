# Spherical harmonic I/O, storage, and conversions

# Spherical harmonic I/O

# SHRead
# SHRead2
# SHReadJPL

# Spherical harmonic storage

export SHCilmToCindex!
"""
    cindex = SHCilmToCindex!(cindex::AbstractArray{Cdouble,2},
                             cilm::AbstractArray{Cdouble,3};
                             degmax::Optional{Cint}=nothing,
                             exitstatus::Optional{Ref{<:Integer}}=nothing)
    cindex::AbstractArray{Cdouble,2}

Convert a three-dimensional array of spherical harmonic coefficients
to a two-dimensional indexed array.

See also: [`SHCilmToCindex`](@ref), [`SHCindexToCilm!`](@ref),
[`SHCilmToVector!`](@ref)
"""
function SHCilmToCindex!(cindex::AbstractArray{Cdouble,2},
                         cilm::AbstractArray{Cdouble,3};
                         degmax::Optional{Cint}=nothing,
                         exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert size(cilm, 1) == 2
    lmaxin = size(cilm, 2) - 1
    @assert lmaxin ≥ 0
    @assert size(cilm, 3) == size(cilm, 2)
    degmax′ = optional(degmax, lmaxin)
    @assert size(cindex, 1) == 2
    @assert size(cindex, 2) ≥ (degmax′ + 1) * (degmax′ + 2) ÷ 2
    exitstatus′ = Ref{Cint}()
    ccall((:SHCilmToCindex, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ref{Cint}, Ref{Cint}), cilm,
          size(cilm, 2), cindex, size(cindex, 2), degmax′, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 &&
            error("SHCilmToCindex!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return cindex
end

export SHCilmToCindex
"""
    cindex = SHCilmToCindex(cilm::AbstractArray{Cdouble,3};
                            degmax::Optional{Cint}=nothing,
                            exitstatus::Optional{Ref{<:Integer}}=nothing)
    cindex::Array{Cdouble,2}

Convert a three-dimensional array of spherical harmonic coefficients
to a two-dimensional indexed array.

See also: [`SHCilmToCindex!`](@ref), [`SHCindexToCilm`](@ref),
[`SHCilmToVector`](@ref)
"""
function SHCilmToCindex(cilm::AbstractArray{Cdouble,3};
                        degmax::Optional{Cint}=nothing,
                        exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert size(cilm, 1) == 2
    lmaxin = size(cilm, 2) - 1
    @assert lmaxin ≥ 0
    @assert size(cilm, 3) == size(cilm, 2)
    degmax′ = optional(degmax, lmaxin)
    cindex = Array{Cdouble}(undef, 2, (degmax′ + 1) * (degmax′ + 2) ÷ 2)
    SHCilmToCindex!(cindex, cilm; degmax=degmax, exitstatus=exitstatus)
    return cindex
end

export SHCindexToCilm!
"""
    cilm = SHCindexToCilm!(cilm::AbstractArray{Cdouble,3},
                           cindex::AbstractArray{Cdouble,2};
                           degmax::Optional{Cint}=nothing,
                           exitstatus::Optional{Ref{<:Integer}}=nothing)
    cilm::AbstractArray{Cdouble,3},

Convert a two-dimensional indexed array of spherical harmonic
coefficients to a three-dimensional array.

See also: [`SHCindexToCilm`](@ref), [`SHCilmToCindex!`](@ref),
[`SHVectorToCilm!`](@ref)
"""
function SHCindexToCilm!(cilm::AbstractArray{Cdouble,3},
                         cindex::AbstractArray{Cdouble,2};
                         degmax::Optional{Cint}=nothing,
                         exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert size(cindex, 1) == 2
    # (lmaxin + 1) * (lmaxin + 2) ÷ 2 = size(cindex, 2)
    lmaxin = round(Int, (sqrt(8 * size(cindex, 2) + 1) - 3) / 2)
    @assert size(cindex, 2) == (lmaxin + 1) * (lmaxin + 2) ÷ 2
    degmax′ = optional(degmax, lmaxin)
    @assert size(cindex, 2) ≥ degmax′
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ degmax′ + 1
    @assert size(cilm, 3) == size(cilm, 2)
    exitstatus′ = Ref{Cint}()
    ccall((:SHCindexToCilm, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ref{Cint}, Ref{Cint}),
          cindex, size(cindex, 2), cilm, size(cilm, 2), degmax′, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 &&
            error("SHCindexToCilm!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return cilm
end

export SHCindexToCilm
"""
    cilm = SHCindexToCilm(cindex::AbstractArray{Cdouble,2};
                          degmax::Optional{Cint}=nothing,
                          exitstatus::Optional{Ref{<:Integer}}=nothing)
    cilm::AbstractArray{Cdouble,3},

Convert a two-dimensional indexed array of spherical harmonic
coefficients to a three-dimensional array.

See also: [`SHCindexToCilm!`](@ref), [`SHCilmToCindex`](@ref),
[`SHVectorToCilm`](@ref)
"""
function SHCindexToCilm(cindex::AbstractArray{Cdouble,2},
                        degmax::Optional{Cint}=nothing,
                        exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert size(cindex, 1) == 2
    # (lmaxin + 1) * (lmaxin + 2) ÷ 2 = size(cindex, 2)
    lmaxin = round(Int, (sqrt(8 * size(cindex, 2) + 1) - 3) / 2)
    @assert size(cindex, 2) == (lmaxin + 1) * (lmaxin + 2) ÷ 2
    degmax′ = optional(degmax, lmaxin)
    @assert size(cindex, 2) ≥ degmax′
    cilm = Array{Cdouble}(undef, 2, degmax′ + 1, degmax′ + 1)
    SHCindexToCilm!(cilm, cindex; degmax=degmax, exitstatus=exitstatus)
    return cilm
end

export SHCilmToVector!
"""
    SHCilmToVector!(vector::AbstractVector{Cdouble},
                    cilm::AbstractArray{Cdouble,3},
                    lmax::Integer;
                    exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Convert a three-dimensional array of real spherical harmonic
coefficients to a 1-dimensional indexed vector.

See also: [`SHCilmToVector`](@ref)
"""
function SHCilmToVector!(vector::AbstractVector{Cdouble},
                         cilm::AbstractArray{Cdouble,3}, lmax::Integer;
                         exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert length(vector) ≥ (lmax + 1)^2
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax + 1
    @assert size(cilm, 3) == size(cilm, 2)
    exitstatus′ = Ref{Cint}()
    ccall((:SHCilmToVector, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ref{Cint}), cilm,
          size(cilm, 2), vector, lmax, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 &&
            error("SHCilmToVector!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return vector
end

export SHCilmToVector
"""
    vector = SHCilmToVector(cilm::AbstractArray{Cdouble,3},
                            lmax::Integer;
                            exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    vector::Vector{Cdouble}

Convert a three-dimensional array of real spherical harmonic
coefficients to a 1-dimensional indexed array.

See also: [`SHCilmToVector!`](@ref)
"""
function SHCilmToVector(cilm::AbstractArray{Cdouble,3}, lmax::Integer;
                        exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    vector = Array{Cdouble}(undef, (lmax + 1)^2)
    SHCilmToVector!(vector, cilm, lmax; exitstatus=exitstatus)
    return vector
end

export SHVectorToCilm!
"""
    SHVectorToCilm!(cilm::AbstractArray{Cdouble,3},
                    vector::AbstractVector{Cdouble},
                    lmax::Integer;
                    exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)

Convert a 1-dimensional indexed vector of real spherical harmonic
coefficients to a 3-dimensional array.

See also: [`SHVectorToCilm`](@ref)
"""
function SHVectorToCilm!(cilm::AbstractArray{Cdouble,3},
                         vector::AbstractVector{Cdouble}, lmax::Integer;
                         exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    @assert size(cilm, 1) == 2
    @assert size(cilm, 2) ≥ lmax + 1
    @assert size(cilm, 3) == size(cilm, 2)
    @assert length(vector) ≥ (lmax + 1)^2
    exitstatus′ = Ref{Cint}()
    ccall((:SHVectorToCilm, libSHTOOLS), Cvoid,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ref{Cint}), vector, cilm,
          size(cilm, 2), lmax, exitstatus′)
    if exitstatus === nothing
        exitstatus′[] ≠ 0 &&
            error("SHVectorToCilm!: Error code $(exitstatus′[])")
    else
        exitstatus[] = exitstatus′[]
    end
    return cilm
end

export SHVectorToCilm
"""
    cilm = SHVectorToCilm(cilm::AbstractArray{Cdouble,3},
                          vector::AbstractVector{Cdouble},
                          lmax::Integer;
                          exitstatus::Union{Nothing,Ref{<:Integer}}=nothing)
    cilm::Array{Cdouble,3}

Convert a 1-dimensional indexed vector of real spherical harmonic
coefficients to a 3-dimensional array.

See also: [`SHVectorToCilm!`](@ref)
"""
function SHVectorToCilm(vector::AbstractVector{Cdouble}, lmax::Integer;
                        exitstatus::Optional{Ref{<:Integer}}=nothing)
    @assert lmax ≥ 0
    cilm = Array{Cdouble}(undef, 2, lmax + 1, lmax + 1)
    @assert length(vector) ≥ (lmax + 1)^2
    exitstatus′ = Ref{Cint}()
    SHVectorToCilm!(cilm, vector, lmax; exitstatus=exitstatus)
    return cilm
end

export YlmIndexVector
"""
    index = YlmIndexVector(i, l, m)

Determine the index of a 1D ordered vector of spherical harmonic
coefficients corresponding to `i`, `l`, and `m`.

See also: [`SHCilmToVector!`](@ref), [`SHCilmToVector`](@ref),
[`SHVectorToCilm!`](@ref), [`SHVectorToCilm`](@ref)
"""
function YlmIndexVector(i::Integer, l::Integer, m::Integer)
    @assert 1 ≤ i ≤ 2
    @assert 0 ≤ l
    @assert 0 ≤ m ≤ l
    return l^2 + (i - 1) * l + m + 1
end
