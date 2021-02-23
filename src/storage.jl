# Spherical harmonic I/O, storage, and conversions

# Spherical harmonic I/O

# SHRead
# SHRead2
# SHReadJPL

# Spherical harmonic storage

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
