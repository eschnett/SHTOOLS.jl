module SHTOOLS

using SHTOOLS_jll

const Optional{T} = Union{Nothing,T}
# For optional arguments
optional(x::Optional, x0) = x !== nothing ? x : x0

include("legendre.jl")
include("transforms.jl")
include("storage.jl")

end
