# Generate documentation with this command:
# (cd docs && julia --color=yes make.jl)

push!(LOAD_PATH, "..")

using Documenter
using SHTOOLS

makedocs(; sitename="SHTOOLS", format=Documenter.HTML(), modules=[SHTOOLS])

deploydocs(; repo="github.com/eschnett/JuliaFormatter.jl.git", devbranch="main",
           push_preview=true)
