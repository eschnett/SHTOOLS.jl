# Spherical harmonic I/O, storage, and conversions

@testset "SHCilmToVector and SHVectorToCilm" begin
    Random.seed!(0)
    # Choose random values
    lmax = 2
    cilm = randn(2, lmax + 1, lmax + 1)

    vector = SHCilmToVector(cilm, lmax)
    @test length(vector) == (lmax + 1)^2
    # These loops test all values in the vector
    for l in 0:lmax, m in 0:l, i in 1:(1 + (m > 0))
        @test vector[l^2 + (i - 1) * l + m + 1] == cilm[i, l + 1, m + 1]
    end

    cilm′ = SHVectorToCilm(vector, lmax)
    for l in 0:lmax, m in 0:l, i in 1:(1 + (m > 0))
        @test isapprox(cilm′[i, l + 1, m + 1], cilm[i, l + 1, m + 1])
    end

    vector′ = SHCilmToVector(cilm′, lmax)
    @test isapprox(vector′, vector)
end
