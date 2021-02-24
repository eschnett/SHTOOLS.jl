# Spherical harmonic I/O, storage, and conversions

@testset "SHCilmToCindex and SHCindexToCilm" begin
    Random.seed!(0)
    # Choose random values
    lmax = 2
    cilm = randn(2, lmax + 1, lmax + 1)

    cindex = SHCilmToCindex(cilm)
    @test size(cindex) == (2, (lmax + 1) * (lmax + 2) ÷ 2)
    # These loops test all values in cindex
    for l in 0:lmax, m in 0:l, i in 1:(1 + (m > 0))
        @test cindex[i, l * (l + 1) ÷ 2 + m + 1] == cilm[i, l + 1, m + 1]
    end

    cilm′ = SHCindexToCilm(cindex)
    for l in 0:lmax, m in 0:l, i in 1:(1 + (m > 0))
        @test isapprox(cilm′[i, l + 1, m + 1], cilm[i, l + 1, m + 1])
    end

    cindex′ = SHCilmToCindex(cilm′)
    @test isapprox(cindex′, cindex)
end

@testset "SHCilmToVector SHVectorToCilm, and YlmIndexVector" begin
    Random.seed!(0)
    # Choose random values
    lmax = 2
    cilm = randn(2, lmax + 1, lmax + 1)

    vector = SHCilmToVector(cilm, lmax)
    @test length(vector) == (lmax + 1)^2
    coverage = falses(length(vector))
    for l in 0:lmax, m in 0:l, i in 1:(1 + (m > 0))
        index = YlmIndexVector(i, l, m)
        @test !coverage[index]
        coverage[index] = true
        @test vector[index] == cilm[i, l + 1, m + 1]
    end
    @test all(coverage)

    cilm′ = SHVectorToCilm(vector, lmax)
    for l in 0:lmax, m in 0:l, i in 1:(1 + (m > 0))
        @test isapprox(cilm′[i, l + 1, m + 1], cilm[i, l + 1, m + 1])
    end

    vector′ = SHCilmToVector(cilm′, lmax)
    @test isapprox(vector′, vector)
end
