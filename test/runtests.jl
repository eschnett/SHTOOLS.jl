using Random
using SHTOOLS
using Test

################################################################################

# Legendre Polynomials

@testset "PlmBar" begin
    p = PlmBar(4, 0)
    p′ = zeros(15)
    PlmBar!(p′, 4, 0)
    @test p′ == p

    p, dp = PlmBar_d1(4, 0)
    dp′ = zeros(15)
    PlmBar_d1!(p′, dp′, 4, 0)
    @test p′ == p
    @test dp′ == dp

    p = PlBar(4, 0)
    p′ = zeros(5)
    PlBar!(p′, 4, 0)
    @test p′ == p

    p, dp = PlBar_d1(4, 0)
    dp′ = zeros(5)
    PlBar_d1!(p′, dp′, 4, 0)
    @test p′ == p
    @test dp′ == dp
end

@testset "PlmON" begin
    p = PlmON(4, 0)
    p′ = zeros(15)
    PlmON!(p′, 4, 0)
    @test p′ == p

    p, dp = PlmON_d1(4, 0)
    dp′ = zeros(15)
    PlmON_d1!(p′, dp′, 4, 0)
    @test p′ == p
    @test dp′ == dp

    p = PlON(4, 0)
    p′ = zeros(5)
    PlON!(p′, 4, 0)
    @test p′ == p

    p, dp = PlON_d1(4, 0)
    dp′ = zeros(5)
    PlON_d1!(p′, dp′, 4, 0)
    @test p′ == p
    @test dp′ == dp
end

@testset "PlmSchmidt" begin
    p = PlmSchmidt(4, 0)
    p′ = zeros(15)
    PlmSchmidt!(p′, 4, 0)
    @test p′ == p

    p, dp = PlmSchmidt_d1(4, 0)
    dp′ = zeros(15)
    PlmSchmidt_d1!(p′, dp′, 4, 0)
    @test p′ == p
    @test dp′ == dp

    p = PlSchmidt(4, 0)
    p′ = zeros(5)
    PlSchmidt!(p′, 4, 0)
    @test p′ == p

    p, dp = PlSchmidt_d1(4, 0)
    dp′ = zeros(5)
    PlSchmidt_d1!(p′, dp′, 4, 0)
    @test p′ == p
    @test dp′ == dp
end

@testset "PLegendreA" begin
    p = PLegendreA(4, 0)
    p′ = zeros(15)
    PLegendreA!(p′, 4, 0)
    @test p′ == p

    p, dp = PLegendreA_d1(4, 0)
    dp′ = zeros(15)
    PLegendreA_d1!(p′, dp′, 4, 0)
    @test p′ == p
    @test dp′ == dp

    p = PLegendre(4, 0)
    p′ = zeros(5)
    PLegendre!(p′, 4, 0)
    @test p′ == p

    p, dp = PLegendre_d1(4, 0)
    dp′ = zeros(5)
    PLegendre_d1!(p′, dp′, 4, 0)
    @test p′ == p
    @test dp′ == dp
end

@testset "PlmIndex" begin
    index = PlmIndex(4, 3)
    @test index == 14
end

################################################################################

# Spherical harmonic transforms

@testset "MakeGridDH, MakeGradientDH, and SHExpandDH" begin
    # Invent a random grid
    Random.seed!(0)
    n = 10
    griddh = randn(n, n)

    # Calculate coefficients
    cilm, lmax = SHExpandDH(griddh, n)
    @test lmax == 4

    # Go back to the grid (it will be low-pass filtered)
    griddh′, n′ = MakeGridDH(cilm, lmax)
    @test n′ == n

    # Calculate derivatives
    theta, phi = MakeGradientDH(cilm, lmax)
    @test size(theta) == size(phi)

    # Go back to coefficients again (they must be the same)
    cilm′, lmax′ = SHExpandDH(griddh, n)
    @test lmax′ == lmax
    @test cilm′ == cilm

    # Go back to the grid again (it must the be same this time)
    griddh″, n″ = MakeGridDH(cilm′, lmax′)
    @test n″ == n
    @test griddh″ == griddh′
end

@testset "MakeGridDHC and SHExpandDHC" begin
    # Invent a random grid
    Random.seed!(0)
    n = 10
    griddh = randn(Complex{Float64}, n, n)

    # Calculate coefficients
    cilm, lmax = SHExpandDHC(griddh, n)
    @test lmax == 4

    # Go back to the grid (it will be low-pass filtered)
    griddh′, n′ = MakeGridDHC(cilm, lmax)
    @test n′ == n

    # Go back to coefficients again (they must be the same)
    cilm′, lmax′ = SHExpandDHC(griddh, n)
    @test lmax′ == lmax
    @test cilm′ == cilm

    # Go back to the grid again (it must the be same this time)
    griddh″, n″ = MakeGridDHC(cilm′, lmax′)
    @test n″ == n
    @test griddh″ == griddh′
end

@testset "SHExpandGLQ and MakeGridGLQ" begin
    Random.seed!(0)

    lmax = 4
    zero, w = SHGLQ(nothing, lmax)

    # Invent a random grid
    gridglq = randn(lmax + 1, 2 * lmax + 1)

    # Calculate coefficients
    cilm = SHExpandGLQ(lmax, gridglq, w, nothing, zero)

    # Go back to the grid (it will be low-pass filtered)
    gridglq′ = MakeGridGLQ(cilm, lmax, nothing, zero)
    @test size(gridglq′) == size(gridglq)

    # Go back to coefficients again (they must be the same)
    cilm′ = SHExpandGLQ(lmax, gridglq′, w, nothing, zero)
    @test isapprox(cilm′, cilm)

    # Go back to the grid again (it must the be same this time)
    gridglq″ = MakeGridGLQ(cilm′, lmax, nothing, zero)
    @test size(gridglq″) == size(gridglq′)
    @test isapprox(gridglq″, gridglq′)

    latglq, longlq = GLQGridCoord(lmax)
    @test (length(latglq), length(longlq)) == size(gridglq)
end

@testset "SHExpandGLQC and MakeGridGLQC" begin
    Random.seed!(0)

    lmax = 4
    zero, w = SHGLQ(nothing, lmax)

    # Invent a random grid
    gridglq = randn(Complex{Cdouble}, lmax + 1, 2 * lmax + 1)

    # Calculate coefficients
    cilm = SHExpandGLQC(lmax, gridglq, w, nothing, zero)

    # Go back to the grid (it will be low-pass filtered)
    gridglq′ = MakeGridGLQC(cilm, lmax, nothing, zero)
    @test size(gridglq′) == size(gridglq)

    # Go back to coefficients again (they must be the same)
    cilm′ = SHExpandGLQC(lmax, gridglq′, w, nothing, zero)
    @test isapprox(cilm′, cilm)

    # Go back to the grid again (it must the be same this time)
    gridglq″ = MakeGridGLQC(cilm′, lmax, nothing, zero)
    @test size(gridglq″) == size(gridglq′)
    @test isapprox(gridglq″, gridglq′)
end

@testset "SHExpandLSQ" begin
    Random.seed!(0)
    lmax = 2
    n = (lmax + 1)^2
    lat = π * rand(n)
    lon = 2π * rand(n)
    d = randn(n)
    weights = Float64[1 for i in 1:n]
    cilm, chi2 = SHExpandLSQ(d, lat, lon, n, lmax; weights=weights)
    @test size(cilm) == (2, lmax + 1, lmax + 1)
    @test abs(chi2) < 10 * eps(Float64)
end

@testset "MakeGrid2d" begin
    Random.seed!(0)
    lmax = 4
    cilm = randn(2, lmax + 1, lmax + 1)
    interval = 15                # grid spacing in degrees
    grid = MakeGrid2d(cilm, lmax, interval)
    @test size(grid) == (180 ÷ interval + 1, 360 ÷ interval + 1)
end

@testset "MakeGridPoint" begin
    Random.seed!(0)
    # Choose random values
    lmax = 2
    nmax = (lmax + 1)^2
    lat = π * rand(nmax)
    lon = 2π * rand(nmax)
    values = randn(nmax)
    weights = Float64[1 for n in 1:nmax]

    # Expand
    cilm, chi2 = SHExpandLSQ(values, lat, lon, nmax, lmax; weights=weights)
    @test size(cilm) == (2, lmax + 1, lmax + 1)
    @test abs(chi2) < 10 * eps(Float64)

    # Convert back to values (they will be low-pass filtered)
    values′ = Float64[MakeGridPoint(cilm, lmax, lat[n], lon[n]) for n in 1:nmax]

    # Expand again (result must be the same)
    cilm′, chi2 = SHExpandLSQ(values′, lat, lon, nmax, lmax; weights=weights)
    @test abs(chi2) < 10 * eps(Float64)
    @test isapprox(cilm′, cilm)

    # Convert back to values (result must be the same this time)
    values″ = Float64[MakeGridPoint(cilm′, lmax, lat[n], lon[n])
                      for n in 1:nmax]
    @test isapprox(values″, values′)
end

@testset "MakeGridPointC" begin
    Random.seed!(0)
    # Choose random values
    lmax = 2
    cilm = randn(Complex{Float64}, 2, lmax + 1, lmax + 1)
    lat = π * rand()
    lon = 2π * rand()
    value = MakeGridPointC(cilm, lmax, lat, lon)
    @test value isa Complex{Float64}
end

@testset "SHMultiply" begin
    Random.seed!(0)
    # Choose random values
    lmax1 = 2
    lmax2 = 3
    cilm1 = randn(2, lmax1 + 1, lmax1 + 1)
    cilm2 = randn(2, lmax2 + 1, lmax2 + 1)
    cilmout = SHMultiply(cilm1, lmax1, cilm2, lmax2)
    @test size(cilmout) == (2, lmax1 + lmax2 + 1, lmax1 + lmax2 + 1)
end

################################################################################

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
