using SHTOOLS
using Test

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
