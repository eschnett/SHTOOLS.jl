@testset "Wigner3j" begin
    w3j_prealloc = zeros(Cdouble, 5)
    for j2 in 0:2, j3 = 0:2, m2 in -j2:j2, m3 in -j3:j3
        jmax = j2 + j3
        for m1 in -jmax:jmax
            w3j, jmin, jmax = Wigner3j(j2, j3, m1, m2, m3)
            Wigner3j!(w3j_prealloc, j2, j3, m1, m2, m3)
            @test w3j_prealloc[1:(jmax - jmin + 1)] == w3j
        end
    end

    @testset "known values" begin
        for j in 0:20, m in -j:j
            w3j, jmin, jmax = Wigner3j(j, 0, m, -m, 0)
            @test length(w3j) == 1
            @test w3j[1] ≈ (-1)^(j-m)/√(2j+1)
        end
    end
end
