@testset "assembly" begin
    material = MaterialParameters(lambda = 1.0, mu = 1.0, gc = 1.0, ell = 0.1)
    K = zeros(2, 2)
    R = zeros(2)

    K_out, R_out = assemble_system!(K, R, nothing, nothing, material)
    @test K_out === K
    @test R_out === R
    @test all(iszero, K)
    @test all(iszero, R)
end
