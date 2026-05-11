@testset "physics" begin
    material = MaterialParameters(lambda = 1.0, mu = 1.0, gc = 1.0, ell = 0.1)
    strain = [1.0 0.0; 0.0 -0.5]

    values, vectors = spectral_decomposition(strain)
    @test length(values) == 2
    @test vectors * Diagonal(values) * vectors' ≈ strain
    @test tensile_energy_density(strain, material) ≈ 1.5
    @test update_history(0.25, strain, material) ≈ 1.5
end
