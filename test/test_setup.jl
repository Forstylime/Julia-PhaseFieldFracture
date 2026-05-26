import Ferrite

@testset "square tension setup" begin
    setup = setup_square_tension(cells = (2, 2), final_displacement = 0.01)

    @test Ferrite.getncells(setup.grid) == 4
    @test Ferrite.ndofs(setup.dh_u) == 18
    @test Ferrite.ndofs(setup.dh_d) == 9
    @test setup.crack_nodes isa Vector{Int}
end
