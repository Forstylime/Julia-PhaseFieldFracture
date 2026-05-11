@testset "active set" begin
    options = MonolithicMEMOptions(active_set_tolerance = 1e-9)
    result = solve_monolithic_mem((; element = Quad4()); options)

    @test result.options === options
    @test result.iterations == 0
    @test result.converged == false
end
