using PhaseFieldFracture

function main(args = ARGS)
    hs = isempty(args) ? [0.1, 0.05, 0.025] : parse.(Float64, args)
    metrics = SimulationMetrics()
    for h in hs
        mesh = make_l_shape_mesh(h = h)
        result = solve_monolithic_mem((; mesh))
        record_step!(metrics; iterations = result.iterations)
        @info "Scalability placeholder run" h result.converged result.iterations
    end
    return metrics
end

abspath(PROGRAM_FILE) == @__FILE__ && main()
