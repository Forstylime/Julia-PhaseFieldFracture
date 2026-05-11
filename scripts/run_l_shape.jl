using PhaseFieldFracture

function main(args = ARGS)
    h = isempty(args) ? 0.05 : parse(Float64, args[1])
    mesh = make_l_shape_mesh(h = h)
    result = solve_monolithic_mem((; mesh))
    @info "Finished L-shape placeholder run" h result.converged result.iterations
    return result
end

abspath(PROGRAM_FILE) == @__FILE__ && main()
