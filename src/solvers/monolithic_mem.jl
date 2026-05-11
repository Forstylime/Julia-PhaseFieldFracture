Base.@kwdef struct MonolithicMEMOptions
    maxiter::Int = 50
    tolerance::Float64 = 1e-8
    active_set_tolerance::Float64 = 1e-10
end

function solve_monolithic_mem(problem; options::MonolithicMEMOptions = MonolithicMEMOptions())
    return (; problem, options, converged = false, iterations = 0)
end
