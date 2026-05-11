Base.@kwdef struct StaggeredOptions
    maxiter::Int = 100
    tolerance::Float64 = 1e-8
end

function solve_staggered(problem; options::StaggeredOptions = StaggeredOptions())
    return (; problem, options, converged = false, iterations = 0)
end
