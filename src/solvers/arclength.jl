Base.@kwdef struct ArcLengthOptions
    maxiter::Int = 50
    tolerance::Float64 = 1e-8
    radius::Float64 = 1.0
    constraint::Symbol = :crisfield
end

function solve_arclength(problem; options::ArcLengthOptions = ArcLengthOptions())
    return (; problem, options, converged = false, iterations = 0)
end
