mutable struct SimulationMetrics
    n_steps::Int
    n_iter::Int
    walltime::Float64
    f_peak::Float64
end

SimulationMetrics() = SimulationMetrics(0, 0, 0.0, -Inf)

function record_step!(metrics::SimulationMetrics; iterations = 0, walltime = 0.0, force = -Inf)
    metrics.n_steps += 1
    metrics.n_iter += iterations
    metrics.walltime += walltime
    metrics.f_peak = max(metrics.f_peak, force)
    return metrics
end
