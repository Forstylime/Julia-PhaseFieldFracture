struct MaterialParameters{T}
    lambda::T
    mu::T
    gc::T
    ell::T
    kappa::T
end

function MaterialParameters(; lambda, mu, gc, ell, kappa = nothing)
    T = promote_type(typeof(lambda), typeof(mu), typeof(gc), typeof(ell))
    kappa_value = isnothing(kappa) ? zero(T) : kappa
    return MaterialParameters{promote_type(T, typeof(kappa_value))}(
        lambda,
        mu,
        gc,
        ell,
        kappa_value,
    )
end

degradation(d) = (1 - d)^2
residual_degradation(d, kappa) = degradation(d) + kappa

function tensile_energy_density(strain, material::MaterialParameters)
    vals = eigvals(Symmetric(strain))
    positive_trace = sum(max(v, zero(v)) for v in vals)
    positive_square = sum(max(v, zero(v))^2 for v in vals)
    return 0.5 * material.lambda * positive_trace^2 + material.mu * positive_square
end

function fracture_energy_density(d, grad_d, material::MaterialParameters)
    return material.gc * (d^2 / (2 * material.ell) + material.ell * dot(grad_d, grad_d) / 2)
end
