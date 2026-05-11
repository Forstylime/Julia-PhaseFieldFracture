function spectral_decomposition(strain)
    F = eigen(Symmetric(strain))
    return F.values, F.vectors
end

function positive_strain(strain)
    values, vectors = spectral_decomposition(strain)
    positive_values = Diagonal(map(v -> max(v, zero(v)), values))
    return vectors * positive_values * vectors'
end

function stress(strain, d, material::MaterialParameters)
    eps_pos = positive_strain(strain)
    tr_pos = tr(eps_pos)
    sigma_pos = material.lambda * tr_pos * I + 2 * material.mu * eps_pos
    return residual_degradation(d, material.kappa) * sigma_pos
end
