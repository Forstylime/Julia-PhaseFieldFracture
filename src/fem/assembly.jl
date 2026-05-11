function assemble_system!(K, R, state, mesh, material::MaterialParameters)
    fill!(K, zero(eltype(K)))
    fill!(R, zero(eltype(R)))
    return K, R
end
