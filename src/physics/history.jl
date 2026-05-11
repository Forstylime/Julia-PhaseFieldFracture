update_history(previous_H, strain, material::MaterialParameters) =
    max(previous_H, tensile_energy_density(strain, material))
