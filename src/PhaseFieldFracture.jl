module PhaseFieldFracture

using LinearAlgebra
using Serialization
using SparseArrays
using Ferrite
using FerriteGmsh
using Tensors

# physics
include("physics/constitutive.jl")
include("physics/energies.jl")

# fem setup
include("fem/setup.jl")
include("fem/assembly.jl")

# solvers
include("solvers/staggered.jl")
include("solvers/sem.jl")
include("solvers/crisfield.jl")

# utilities
include("utils/utils_fun.jl")
include("utils/meshes.jl")

export
    PhaseFieldMaterial,
    strain_spectral_split,
    tensile_energy_density,
    refine_range,
    refine_grid!,
    TensionSetup,
    MonolithicTensionSetup,
    make_square_tension_grid,
    create_l_shape_grid,
    create_staggered_dofhandlers,
    create_monolithic_dofhandler,
    create_displacement_constraints,
    create_phase_field_constraints,
    create_monolithic_constraints,
    initial_crack_nodes,
    setup_square_tension,
    setup_l_tension,
    setup_l_tension_monolithic,
    assemble_u!,
    assemble_d!,
    update_history!,
    elastic_energy,
    surface_energy,
    total_energy,
    elastic_energy_monolithic,
    surface_energy_monolithic,
    total_energy_monolithic,
    solve_staggered,
    solve_sem,
    solve_crisfield,
    assemble_mass_matrix_d!,
    assemble_monolithic!,
    compute_g,
    update_history_mono!,
    adapt_rho!
end
