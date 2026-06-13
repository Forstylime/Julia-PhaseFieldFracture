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
include("solvers/gamma.jl")
include("solvers/h1.jl")
include("solvers/L2.jl")
include("solvers/mem.jl")

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
    # --- setups ---
    setup_square_tension,
    setup_l_tension,
    setup_l_tension_monolithic,
    setup_l_tension_mem,
    # --- energies ---
    elastic_energy,
    surface_energy,
    total_energy,
    elastic_energy_monolithic,
    surface_energy_monolithic,
    total_energy_monolithic,
    # --- solvers ---
    solve_staggered,
    solve_sem,
    solve_crisfield,
    solve_gamma,
    solve_h1,
    solve_l2,
    solve_mem,
    # --- assembly utils ---
    assemble_u!,
    assemble_d!,
    assemble_mass_matrix_d!,
    assemble_monolithic!,
    # --- utils ---
    update_history!,
    update_history_mono!,
    compute_g,
    get_right_dofs,
    compute_reaction_forces,
    adapt_rho!,
    miehe_spectral_decomposition,
    evaluate_gamma_constraint
end
