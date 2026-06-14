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
    # --- struct ---
    PhaseFieldMaterial,
    # --- setups ---
    setup_tension_monolithic,
    setup_l_tension,
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
    compute_g,
    get_right_dofs,
    compute_reaction_forces
end
