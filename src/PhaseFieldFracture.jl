module PhaseFieldFracture

using LinearAlgebra
using SparseArrays
using Ferrite
using Tensors

# include("physics/energies.jl")
include("physics/constitutive.jl")
include("physics/history.jl")

include("fem/elements.jl")
include("fem/quadrature.jl")
include("fem/setup.jl")
include("fem/assembly.jl")

include("solvers/monolithic_mem.jl")
include("solvers/staggered.jl")
include("solvers/arclength.jl")
include("solvers/block_solver.jl")

include("utils/meshes.jl")
include("utils/io.jl")
include("utils/metrics.jl")

export MaterialParameters,
    PhaseFieldMaterial,
    strain_spectral_split,
    tensile_energy_density,
    SquareTensionSetup,
    refine_range,
    refine_grid!,
    make_square_tension_grid,
    create_staggered_dofhandlers,
    create_displacement_constraints,
    create_phase_field_constraints,
    initial_crack_nodes,
    setup_square_tension,
    assemble_u!,
    assemble_d!,
    update_history!,
    solve_staggered
end
