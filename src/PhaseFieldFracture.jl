module PhaseFieldFracture

using LinearAlgebra
using SparseArrays
import Ferrite

include("physics/energies.jl")
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
    spectral_decomposition,
    stress,
    tensile_energy_density,
    fracture_energy_density,
    update_history,
    Quad4,
    Tri3,
    gauss_rule,
    SquareTensionSetup,
    make_square_tension_grid,
    create_staggered_dofhandlers,
    create_displacement_constraints,
    create_phase_field_constraints,
    initial_crack_nodes,
    setup_square_tension,
    assemble_system!,
    MonolithicMEMOptions,
    solve_monolithic_mem,
    StaggeredOptions,
    solve_staggered,
    ArcLengthOptions,
    solve_arclength,
    solve_block_system,
    make_l_shape_mesh,
    make_ct_specimen_mesh,
    write_results,
    SimulationMetrics,
    record_step!

end
