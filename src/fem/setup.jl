Base.@kwdef struct SquareTensionSetup{G,DHU,DHD,CHU,CHD,N}
    grid::G
    dh_u::DHU
    dh_d::DHD
    ch_u::CHU
    ch_d::CHD
    crack_nodes::N
end

function make_square_tension_grid(cells::NTuple{2,Int} = (50, 50))
    return Ferrite.generate_grid(Ferrite.Quadrilateral, cells)
end

function create_staggered_dofhandlers(grid)
    dh_u = Ferrite.DofHandler(grid)
    ip_u = Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2
    Ferrite.add!(dh_u, :u, ip_u)
    Ferrite.close!(dh_u)

    dh_d = Ferrite.DofHandler(grid)
    ip_d = Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()
    Ferrite.add!(dh_d, :d, ip_d)
    Ferrite.close!(dh_d)

    return dh_u, dh_d
end

function create_displacement_constraints(dh_u, grid; top_displacement = 0.0)
    ch_u = Ferrite.ConstraintHandler(dh_u)

    bottom = Ferrite.getfacetset(grid, "bottom")
    top = Ferrite.getfacetset(grid, "top")

    Ferrite.add!(
        ch_u,
        Ferrite.Dirichlet(:u, bottom, (x, t) -> zeros(2), [1, 2]),
    )
    Ferrite.add!(
        ch_u,
        Ferrite.Dirichlet(:u, top, (x, t) -> t * top_displacement, 2),
    )

    Ferrite.close!(ch_u)
    Ferrite.update!(ch_u, 0.0)
    return ch_u
end

function create_phase_field_constraints(dh_d)
    ch_d = Ferrite.ConstraintHandler(dh_d)
    Ferrite.close!(ch_d)
    Ferrite.update!(ch_d, 0.0)
    return ch_d
end

function initial_crack_nodes(
    grid;
    y = 0.0,
    x_min = -1.0,
    x_max = 0.0,
    half_width = 1e-10,
)
    nodes = Int[]
    for (i, node) in pairs(grid.nodes)
        xcoord = node.x[1]
        ycoord = node.x[2]
        if x_min <= xcoord <= x_max && abs(ycoord - y) <= half_width
            push!(nodes, i)
        end
    end
    return nodes
end

function setup_square_tension(;
    cells::NTuple{2,Int} = (50, 50),
    top_displacement = 0.0,
    crack_y = 0.0,
    crack_x_min = -1.0,
    crack_x_max = 0.0,
    crack_half_width = 1e-10,
)
    grid = make_square_tension_grid(cells)
    dh_u, dh_d = create_staggered_dofhandlers(grid)
    ch_u = create_displacement_constraints(dh_u, grid; top_displacement)
    ch_d = create_phase_field_constraints(dh_d)
    crack_nodes = initial_crack_nodes(
        grid;
        y = crack_y,
        x_min = crack_x_min,
        x_max = crack_x_max,
        half_width = crack_half_width,
    )

    return SquareTensionSetup(; grid, dh_u, dh_d, ch_u, ch_d, crack_nodes)
end
