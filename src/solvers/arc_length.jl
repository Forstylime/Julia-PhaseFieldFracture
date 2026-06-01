using Ferrite
using Tensors
using LinearAlgebra
using SparseArrays

function update_history_coupled!(
    H::Vector{Float64}, dh::DofHandler, a::Vector{Float64},
    mat::PhaseFieldMaterial, cv_u::CellValues
)
    qp_count = 1
    dofs_u_range = dof_range(dh, :u)
    for cell in CellIterator(dh)
        reinit!(cv_u, cell)
        a_loc = a[celldofs(cell)]
        u_loc = a_loc[dofs_u_range]
        for q_point in 1:getnquadpoints(cv_u)
            ε_q = function_symmetric_gradient(cv_u, q_point, u_loc)
            ψ_plus = tensile_energy_density(ε_q, mat)
            H[qp_count] = max(H[qp_count], ψ_plus)
            qp_count += 1
        end
    end
end

function elastic_energy_coupled(
    dh::DofHandler, a::Vector{Float64},
    mat::PhaseFieldMaterial, cv_u::CellValues, cv_d::CellValues
)
    energy = 0.0
    dofs_u_range = dof_range(dh, :u)
    dofs_d_range = dof_range(dh, :d)
    for cell in CellIterator(dh)
        reinit!(cv_u, cell)
        reinit!(cv_d, cell)
        a_loc = a[celldofs(cell)]
        u_loc = a_loc[dofs_u_range]
        d_loc = a_loc[dofs_d_range]
        for q_point in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, q_point)
            ε_q = function_symmetric_gradient(cv_u, q_point, u_loc)
            d_q = function_value(cv_d, q_point, d_loc)
            energy += elastic_energy_density(ε_q, d_q, mat) * dΩ
        end
    end
    return energy
end

function surface_energy_coupled(
    dh::DofHandler, a::Vector{Float64},
    mat::PhaseFieldMaterial, cv_d::CellValues
)
    energy = 0.0
    dofs_d_range = dof_range(dh, :d)
    for cell in CellIterator(dh)
        reinit!(cv_d, cell)
        a_loc = a[celldofs(cell)]
        d_loc = a_loc[dofs_d_range]
        for q_point in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, q_point)
            d_q = function_value(cv_d, q_point, d_loc)
            ∇d_q = function_gradient(cv_d, q_point, d_loc)
            energy += (mat.gc / (2 * mat.l)) * (d_q^2 + mat.l^2 * (∇d_q ⋅ ∇d_q)) * dΩ
        end
    end
    return energy
end

function assemble_coupled!(
    K::SparseMatrixCSC, R::Vector{Float64},
    dh::DofHandler, a::Vector{Float64},
    H::Vector{Float64}, mat::PhaseFieldMaterial,
    cv_u::CellValues, cv_d::CellValues
)
    assembler = start_assemble(K, R)
    n_basefuncs_u = getnbasefunctions(cv_u)
    n_basefuncs_d = getnbasefunctions(cv_d)
    n_basefuncs = n_basefuncs_u + n_basefuncs_d
    
    Ke = zeros(n_basefuncs, n_basefuncs)
    Re = zeros(n_basefuncs)
    
    dofs_u_range = dof_range(dh, :u)
    dofs_d_range = dof_range(dh, :d)

    ψ0_e(ε) = (mat.λ / 2) * tr(ε)^2 + mat.μ * tr(ε ⋅ ε)
    ψ_plus_fun(ε) = tensile_energy_density(ε, mat)

    qp_count = 1
    for cell in CellIterator(dh)
        reinit!(cv_u, cell)
        reinit!(cv_d, cell)
        
        a_loc = a[celldofs(cell)]
        u_loc = a_loc[dofs_u_range]
        d_loc = a_loc[dofs_d_range]
        
        fill!(Ke, 0.0)
        fill!(Re, 0.0)
        
        for q_point in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, q_point)
            
            ε_q = function_symmetric_gradient(cv_u, q_point, u_loc)
            d_q = function_value(cv_d, q_point, d_loc)
            ∇d_q = function_gradient(cv_d, q_point, d_loc)
            
            # Kinematics & Constitutive
            σ0 = Tensors.gradient(ψ0_e, ε_q)
            ℂ0 = Tensors.hessian(ψ0_e, ε_q)
            
            ψ0_plus = ψ_plus_fun(ε_q)
            σ0_plus = Tensors.gradient(ψ_plus_fun, ε_q)
            
            # History
            H_prev = H[qp_count]
            if ψ0_plus > H_prev
                H_q = ψ0_plus
                ∂H_∂ε = σ0_plus
            else
                H_q = H_prev
                ∂H_∂ε = zero(σ0_plus)
            end
            
            g_q = (1.0 - d_q)^2 + mat.k_tol
            σ = g_q * σ0
            
            # 1. u equations
            for (i, i_u) in enumerate(dofs_u_range)
                δu = shape_symmetric_gradient(cv_u, q_point, i)
                Re[i_u] += (σ ⊡ δu) * dΩ
                
                # u-u block
                for (j, j_u) in enumerate(dofs_u_range)
                    Δu = shape_symmetric_gradient(cv_u, q_point, j)
                    Ke[i_u, j_u] += (δu ⊡ (g_q * ℂ0) ⊡ Δu) * dΩ
                end
                
                # u-d block
                for (j, j_d) in enumerate(dofs_d_range)
                    Δd = shape_value(cv_d, q_point, j)
                    Ke[i_u, j_d] += (δu ⊡ (-2.0 * (1.0 - d_q) * σ0) * Δd) * dΩ
                end
            end
            
            # 2. d equations
            for (i, i_d) in enumerate(dofs_d_range)
                δd = shape_value(cv_d, q_point, i)
                ∇δd = shape_gradient(cv_d, q_point, i)
                
                Re[i_d] += (-2.0 * (1.0 - d_q) * H_q * δd + 
                             (mat.gc / mat.l) * d_q * δd + 
                             mat.gc * mat.l * (∇d_q ⋅ ∇δd)) * dΩ
                             
                # d-u block
                for (j, j_u) in enumerate(dofs_u_range)
                    Δu = shape_symmetric_gradient(cv_u, q_point, j)
                    Ke[i_d, j_u] += (δd * (-2.0 * (1.0 - d_q)) * (∂H_∂ε ⊡ Δu)) * dΩ
                end
                
                # d-d block
                for (j, j_d) in enumerate(dofs_d_range)
                    Δd = shape_value(cv_d, q_point, j)
                    ∇Δd = shape_gradient(cv_d, q_point, j)
                    Ke[i_d, j_d] += ((2.0 * H_q + mat.gc / mat.l) * δd * Δd + 
                                      mat.gc * mat.l * (∇δd ⋅ ∇Δd)) * dΩ
                end
            end
            
            qp_count += 1
        end
        assemble!(assembler, celldofs(cell), Ke, Re)
    end
end

"""
Consistent Monolithic Newton Solver (Bordering Algorithm)
"""
function solve_arc_length(
    setup::TensionSetup, mat::PhaseFieldMaterial;
    n_steps_max::Int = 1000, 
    ρ_init::Float64 = 0.05, 
    tol_newton::Float64 = 1e-6, 
    max_newton_iter::Int = 15,
    t_max::Float64 = 1.0,
    output_freq::Int = 5
)
    grid = setup.grid
    
    # 1. Coupled Setup
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral, 1}()^2)
    add!(dh, :d, Lagrange{RefQuadrilateral, 1}())
    close!(dh)
    
    ch = ConstraintHandler(dh)
    top = Ferrite.getfacetset(grid, "top")
    right = Ferrite.getfacetset(grid, "right")
    
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, top, (x, t) -> zeros(2), [1, 2]))
    # For constraint extraction we map to `t * final_displacement`
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, right, (x, t) -> t * setup.final_displacement, 2))
    
    if !isempty(setup.crack_nodes)
        Ferrite.add!(ch, Ferrite.Dirichlet(:d, Set(setup.crack_nodes), (x, t) -> 1.0))
    end
    close!(ch)
    update!(ch, 1.0)
    
    ndofs_total = ndofs(dh)
    
    # Extract u_dofs exactly
    dofs_u_range = dof_range(dh, :u)
    u_dofs_list = Int[]
    dofs_d_range = dof_range(dh, :d)
    d_dofs_list = Int[]
    
    for cell_id in 1:getncells(grid)
        dofs = celldofs(dh, cell_id)
        append!(u_dofs_list, dofs[dofs_u_range])
        append!(d_dofs_list, dofs[dofs_d_range])
    end
    u_dofs = unique(u_dofs_list)
    d_dofs = unique(d_dofs_list)
    
    # 2. Setup f_hat from ConstraintHandler at t=1.0
    f_hat = zeros(ndofs_total)
    update!(ch, 1.0)
    apply!(f_hat, ch)
    
    # Ensure phase field DOFs don't scale with λ
    for dof in d_dofs
        f_hat[dof] = 0.0
    end
    
    # Now reset ConstraintHandler to t=0.0 for the Newton loop
    update!(ch, 0.0) 
    
    # 3. Allocations
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())
    
    n_qpoints = getncells(grid) * getnquadpoints(cv_u)
    H_history = zeros(n_qpoints)
    
    K_aa = allocate_matrix(dh)
    r_a = zeros(ndofs_total)
    
    a_prev = zeros(ndofs_total)
    update!(ch, 1.0)
    apply!(a_prev, ch) # initial crack and displacement
    fill!(a_prev, 0.0)
    update!(ch, 0.0)
    apply!(a_prev, ch) # initial crack = 1.0, initial disp = 0.0
    
    # Equilibrate initial crack at λ=0
    println("--- Equilibrating initial crack profile ---")
    for iter in 1:15
        assemble_coupled!(K_aa, r_a, dh, a_prev, H_history, mat, cv_u, cv_d)
        r_a_free = copy(r_a)
        apply_zero!(K_aa, r_a_free, ch)
        norm_r = norm(r_a_free)
        println("  Init Iter $iter: norm_r = $norm_r")
        if norm_r < tol_newton
            break
        end
        lu_K = lu(K_aa)
        Δa = lu_K \ (-r_a_free)
        a_prev .+= Δa
    end
    update_history_coupled!(H_history, dh, a_prev, mat, cv_u)
    println("--- Initial crack equilibrated ---")
    
    a_n = copy(a_prev)
    λ_prev = 0.0
    λ_n = 0.0
    
    Δa_u_prev = zeros(length(u_dofs))
    
    ρ = ρ_init
    
    displacements = Float64[0.0]
    reaction_forces = Float64[0.0]
    elastic_energies = Float64[0.0]
    surface_energies = Float64[0.0]
    lambdas = Float64[0.0]
    
    mkpath("data/sims/arc_length")
    
    VTKGridFile("data/sims/arc_length/step_0", dh) do vtk
        write_solution(vtk, dh, a_prev)
    end
    
    step = 1
    println("--- Starting Monolithic Arc-Length Solver ---")
    
    while λ_prev < t_max && step <= n_steps_max
        println("Step $step | λ_prev = $(round(λ_prev, digits=4)) | ρ = $(round(ρ, digits=5))")
        
        # Predictor: we need Δa_II which is the elastic response to f_hat
        assemble_coupled!(K_aa, r_a, dh, a_prev, H_history, mat, cv_u, cv_d)
        
        # Exact solution for prescribed displacement f_hat
        Δa_II = copy(f_hat)
        r_II = -(K_aa * Δa_II)
        apply_zero!(K_aa, r_II, ch)
        
        lu_K = lu(K_aa)
        Δa_II .+= lu_K \ r_II
        
        Δa_II_u = Δa_II[u_dofs]
        
        println("  Predictor: norm(f_hat) = $(norm(f_hat)), norm(Δa_II_u) = $(norm(Δa_II_u))")
        
        if step == 1
            Δλ_pred = abs(ρ / norm(Δa_II_u))
        else
            Δλ_pred = ρ / norm(Δa_II_u)
            if dot(Δa_II_u, Δa_u_prev) < 0
                Δλ_pred = -Δλ_pred
            end
        end
        
        a_n .= a_prev .+ Δλ_pred .* Δa_II
        λ_n = λ_prev + Δλ_pred
        
        # --- Corrector ---
        converged = false
        iter = 0
        
        while iter < max_newton_iter
            iter += 1
            
            assemble_coupled!(K_aa, r_a, dh, a_n, H_history, mat, cv_u, cv_d)
            K_aa_unconstrained = copy(K_aa)
            
            r_a_free = copy(r_a)
            apply_zero!(K_aa, r_a_free, ch)
            
            Δu_n = a_n[u_dofs] .- a_prev[u_dofs]
            fc = dot(Δu_n, Δu_n) - ρ^2
            
            norm_r = norm(r_a_free)
            
            println("  Iter $iter: norm_r = $norm_r, abs(fc) = $(abs(fc))")
            
            if norm_r < tol_newton && (abs(fc) / (ρ^2) < tol_newton)
                converged = true
                break
            end
            
            lu_K_n = lu(K_aa)
            Δa_I = lu_K_n \ (-r_a_free)
            
            Δa_II = copy(f_hat)
            r_II = -(K_aa_unconstrained * Δa_II)
            # Since K_aa is already zeroed, we can't use apply_zero! to r_II using K_aa!
            # Wait, apply_zero! only modifies K and r. If we don't care about K here, we can just zero out r_II at dirichlet dofs!
            for (i, dof) in enumerate(ch.prescribed_dofs)
                r_II[dof] = 0.0
            end
            Δa_II .+= lu_K_n \ r_II
            
            Δa_I_u = Δa_I[u_dofs]
            Δa_II_u = Δa_II[u_dofs]
            
            K_λa_Δa_I = 2.0 * dot(Δu_n, Δa_I_u)
            K_λa_Δa_II = 2.0 * dot(Δu_n, Δa_II_u)
            
            Δλ = (-fc - K_λa_Δa_I) / K_λa_Δa_II
            println("    Δλ = $Δλ (fc=$fc, K_λa_Δa_I=$(K_λa_Δa_I), K_λa_Δa_II=$(K_λa_Δa_II))")
            
            a_n .+= Δa_I .+ Δλ .* Δa_II
            λ_n += Δλ
        end
        
        if converged
            println(" -> Converged in $iter iter. λ = $(round(λ_n, digits=4))")
            Δa_u_prev .= a_n[u_dofs] .- a_prev[u_dofs]
            a_prev .= a_n
            λ_prev = λ_n
            
            update_history_coupled!(H_history, dh, a_prev, mat, cv_u)
            
            assemble_coupled!(K_aa, r_a, dh, a_prev, H_history, mat, cv_u, cv_d)
            f_total = 0.0
            for (i, dof) in enumerate(ch.prescribed_dofs)
                if f_hat[dof] != 0.0
                    f_total += r_a[dof]
                end
            end
            
            push!(lambdas, λ_prev)
            push!(displacements, λ_prev * setup.final_displacement)
            push!(reaction_forces, f_total)
            push!(elastic_energies, elastic_energy_coupled(dh, a_prev, mat, cv_u, cv_d))
            push!(surface_energies, surface_energy_coupled(dh, a_prev, mat, cv_d))
            
            if step % output_freq == 0 || λ_prev >= t_max
                VTKGridFile("data/sims/arc_length/step_$step", dh) do vtk
                    write_solution(vtk, dh, a_prev)
                end
            end
            
            if iter < 4
                ρ *= 1.2
            elseif iter > 10
                ρ *= 0.5
            end
            
            step += 1
        else
            println(" -> Failed to converge! Cutting step size.")
            ρ *= 0.5
            if ρ < 1e-4
                println("ρ is too small. Aborting.")
                break
            end
        end
    end
    
    return lambdas, displacements, reaction_forces, elastic_energies, surface_energies
end
