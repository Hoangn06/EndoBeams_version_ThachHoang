# ---------------------------------------------------------------
# This function ensures displacement, velocity, and acceleration states are 
# correctly set up for the next time step.
# ---------------------------------------------------------------

function initialize_global_corrector!(conf::BeamsConfiguration, state::SimulationState)
    for n in LazyRows(conf.nodes)
        # Reset incremental displacement for translational DOFs
        state.solⁿ⁺¹.ΔD[n.global_dofs_disp] .= 0  
        
        # Set displacement, velocity, and acceleration from nodal values
        state.solⁿ⁺¹.D[n.global_dofs_disp] .= n.u
        state.solⁿ⁺¹.Ḋ[n.global_dofs_disp] .= n.u̇
        state.solⁿ⁺¹.D̈[n.global_dofs_disp] .= n.ü
        
        # Reset incremental displacement for rotational DOFs
        state.solⁿ⁺¹.ΔD[n.global_dofs_rot] .= 0  
        
        # Set rotational displacement, velocity, and acceleration
        state.solⁿ⁺¹.D[n.global_dofs_rot] .= n.w
        state.solⁿ⁺¹.Ḋ[n.global_dofs_rot] .= n.ẇ
        state.solⁿ⁺¹.D̈[n.global_dofs_rot] .= n.ẅ
    end
end

# ---------------------------------------------------------------
# This function sets initial guesses for displacements, velocities, and accelerations.
# ---------------------------------------------------------------

function initialize_global_predictor!(conf::BeamsConfiguration, state::SimulationState)
    for n in LazyRows(conf.nodes)
        # Reset incremental displacement for translational DOFs
        state.solⁿ⁺¹.ΔD[n.global_dofs_disp] .= 0  
        
        # Set initial displacement, velocity, and acceleration for translational DOFs
        state.solⁿ⁺¹.D[n.global_dofs_disp] .= n.u
        state.solⁿ⁺¹.Ḋ[n.global_dofs_disp] .= n.u̇
        state.solⁿ⁺¹.D̈[n.global_dofs_disp] .= n.ü
        
        # Reset incremental displacement for rotational DOFs
        state.solⁿ⁺¹.ΔD[n.global_dofs_rot] .= 0  
        
        # Predict rotational displacement as zero but keep velocity and acceleration
        state.solⁿ⁺¹.D[n.global_dofs_rot] .= 0
        state.solⁿ⁺¹.Ḋ[n.global_dofs_rot] .= n.ẇ
        state.solⁿ⁺¹.D̈[n.global_dofs_rot] .= n.ẅ
    end
end 

# ---------------------------------------------------------------
# Computes the tangent matrix and residual forces in the predictor step.
# Uses Newmark parameters (α, β, γ) to construct the system.
# ---------------------------------------------------------------
function compute_tangent_and_residuals_predictor!(state::SimulationState, Δt, α, β, γ)

    # Compute tangent stiffness matrix (Ktan) at predictor step
    # Ktan = (1+α) * K + (γ/(βΔt)) * C + (1/(βΔt²)) * M
    @. state.solⁿ⁺¹.Ktan.nzval = (1 + α) * state.matricesⁿ.K.nzval + (γ / (β * Δt)) * state.matricesⁿ.C.nzval + (1 / (β * Δt^2)) * state.matricesⁿ.M.nzval
    
    # Compute residual vector (r) at predictor step
    # r = (1+α) * external forces - (internal forces + contact forces) - α * previous external forces
    @. state.solⁿ⁺¹.r = (1 + α) * state.forcesⁿ⁺¹.fᵉˣᵗ - state.forcesⁿ.Tⁱⁿᵗ - state.forcesⁿ.Tᶜ - state.forcesⁿ.Tᵏ - α * state.forcesⁿ.fᵉˣᵗ

    # Add contribution from damping matrix (C)
    # temp = (γ/β) * Ḋ - (Δt/2 * (2β - γ) / β) * D̈
    @. state.solⁿ⁺¹.temp = (γ / β) * state.solⁿ⁺¹.Ḋ - (Δt / 2 * (2β - γ) / β) * state.solⁿ⁺¹.D̈
    mul!(state.solⁿ⁺¹.r, state.matricesⁿ.C, state.solⁿ⁺¹.temp, 1, 1)  # r += C * temp

    # Add contribution from mass matrix (M)
    # temp = Δt * Ḋ + (Δt² / 2) * D̈
    @. state.solⁿ⁺¹.temp = Δt * state.solⁿ⁺¹.Ḋ + (Δt^2 / 2) * state.solⁿ⁺¹.D̈
    mul!(state.solⁿ⁺¹.r, state.matricesⁿ.M, state.solⁿ⁺¹.temp, 1 / (β * Δt^2), 1)  # r += (1 / (βΔt²)) * M * temp

end

# ---------------------------------------------------------------
# Computes the tangent matrix and residual forces in the corrector step.
# Used to refine the solution after the predictor step.
# ---------------------------------------------------------------
function compute_tangent_and_residuals_corrector!(state::SimulationState, Δt, α, β, γ)

    # Compute tangent stiffness matrix (Ktan) at corrector step (same as in predictor)
    @. state.solⁿ⁺¹.Ktan.nzval = (1 + α) * state.matricesⁿ⁺¹.K.nzval + (γ / (β * Δt)) * state.matricesⁿ⁺¹.C.nzval + (1 / (β * Δt^2)) * state.matricesⁿ⁺¹.M.nzval

    # Compute residual vector (r) at corrector step
    # r = (1+α) * (external + contact - internal forces) - α * (previous external + contact - internal forces) - previous contact force
    @. state.solⁿ⁺¹.r = (1 + α) * (state.forcesⁿ⁺¹.fᵉˣᵗ + state.forcesⁿ⁺¹.Tᶜ - state.forcesⁿ⁺¹.Tⁱⁿᵗ) - α * (state.forcesⁿ.fᵉˣᵗ + state.forcesⁿ.Tᶜ - state.forcesⁿ.Tⁱⁿᵗ) - state.forcesⁿ⁺¹.Tᵏ
    
end

# ---------------------------------------------------------------
# Extracts the residual and tangent stiffness matrix for the free DOFs.
# ---------------------------------------------------------------
function extract_free_dofs!(conf::SimulationConfiguration, state::SimulationState)
    
    # Extract the residual vector (r) for free DOFs
    @views state.solⁿ⁺¹.r_free .= state.solⁿ⁺¹.r[conf.bcs.free_dofs]

    # Extract the corresponding entries of the tangent stiffness matrix (Ktan)
    @views state.solⁿ⁺¹.Ktan_free.nzval .= state.solⁿ⁺¹.Ktan.nzval[state.matricesⁿ⁺¹.sparsity_free]
    
end

# ---------------------------------------------------------------
# Solves for the displacement increments at the free DOFs.
# ---------------------------------------------------------------
function solve_free_dofs!(conf::SimulationConfiguration, state::SimulationState, solver)
    
    # 🔹 Perform symbolic analysis on the reduced stiffness matrix
    # This step prepares the solver for the system K ΔD = r.
    analyze!(solver, state.solⁿ⁺¹.Ktan_free, state.solⁿ⁺¹.r_free)
    
    # Solve for displacement increments at free DOFs
    # The solver finds ΔD_free by solving Ktan_free * ΔD_free = r_free
    solve!(solver, state.solⁿ⁺¹.ΔD_free, state.solⁿ⁺¹.Ktan_free, state.solⁿ⁺¹.r_free)

    # Store the computed displacement increments in the full displacement vector
    # Only the free DOFs are updated, keeping fixed DOFs unchanged.
    state.solⁿ⁺¹.ΔD[conf.bcs.free_dofs] .= state.solⁿ⁺¹.ΔD_free
end

# ---------------------------------------------------------------
# Computes the norms of the residual and displacement increment
# during the corrector loop for convergence checking.
# ---------------------------------------------------------------
function compute_norms_corrector(conf::SimulationConfiguration, state::SimulationState)
    
    # Extract the sets of free and fixed DOFs
    free_dofs = conf.bcs.free_dofs
    fixed_dofs = conf.bcs.fixed_dofs
    
    # Compute the norm of the displacement increment (ΔD)
    ΔD_norm = norm(state.solⁿ⁺¹.ΔD_free)
    
    # Compute the residual norm for free DOFs
    res_norm = norm(state.solⁿ⁺¹.r_free)
    
    # Compute the total external force norm for free DOFs
    f_norm = norm(state.forcesⁿ⁺¹.fᵉˣᵗ[free_dofs] .+ state.forcesⁿ⁺¹.Tᶜ[free_dofs] .+ state.forcesⁿ⁺¹.Tᵏ[free_dofs])
    
    # Compute the total reaction force norm for fixed DOFs
    e_norm = norm(state.forcesⁿ⁺¹.fᵉˣᵗ[fixed_dofs] .+ state.forcesⁿ⁺¹.Tᶜ[fixed_dofs] .+ state.forcesⁿ⁺¹.Tᵏ[fixed_dofs])
    
    # Normalize the residual norm for better convergence checks
    if f_norm + e_norm > 1e-12    
        res_norm = res_norm / (f_norm + e_norm)
    end
    
    # Return the computed norms for convergence assessment
    return res_norm, ΔD_norm   
end

# -------------------------------------------------------------------
# Updates the global displacement, velocity, and acceleration 
# at the end of the predictor step based on the computed increments.
# -------------------------------------------------------------------
function compute_global_predictor!(state::SimulationState, Δt, β, γ)
    
    # Precompute frequently used coefficients
    γ_over_β = γ / β
    γ_over_βΔt = γ / (β * Δt)
    inv_βΔt² = 1 / (β * Δt^2)
    Δt_over_2β = Δt * (2 * β - γ) / (2 * β)
    
    # Update global predictor values using Newmark-beta formulas
    @inbounds for i in eachindex(state.solⁿ⁺¹.ΔD)
        Ḋⁿ = deepcopy(state.solⁿ⁺¹.Ḋ[i])  # Store previous velocity to avoid overwriting during updates
        state.solⁿ⁺¹.D[i] += state.solⁿ⁺¹.ΔD[i] 
        state.solⁿ⁺¹.Ḋ[i] += γ_over_βΔt * state.solⁿ⁺¹.ΔD[i] - γ_over_β * Ḋⁿ + Δt_over_2β * state.solⁿ⁺¹.D̈[i]  
        state.solⁿ⁺¹.D̈[i] += inv_βΔt² * (state.solⁿ⁺¹.ΔD[i] - Δt * Ḋⁿ - (Δt^2 / 2) * state.solⁿ⁺¹.D̈[i])  
    end
end

# -------------------------------------------------------------------
# Updates the global displacement, velocity, and acceleration 
# at the end of the corrector step based on the computed increments.
# -------------------------------------------------------------------
function compute_global_corrector!(state::SimulationState, Δt, β, γ)
    
    # Precompute frequently used coefficients 
    inv_βΔt² = 1 / (β * Δt^2)
    γ_over_βΔt = γ / (β * Δt)

    # 🔹 Update global corrector values
    @inbounds for i in eachindex(state.solⁿ⁺¹.ΔD)
        state.solⁿ⁺¹.D[i] += state.solⁿ⁺¹.ΔD[i]
        state.solⁿ⁺¹.Ḋ[i] += γ_over_βΔt * state.solⁿ⁺¹.ΔD[i]  
        state.solⁿ⁺¹.D̈[i] += inv_βΔt² * state.solⁿ⁺¹.ΔD[i]  
    end
end 

# -------------------------------------------------------------------
# Updates local vectors for different configurations 
# based on the global predictor solution.
# -------------------------------------------------------------------

function compute_local_predictor!(conf::BeamsConfiguration, state::SimulationState, γ=nothing, β=nothing, Δt=nothing)
    
    @inbounds for i in eachindex(conf.nodes)
        # Update displacement, velocity, and acceleration
        conf.nodes.u[i] = state.solⁿ⁺¹.D[conf.nodes.global_dofs_disp[i]]
        conf.nodes.u̇[i] = state.solⁿ⁺¹.Ḋ[conf.nodes.global_dofs_disp[i]]
        conf.nodes.ü[i] = state.solⁿ⁺¹.D̈[conf.nodes.global_dofs_disp[i]]
        
        # Update rotational degrees of freedom
        θ̃ = state.solⁿ⁺¹.D[conf.nodes.global_dofs_rot[i]]
        conf.nodes.w[i] = θ̃
        conf.nodes.ΔR[i] = rotation_matrix(θ̃)  # Compute incremental rotation matrix
        conf.nodes.R[i] = conf.nodes.ΔR[i] * conf.nodes.Rⁿ[i]  # Update total rotation

        # Compute angular velocity and acceleration
        ẇⁿ = conf.nodes.ẇⁿ[i]
        ẅⁿ = conf.nodes.ẅⁿ[i]
        conf.nodes.ẇ[i] = conf.nodes.ΔR[i] * (γ/(β*Δt)*θ̃ + (β-γ)/β*ẇⁿ + Δt*(β-γ/2)/β*ẅⁿ)     
        conf.nodes.ẅ[i] = conf.nodes.ΔR[i] * (1/(β*Δt^2)*θ̃ - 1/(β*Δt)*ẇⁿ - (1/(2*β)-1)*ẅⁿ) 
    end  
end

# -------------------------------------------------------------------
# Updates local vectors for different configurations 
# based on the computed corrector solution.
# -------------------------------------------------------------------

function compute_local_corrector!(conf::BeamsConfiguration, state::SimulationState, γ=nothing, β=nothing, Δt=nothing)
    
    @inbounds for i in eachindex(conf.nodes)
        # Update displacement, velocity, and acceleration
        conf.nodes.u[i] = state.solⁿ⁺¹.D[conf.nodes.global_dofs_disp[i]]
        conf.nodes.u̇[i] = state.solⁿ⁺¹.Ḋ[conf.nodes.global_dofs_disp[i]]
        conf.nodes.ü[i] = state.solⁿ⁺¹.D̈[conf.nodes.global_dofs_disp[i]]
        
        # Compute the incremental rotation and apply it
        conf.nodes.ΔR[i] = rotation_matrix(state.solⁿ⁺¹.ΔD[conf.nodes.global_dofs_rot[i]]) * conf.nodes.ΔR[i]

        # Compute new angular velocity and acceleration
        ẇⁿ = conf.nodes.ẇⁿ[i]
        ẅⁿ = conf.nodes.ẅⁿ[i]
        wⁿ⁺¹ = toangle(conf.nodes.ΔR[i])
        conf.nodes.ẇ[i] = conf.nodes.ΔR[i] * (γ/(β*Δt)*wⁿ⁺¹ + (β-γ)/β*ẇⁿ + Δt*(β-γ/2)/β*ẅⁿ)     
        conf.nodes.ẅ[i] = conf.nodes.ΔR[i] * (1/(β*Δt^2)*wⁿ⁺¹ - 1/(β*Δt)*ẇⁿ - (1/(2*β)-1)*ẅⁿ)
        
        # Update total rotation
        conf.nodes.R[i] = conf.nodes.ΔR[i] * conf.nodes.Rⁿ[i] 
    end  
end 

# -------------------------------------------------------------------
# Updates nodal values and forces for the next time step 
# upon convergence for different configurations.
# -------------------------------------------------------------------

@inline function copy_superelastic_state!(dst, src)
    dst.xi_S .= src.xi_S
    dst.F .= src.F
    dst.H_AS .= src.H_AS
    dst.H_SA .= src.H_SA
    dst.H_SS .= src.H_SS
    dst.eps .= src.eps
    dst.sigma .= src.sigma
    return nothing
end

@inline function copy_plastic_state!(dst, src)
    dst.sigma .= src.sigma
    dst.epsilon .= src.epsilon
    dst.epsilon_pl .= src.epsilon_pl
    dst.e_pl .= src.e_pl
    dst.F .= src.F
    return nothing
end

function compute_vm_superelastic(sigma, eps)
    n_gauss_points = size(sigma, 2)
    stress_VM = Vector{Float64}(undef, n_gauss_points)
    strain_VM = Vector{Float64}(undef, n_gauss_points)

    for gp in 1:n_gauss_points
        sigma_gp = sigma[:, gp]
        e_gp = eps[:, gp]
        stress_VM[gp] = sqrt(sigma_gp[1]^2 + 3 * sigma_gp[2]^2 + 3 * sigma_gp[3]^2)
        strain_VM[gp] = sqrt(e_gp[1]^2 + 3 * e_gp[2]^2 + 3 * e_gp[3]^2)
    end

    return stress_VM, strain_VM
end

function compute_vm_plastic(sigma, epsilon)
    n_gauss_points = size(sigma, 2)
    stress_VM = Vector{Float64}(undef, n_gauss_points)
    strain_VM = Vector{Float64}(undef, n_gauss_points)

    for gp in 1:n_gauss_points
        sigma_gp = sigma[:, gp]
        eps_gp = epsilon[:, gp]
        stress_VM[gp] = sqrt(sigma_gp[1]^2 + 3 * sigma_gp[2]^2 + 3 * sigma_gp[3]^2)
        strain_VM[gp] = sqrt(eps_gp[1]^2 + 3 * eps_gp[2]^2 + 3 * eps_gp[3]^2)
    end

    return stress_VM, strain_VM
end

function update_converged!(conf::BeamsConfiguration, state::SimulationState;
    t::Float64=0.0,
    rotation_storage=nothing,
    displacement_storage=nothing,
    stress_strain_storage=nothing,
)
    @inbounds for i in eachindex(conf.nodes)
        # Update nodal displacements, velocities, and accelerations to the converged values
        conf.nodes.uⁿ[i] = conf.nodes.u[i]
        conf.nodes.u̇ⁿ[i] = conf.nodes.u̇[i]
        conf.nodes.üⁿ[i] = conf.nodes.ü[i]
        conf.nodes.wⁿ[i] = conf.nodes.w[i]
        conf.nodes.ẇⁿ[i] = conf.nodes.ẇ[i]
        conf.nodes.ẅⁿ[i] = conf.nodes.ẅ[i]
        conf.nodes.Rⁿ[i] = conf.nodes.R[i]
        conf.nodes.ΔRⁿ[i] = conf.nodes.ΔR[i]
    end

    # Accumulate rotation data if storage is provided
    if rotation_storage !== nothing
        accumulate_rotation_data!(rotation_storage, conf, t)
    end
    
    # Accumulate displacement data if storage is provided
    if displacement_storage !== nothing
        accumulate_displacement_data!(displacement_storage, conf, t)
    end
    
    # Update the force vectors for the current time step
    state.forcesⁿ.Tⁱⁿᵗ .= state.forcesⁿ⁺¹.Tⁱⁿᵗ
    state.forcesⁿ.Tᵏ .= state.forcesⁿ⁺¹.Tᵏ
    state.forcesⁿ.Tᶜ .= state.forcesⁿ⁺¹.Tᶜ
    state.forcesⁿ.fᵉˣᵗ .= state.forcesⁿ⁺¹.fᵉˣᵗ
    state.matricesⁿ.K .= state.matricesⁿ⁺¹.K
    state.matricesⁿ.C .= state.matricesⁿ⁺¹.C 
    state.matricesⁿ.M .= state.matricesⁿ⁺¹.M

    # Update constitutive states and optionally store VM stress/strain
    ms = state.material_states.by_beam
    for beam in LazyRows(conf.beams)
        if beam.material == :superelastic
            se_state = ms[beam.ind]::SuperelasticStates
            copy_superelastic_state!(se_state.superelasticⁿ, se_state.superelasticⁿ⁺¹)

            # Accumulate stress/strain VM data if storage is provided
            if stress_strain_storage !== nothing
                stress_VM, strain_VM = compute_vm_superelastic(se_state.superelasticⁿ.sigma, se_state.superelasticⁿ.eps)
                accumulate_stress_strain_data!(stress_strain_storage, state, t, stress_VM, strain_VM)
            end

        elseif beam.material == :plastic
            pl_state = ms[beam.ind]::PlasticStates
            copy_plastic_state!(pl_state.Plasticⁿ, pl_state.Plasticⁿ⁺¹)

            # Accumulate stress/strain VM data if storage is provided
            if stress_strain_storage !== nothing
                stress_VM, strain_VM = compute_vm_plastic(pl_state.Plasticⁿ.sigma, pl_state.Plasticⁿ.epsilon)
                accumulate_stress_strain_data!(stress_strain_storage, state, t, stress_VM, strain_VM)
            end

        elseif beam.material == :DFT_SEP
            dft_sep_state = ms[beam.ind]::DFT_SEP_States
            copy_superelastic_state!(dft_sep_state.superelasticⁿ, dft_sep_state.superelasticⁿ⁺¹)
            copy_plastic_state!(dft_sep_state.Plasticⁿ, dft_sep_state.Plasticⁿ⁺¹)

            # Accumulate stress/strain VM data if storage is provided
            if stress_strain_storage !== nothing
                stress_outer, strain_outer = compute_vm_superelastic(dft_sep_state.superelasticⁿ.sigma, dft_sep_state.superelasticⁿ.eps)
                stress_inner, strain_inner = compute_vm_plastic(dft_sep_state.Plasticⁿ.sigma, dft_sep_state.Plasticⁿ.epsilon)
                stress_VM = vcat(stress_outer, stress_inner)
                strain_VM = vcat(strain_outer, strain_inner)
                accumulate_stress_strain_data!(stress_strain_storage, state, t, stress_VM, strain_VM)
            end
        elseif beam.material == :DFT_EP
            dft_ep_state = ms[beam.ind]::DFT_EP_States
            copy_plastic_state!(dft_ep_state.Plasticⁿ, dft_ep_state.Plasticⁿ⁺¹)

            # Accumulate stress/strain VM data if storage is provided
            if stress_strain_storage !== nothing
                stress_VM, strain_VM = compute_vm_plastic(dft_ep_state.Plasticⁿ.sigma, dft_ep_state.Plasticⁿ.epsilon)
                accumulate_stress_strain_data!(stress_strain_storage, state, t, stress_VM, strain_VM)
            end
        elseif beam.material == :DFT_PSE
            dft_pse_state = ms[beam.ind]::DFT_PSE_States
            copy_plastic_state!(dft_pse_state.Plasticⁿ, dft_pse_state.Plasticⁿ⁺¹)
            copy_superelastic_state!(dft_pse_state.superelasticⁿ, dft_pse_state.superelasticⁿ⁺¹)

            # Accumulate stress/strain VM data if storage is provided
            if stress_strain_storage !== nothing
                stress_outer, strain_outer = compute_vm_plastic(dft_pse_state.Plasticⁿ.sigma, dft_pse_state.Plasticⁿ.epsilon)
                stress_inner, strain_inner = compute_vm_superelastic(dft_pse_state.superelasticⁿ.sigma, dft_pse_state.superelasticⁿ.eps)
                stress_VM = vcat(stress_outer, stress_inner)
                strain_VM = vcat(strain_outer, strain_inner)
                accumulate_stress_strain_data!(stress_strain_storage, state, t, stress_VM, strain_VM)
            end
        elseif beam.material == :DFT_PP
            dft_pp_state = ms[beam.ind]::DFT_PP_States
            copy_plastic_state!(dft_pp_state.Plasticⁿ_outer, dft_pp_state.Plasticⁿ⁺¹_outer)
            copy_plastic_state!(dft_pp_state.Plasticⁿ_inner, dft_pp_state.Plasticⁿ⁺¹_inner)

            # Accumulate stress/strain VM data if storage is provided
            if stress_strain_storage !== nothing
                stress_outer, strain_outer = compute_vm_plastic(dft_pp_state.Plasticⁿ_outer.sigma, dft_pp_state.Plasticⁿ_outer.epsilon)
                stress_inner, strain_inner = compute_vm_plastic(dft_pp_state.Plasticⁿ_inner.sigma, dft_pp_state.Plasticⁿ_inner.epsilon)
                stress_VM = vcat(stress_outer, stress_inner)
                strain_VM = vcat(strain_outer, strain_inner)
                accumulate_stress_strain_data!(stress_strain_storage, state, t, stress_VM, strain_VM)
            end
        elseif beam.material == :DFT_PE
            dft_pe_state = ms[beam.ind]::DFT_PE_States
            copy_plastic_state!(dft_pe_state.Plasticⁿ, dft_pe_state.Plasticⁿ⁺¹)

            # Accumulate stress/strain VM data if storage is provided
            if stress_strain_storage !== nothing
                stress_VM, strain_VM = compute_vm_plastic(dft_pe_state.Plasticⁿ.sigma, dft_pe_state.Plasticⁿ.epsilon)
                accumulate_stress_strain_data!(stress_strain_storage, state, t, stress_VM, strain_VM)
            end
        end
    end

    # Update beam-to-beam contacts if they are active
    if state.beam2beam_contactsⁿ⁺¹ !== nothing
        state.beam2beam_contactsⁿ.contacts = deepcopy(state.beam2beam_contactsⁿ⁺¹.contacts)
    end
    
end

# -------------------------------------------------------------------
# Reverts nodal values and forces for the next time step 
# if convergence is not achieved.
# -------------------------------------------------------------------

function update_not_converged!(conf::BeamsConfiguration, state::SimulationState)
    
    @inbounds for i in eachindex(conf.nodes)
        # Revert nodal displacements, velocities, and accelerations
        conf.nodes.u[i] = conf.nodes.uⁿ[i]
        conf.nodes.u̇[i] = conf.nodes.u̇ⁿ[i]
        conf.nodes.ü[i] = conf.nodes.üⁿ[i]
        conf.nodes.w[i] = conf.nodes.wⁿ[i]
        conf.nodes.ẇ[i] = conf.nodes.ẇⁿ[i]
        conf.nodes.ẅ[i] = conf.nodes.ẅⁿ[i]
        conf.nodes.R[i] = conf.nodes.Rⁿ[i]
        conf.nodes.ΔR[i] = conf.nodes.ΔRⁿ[i]
    end
    
    # Revert force vectors
    state.forcesⁿ⁺¹.Tⁱⁿᵗ .= state.forcesⁿ.Tⁱⁿᵗ
    state.forcesⁿ⁺¹.Tᵏ .= state.forcesⁿ.Tᵏ
    state.forcesⁿ⁺¹.Tᶜ .= state.forcesⁿ.Tᶜ
    state.forcesⁿ⁺¹.fᵉˣᵗ .= state.forcesⁿ.fᵉˣᵗ
    state.matricesⁿ⁺¹.K .= state.matricesⁿ.K
    state.matricesⁿ⁺¹.C .= state.matricesⁿ.C    
    state.matricesⁿ⁺¹.M .= state.matricesⁿ.M

    # Revert constitutive states
    ms = state.material_states.by_beam
    for beam in LazyRows(conf.beams)
        if beam.material == :superelastic
            se_state = ms[beam.ind]::SuperelasticStates
            copy_superelastic_state!(se_state.superelasticⁿ⁺¹, se_state.superelasticⁿ)
        elseif beam.material == :plastic
            pl_state = ms[beam.ind]::PlasticStates
            copy_plastic_state!(pl_state.Plasticⁿ⁺¹, pl_state.Plasticⁿ)
        elseif beam.material == :DFT_SEP
            dft_sep_state = ms[beam.ind]::DFT_SEP_States
            copy_superelastic_state!(dft_sep_state.superelasticⁿ⁺¹, dft_sep_state.superelasticⁿ)
            copy_plastic_state!(dft_sep_state.Plasticⁿ⁺¹, dft_sep_state.Plasticⁿ)
        elseif beam.material == :DFT_EP
            dft_ep_state = ms[beam.ind]::DFT_EP_States
            copy_plastic_state!(dft_ep_state.Plasticⁿ⁺¹, dft_ep_state.Plasticⁿ)
        elseif beam.material == :DFT_PSE
            dft_pse_state = ms[beam.ind]::DFT_PSE_States
            copy_plastic_state!(dft_pse_state.Plasticⁿ⁺¹, dft_pse_state.Plasticⁿ)
            copy_superelastic_state!(dft_pse_state.superelasticⁿ⁺¹, dft_pse_state.superelasticⁿ)
        elseif beam.material == :DFT_PP
            dft_pp_state = ms[beam.ind]::DFT_PP_States
            copy_plastic_state!(dft_pp_state.Plasticⁿ⁺¹_outer, dft_pp_state.Plasticⁿ_outer)
            copy_plastic_state!(dft_pp_state.Plasticⁿ⁺¹_inner, dft_pp_state.Plasticⁿ_inner)
        elseif beam.material == :DFT_PE
            dft_pe_state = ms[beam.ind]::DFT_PE_States
            copy_plastic_state!(dft_pe_state.Plasticⁿ⁺¹, dft_pe_state.Plasticⁿ)
        end
    end

    # Revert beam-to-beam contacts if they are active
    if state.beam2beam_contactsⁿ⁺¹ !== nothing
        state.beam2beam_contactsⁿ⁺¹.contacts = deepcopy(state.beam2beam_contactsⁿ.contacts)
    end
    
end