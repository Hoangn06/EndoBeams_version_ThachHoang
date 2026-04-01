function plasticity_return_mapping(E::Float64, G::Float64, ν::Float64, y0::Float64, K::Float64,
    plastic_state::Union{PlasticStates, DFT_SEP_States, DFT_EP_States, DFT_PSE_States, DFT_PE_States},
    ind_G::Int, x1::Float64, x2::Float64, x3::Float64, l0::Float64,
    ū::Float64, Θ̅₁::AbstractVector{Float64}, Θ̅₂::AbstractVector{Float64}, beam_ind::Int=0)

    # Step 1: Initialize variables from converged state
    epsilon_pl_n = plastic_state.Plasticⁿ.epsilon_pl[:,ind_G]
    e_pl_n = plastic_state.Plasticⁿ.e_pl[ind_G]

    # Step 2: Compute total strain at the current Gauss point
    dofs = @SVector [ū, Θ̅₁[1], Θ̅₁[2], Θ̅₁[3], Θ̅₂[1], Θ̅₂[2], Θ̅₂[3]]
    Be = matrix_Be(x1, x2, x3, l0)
    epsilon_total = Be * dofs

    # Step 3: Elastic predictor
    C = matrix_constitutive(E, G)
    epsilon_e_trial = epsilon_total - epsilon_pl_n
    sigma_trial = C * epsilon_e_trial

    P = @SMatrix [1.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 3.0]
    h_trial = sqrt(sigma_trial' * P * sigma_trial)
    y_trial = y0 + K * e_pl_n
    F_trial = h_trial - y_trial

    # Step 4: Check yield condition
    if F_trial <= 0.0
        # Elastic step — store trial state and return
        plastic_state.Plasticⁿ⁺¹.sigma[:,ind_G]      = sigma_trial
        plastic_state.Plasticⁿ⁺¹.epsilon[:,ind_G]    = epsilon_total
        plastic_state.Plasticⁿ⁺¹.epsilon_pl[:,ind_G] = epsilon_pl_n
        plastic_state.Plasticⁿ⁺¹.e_pl[ind_G]         = e_pl_n
        plastic_state.Plasticⁿ⁺¹.F[ind_G]            = F_trial

        return sigma_trial, C, true
    else
        # Plastic step — local Newton return mapping
        γ = 0.0
        sigma  = sigma_trial
        e_pl   = e_pl_n
        C_bar  = C
        h      = h_trial
        y      = y_trial
        F      = F_trial

        tol    = 1e-9
        max_it = 100
        k      = 0

        while abs(F) > tol && k < max_it
            N     = (1.0 / h) * P * sigma
            dF_dγ = -((h/y) * (1 - (γ/y)*K) * N' * C_bar * N + K)
            γ     = γ - F / dF_dγ

            e_pl  = e_pl_n + γ
            y     = y0 + K * e_pl

            Cbar1 = E / (1 + (E * γ)/y)
            Cbar2 = G / (1 + (3*G * γ)/y)
            Cbar3 = G / (1 + (3*G * γ)/y)
            C_bar = @SMatrix [Cbar1 0.0 0.0; 0.0 Cbar2 0.0; 0.0 0.0 Cbar3]

            sigma = C_bar * epsilon_e_trial
            h     = sqrt(sigma' * P * sigma)
            F     = h - y
            k    += 1
        end

        if abs(F) > tol
            # Keep return types consistent on failure to avoid shape errors upstream
            return zeros(Vec3{Float64}), zeros(Mat33{Float64}), false
        end

        N          = (1.0 / h) * P * sigma
        epsilon_pl = epsilon_pl_n + γ * N

        # Update plastic state
        plastic_state.Plasticⁿ⁺¹.sigma[:,ind_G]      = sigma
        plastic_state.Plasticⁿ⁺¹.epsilon[:,ind_G]    = epsilon_total
        plastic_state.Plasticⁿ⁺¹.epsilon_pl[:,ind_G] = epsilon_pl
        plastic_state.Plasticⁿ⁺¹.e_pl[ind_G]         = e_pl
        plastic_state.Plasticⁿ⁺¹.F[ind_G]            = F

        beta    = K / (1 - γ * K / y)
        C_const = C_bar - (C_bar * N * N' * C_bar) / (N' * C_bar * N + beta)

        return sigma, C_const, true
    end
end

function plasticity_return_mapping(E::Float64, G::Float64, ν::Float64, y0::Float64, K::Float64,
    plastic_state::DFT_PP_States, layer::Symbol,
    ind_G::Int, x1::Float64, x2::Float64, x3::Float64, l0::Float64,
    ū::Float64, Θ̅₁::AbstractVector{Float64}, Θ̅₂::AbstractVector{Float64}, beam_ind::Int=0)

    previous, current = if layer === :outer
        (plastic_state.Plasticⁿ_outer, plastic_state.Plasticⁿ⁺¹_outer)
    elseif layer === :inner
        (plastic_state.Plasticⁿ_inner, plastic_state.Plasticⁿ⁺¹_inner)
    else
        error("Invalid DFT_PP layer: $(layer). Use :outer or :inner.")
    end

    # Step 1: Initialize variables from converged state
    epsilon_pl_n = previous.epsilon_pl[:,ind_G]
    e_pl_n = previous.e_pl[ind_G]

    # Step 2: Compute total strain at the current Gauss point
    dofs = @SVector [ū, Θ̅₁[1], Θ̅₁[2], Θ̅₁[3], Θ̅₂[1], Θ̅₂[2], Θ̅₂[3]]
    Be = matrix_Be(x1, x2, x3, l0)
    epsilon_total = Be * dofs

    # Step 3: Elastic predictor
    C = matrix_constitutive(E, G)
    epsilon_e_trial = epsilon_total - epsilon_pl_n
    sigma_trial = C * epsilon_e_trial

    P = @SMatrix [1.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 3.0]
    h_trial = sqrt(sigma_trial' * P * sigma_trial)
    y_trial = y0 + K * e_pl_n
    F_trial = h_trial - y_trial

    # Step 4: Check yield condition
    if F_trial <= 0.0
        # Elastic step — store trial state and return
        current.sigma[:,ind_G]      = sigma_trial
        current.epsilon[:,ind_G]    = epsilon_total
        current.epsilon_pl[:,ind_G] = epsilon_pl_n
        current.e_pl[ind_G]         = e_pl_n
        current.F[ind_G]            = F_trial

        return sigma_trial, C, true
    else
        # Plastic step — local Newton return mapping
        γ = 0.0
        sigma  = sigma_trial
        e_pl   = e_pl_n
        C_bar  = C
        h      = h_trial
        y      = y_trial
        F      = F_trial

        tol    = 1e-9
        max_it = 100
        k      = 0

        while abs(F) > tol && k < max_it
            N     = (1.0 / h) * P * sigma
            dF_dγ = -((h/y) * (1 - (γ/y)*K) * N' * C_bar * N + K)
            γ     = γ - F / dF_dγ

            e_pl  = e_pl_n + γ
            y     = y0 + K * e_pl

            Cbar1 = E / (1 + (E * γ)/y)
            Cbar2 = G / (1 + (3*G * γ)/y)
            Cbar3 = G / (1 + (3*G * γ)/y)
            C_bar = @SMatrix [Cbar1 0.0 0.0; 0.0 Cbar2 0.0; 0.0 0.0 Cbar3]

            sigma = C_bar * epsilon_e_trial
            h     = sqrt(sigma' * P * sigma)
            F     = h - y
            k    += 1
        end

        if abs(F) > tol
            # Keep return types consistent on failure to avoid shape errors upstream
            return zeros(Vec3{Float64}), zeros(Mat33{Float64}), false
        end

        N          = (1.0 / h) * P * sigma
        epsilon_pl = epsilon_pl_n + γ * N

        # Update plastic state
        current.sigma[:,ind_G]      = sigma
        current.epsilon[:,ind_G]    = epsilon_total
        current.epsilon_pl[:,ind_G] = epsilon_pl
        current.e_pl[ind_G]         = e_pl
        current.F[ind_G]            = F

        beta    = K / (1 - γ * K / y)
        C_const = C_bar - (C_bar * N * N' * C_bar) / (N' * C_bar * N + beta)

        return sigma, C_const, true
    end
end
