#Compute Residual for forward A→S transformation
function residual_linear( F::Float64, F_prev::Float64, 
    F_f_AS::Float64, lambda_AS::Float64, H_AS::Float64,
    F_f_SA::Float64, lambda_SA::Float64, H_SA::Float64,
    xi::Float64)
    R_AS = F_f_AS * lambda_AS + H_AS * (1.0 - xi) * (F - F_prev)
    R_SA = F_f_SA * lambda_SA - H_SA * xi * (F - F_prev)
    
   return R_AS, R_SA
end


#Compute Jacobian matrix for the residuals of the transformations
function jacobian_linear(F::Float64, F_prev::Float64,
    F_f_AS::Float64, lambda_AS::Float64,H_AS::Float64, 
    F_f_SA::Float64, lambda_SA::Float64,H_SA::Float64, 
    n::AbstractVector{Float64}, xi::Float64, Eᴬ::Float64, Eᴹ::Float64, Gᴬ::Float64, Gᴹ::Float64, εL::Float64, α::Float64, 
    eps_total_elastic::AbstractVector{Float64},
    beam_ind::Int=0, ind_G::Int=0)

    # Compute G1
    Cᴱᴬ = matrix_constitutive(Eᴬ, Gᴬ)
    Cᴱᴹ = matrix_constitutive(Eᴹ, Gᴹ)

    # Voigt model
    Cᴱ = (1 - xi) * Cᴱᴬ + xi * Cᴱᴹ
    E = (1 - xi) * Eᴬ + xi * Eᴹ

    G1 = n' * ((Cᴱᴹ - Cᴱᴬ) * eps_total_elastic - Cᴱ * (n + α * I) * εL) + α * ((Eᴹ - Eᴬ) * eps_total_elastic[1] - E * (n[1] + α ) * εL)

    # Reuss Model
    #Cᴱ = inv((1-xi)*inv(Cᴱᴬ) + xi*inv(Cᴱᴹ))
    #E = (Eᴬ)/(1 + (Eᴬ/Eᴹ -1) * xi)

    #dC = -Cᴱ * (inv(Cᴱᴹ) - inv(Cᴱᴬ)) * Cᴱ
    #dE =  (Eᴬ*(Eᴬ/Eᴹ - 1)) / ((1 + (Eᴬ/Eᴹ -1) * xi)^2)
    #G1 = n' * (dC * eps_total_elastic - Cᴱ * (n + α * I) * εL) + α * (dE * eps_total_elastic[1] - E * (n[1] + α ) * εL)

    # Compute the jacobian component
    a = G1 * lambda_AS - H_AS*((F - F_prev) - (1 - xi) * G1) + F_f_AS
    b = G1 * lambda_AS - H_AS*((F - F_prev) - (1 - xi) * G1)
    c = G1 * lambda_SA - H_SA*((F - F_prev) + xi * G1)
    d = G1 * lambda_SA - H_SA*((F - F_prev) + xi * G1) + F_f_SA

    return a, b, c, d
end

# Function trial_state will be called at each Gauss point to compute the trial state
# This function is used to compute the trial state of the stress and view if there is transformation
function trial_state(Eᴬ::Float64, Eᴹ::Float64, Gᴬ::Float64, Gᴹ::Float64,
    εL::Float64, α::Float64, 
    R_s_AS::Float64, R_s_SA::Float64,
    superelastic_state::Union{SuperelasticStates, DFT_SEP_States, DFT_PSE_States},
    ind_G::Int, x1::Float64, x2::Float64, x3::Float64, l0::Float64,
    ū::Float64, Θ̅₁::AbstractVector{Float64}, Θ̅₂::AbstractVector{Float64}, beam_ind::Int=0)

    # Step 0: Initialize the variables for the trial state
    xi_S = superelastic_state.superelasticⁿ.xi_S[ind_G] # The martensite fraction is the same as the previous step
 
    # Step 0.2: Compute the total strain at the current Gauss points
    # Concatenate DOFs: [ū, Θ̅₁[1], Θ̅₁[2], Θ̅₁[3], Θ̅₂[1], Θ̅₂[2], Θ̅₂[3]]
    dofs = @SVector [ū, Θ̅₁[1], Θ̅₁[2], Θ̅₁[3], Θ̅₂[1], Θ̅₂[2], Θ̅₂[3]] #Compute the vector dofs
    Be = matrix_Be(x1, x2, x3, l0)
    eps_total = Be * dofs

    #Step 1: Compute the deviatoric / volumetric parts of total strain
    e = deviatoric_tensor(eps_total)
    e_norm = norm_strain_deviatoric(e)

    # Compute the direction of the transformation/deviatoric strain/deviatoric stress
    # n = t/||t|| = e/||e||
    # n = [n11 2n12 2n13]
    if e_norm < 1e-10
        n = @SVector [0.0, 0.0, 0.0]
    else
        n = e / e_norm
    end

    #Step 2: Compute the trial state of the stress
    # Supposing that the xi_S is the same as the previous step

    # Compute the trial transformation strain
    u_trial = xi_S * (n + α * I) # Compute the trial transformation strain

    # Compute the trial deviatoric / volumetric parts of stress
    Cᵀᴬ = matrix_constitutive_deviatoric(Eᴬ, Gᴬ) # Construct the matrix constitutive the deviatoric part
    Cᵀᴹ = matrix_constitutive_deviatoric(Eᴹ, Gᴹ) # Construct the matrix constitutive the deviatoric part
    # Voigt Model
    Cᵀ = (1 - xi_S) * Cᵀᴬ + xi_S * Cᵀᴹ
    E = (1 - xi_S) * Eᴬ + xi_S * Eᴹ
    # Reuss Model
    #Cᵀ = inv((1-xi_S)*inv(Cᵀᴬ) + xi_S*inv(Cᵀᴹ))
    #E = (Eᴬ)/(1 + (Eᴬ/Eᴹ -1) * xi_S)

    # Compute the trial deviatoric / volumetric parts of stress
    t_trial = Cᵀ * (eps_total - εL * u_trial)
    p_trial = (1/3)* E * (eps_total[1] - εL * u_trial[1])
    sigma_trial = t_trial + p_trial * I

    tnorm_trial = norm_stress_deviatoric(t_trial)

    # Check for superelastic transformation

    # Compute the DP trial function in the current step
    F_trial = tnorm_trial + 3.0 * α * p_trial

    if F_trial > R_s_AS && F_trial > superelastic_state.superelasticⁿ.F[ind_G] && xi_S < 1.0
        H_AS = 1.0
    else
        H_AS = 0.0
    end

    if F_trial < R_s_SA && F_trial < superelastic_state.superelasticⁿ.F[ind_G] && xi_S > 0.0
        H_SA = 1.0
    else
        H_SA = 0.0
    end

    if xi_S > 0.0
        H_SS = 1.0
    else
        H_SS = 0.0
    end

    # Save the variables into the test structures
    superelastic_state.superelasticⁿ⁺¹.xi_S[ind_G] = xi_S
    superelastic_state.superelasticⁿ⁺¹.F[ind_G] = F_trial

    # Save the variables into the current structures
    superelastic_state.superelasticⁿ⁺¹.H_AS[ind_G] = H_AS
    superelastic_state.superelasticⁿ⁺¹.H_SA[ind_G] = H_SA
    superelastic_state.superelasticⁿ⁺¹.H_SS[ind_G] = H_SS
    superelastic_state.superelasticⁿ⁺¹.eps[:,ind_G] = eps_total
    superelastic_state.superelasticⁿ⁺¹.sigma[:,ind_G] = sigma_trial
    
end

function superelasticity_newton_local(Eᴬ::Float64, Eᴹ::Float64, Gᴬ::Float64, Gᴹ::Float64,
    εL::Float64, α::Float64, 
    R_f_AS::Float64, R_f_SA::Float64, 
    superelastic_state::Union{SuperelasticStates, DFT_SEP_States, DFT_PSE_States},
    ind_G::Int, beam_ind::Int=0)


    if superelastic_state.superelasticⁿ⁺¹.H_AS[ind_G] == 1.0 || superelastic_state.superelasticⁿ⁺¹.H_SA[ind_G] == 1.0
        

        # Extract the variables from the structures for Newton local 
        xi_S = superelastic_state.superelasticⁿ⁺¹.xi_S[ind_G]
        F_trial = superelastic_state.superelasticⁿ⁺¹.F[ind_G]
        F_prev = superelastic_state.superelasticⁿ.F[ind_G]
        eps_total = Vec3(superelastic_state.superelasticⁿ⁺¹.eps[:,ind_G])
        sigma_trial = Vec3(superelastic_state.superelasticⁿ⁺¹.sigma[:,ind_G])
        H_AS = superelastic_state.superelasticⁿ⁺¹.H_AS[ind_G]
        H_SA = superelastic_state.superelasticⁿ⁺¹.H_SA[ind_G]

        # Compute the values for Newton local
        # Compute the direction of the transformation
        e = deviatoric_tensor(eps_total)
        e_norm = norm_strain_deviatoric(e)
        if e_norm < 1e-10
            n = @SVector [0.0, 0.0, 0.0]
        else
            n = e / e_norm
        end

        # Compute the normalised transformation strain
        u_trial = xi_S * (n + α * I)

        # Compute the deviatoric / volumetric parts of stress
        p_trial = (1/3) * sigma_trial[1]
        t_trial = deviatoric_tensor(sigma_trial)

        # Compute the trial deviatoric / volumetric parts of stress
        tnorm_trial = norm_stress_deviatoric(t_trial)

        # Compute the trial DP function
        F_trial = tnorm_trial + 3.0 * α * p_trial
    
        # Compute the current matrix constitutive and current total Young's modulus
        Cᵀᴬ = matrix_constitutive_deviatoric(Eᴬ, Gᴬ)
        Cᵀᴹ = matrix_constitutive_deviatoric(Eᴹ, Gᴹ)

        # Compute the current elastic strain
        eps_total_elastic = eps_total - εL * u_trial
        
        # Compute the F_f_AS and F_f_SA at the first iteration
        F_f_AS_trial = F_trial - R_f_AS
        F_f_SA_trial = F_trial - R_f_SA

        #SOLVE THE NONLINEAR SYSTEM OF EQUATIONS
        tol = 1e-4
        norm_res = 1e6
        max_it = 100
        k = 0
        norm_R_0 = 0.0
        force_out = false
        lambda_AS = 0.0
        lambda_SA = 0.0
        
        #SOLVE THE NONLINEAR SYSTEM OF EQUATIONS
        while norm_res > tol && k < max_it

            if k != 0 # Recompute the trial deviatoric / volumetric parts of stress 

                # Update the matrix constitutive and current total Young's modulus
                # Voigt Model
                Cᵀ = (1 - xi_S) * Cᵀᴬ + xi_S * Cᵀᴹ
                E = (1 - xi_S) * Eᴬ + xi_S * Eᴹ

                # Reuss Model
                #Cᵀ = inv((1-xi_S)*inv(Cᵀᴬ) + xi_S*inv(Cᵀᴹ))
                #E = (Eᴬ)/(1 + (Eᴬ/Eᴹ -1) * xi_S)

                # With the updated xi
                u_trial = xi_S * (n + α * I)
                eps_total_elastic = eps_total - εL * u_trial

                t_trial = Cᵀ * eps_total_elastic
                p_trial = (1/3)* E * (eps_total_elastic[1])

                tnorm_trial = norm_stress_deviatoric(t_trial)

                F_trial = tnorm_trial + 3.0*α * p_trial

                F_f_AS_trial = F_trial - R_f_AS
                F_f_SA_trial = F_trial - R_f_SA
            end

            a, b, c, d = jacobian_linear(F_trial, F_prev, 
                F_f_AS_trial, lambda_AS, H_AS,
                F_f_SA_trial, lambda_SA, H_SA, 
                n, xi_S, Eᴬ, Eᴹ, Gᴬ, Gᴹ, εL, α, eps_total_elastic, beam_ind, ind_G)
    
            res1, res2 = residual_linear(F_trial, F_prev,
                F_f_AS_trial, lambda_AS, H_AS,
                F_f_SA_trial, lambda_SA, H_SA,
                xi_S)
            
            dR = Mat22(a, c, b, d)
            R = Vec2(-res1, -res2)

            
            if det(dR) != 0.0
                Δlambda = dR\R
            else
                # Jacobian is not invertible, newton local failed to converge
                return false
            end

            Δlambda_AS = Δlambda[1]
            Δlambda_SA = Δlambda[2]
        
            lambda_AS = lambda_AS + Δlambda_AS
            lambda_SA = lambda_SA + Δlambda_SA

            xi_S = xi_S + Δlambda_AS + Δlambda_SA

            norm_res = norm(R) 

            #if ind_G == 1 && beam_ind == 1
            #    println("--------------------------------")
            #    println("k: ", k)
            #    println("xi_S: ", xi_S)
            #    println("lambda_AS: ", lambda_AS)
            #    println("lambda_SA: ", lambda_SA)
            #end


            if k == 0
                norm_R_0 = norm(R)
            end
        
            if k != 0
                norm_res = norm(R) / norm_R_0
            end

            if xi_S > (1.0)
                xi_S = 1.0
                force_out = true
                break
            end

            if xi_S < 0.0
                xi_S = 0.0
                force_out = true
                break
            end
            
            k = k+1
        end

        # Check the convergence criteria if there is no force_out
        if force_out == false
            if norm_res > tol
                # Residual is not converged, newton local failed to converge
                return false
            end
        end

        # Compute the current matrix constitutive and current total Young's modulus

        Cᵀ = (1 - xi_S) * Cᵀᴬ + xi_S * Cᵀᴹ
        E = (1 - xi_S) * Eᴬ + xi_S * Eᴹ

        # Reuss Model
        #Cᵀ = inv((1-xi_S)*inv(Cᵀᴬ) + xi_S*inv(Cᵀᴹ))
        #E = (Eᴬ)/(1 + (Eᴬ/Eᴹ -1) * xi_S)

        # Update the stress in this state with the new xi_S
        u = xi_S * (n + α * I)
        t = Cᵀ * (eps_total - εL * u)
        p = (1/3)* E * (eps_total[1] - εL * u[1])
        sigma = t + p * I
        F = norm_stress_deviatoric(t) + 3.0*α * p

        # Update the current structures with the new xi_S
        superelastic_state.superelasticⁿ⁺¹.xi_S[ind_G] = xi_S
        superelastic_state.superelasticⁿ⁺¹.F[ind_G] = F
        superelastic_state.superelasticⁿ⁺¹.sigma[:,ind_G] = sigma

        #if ind_G == 1 && beam_ind == 1
        #    println("xi_S: ", xi_S)
        #    println("Cᵀᴹ: ", Cᵀᴹ)
        #    println("Cᵀ: ", Cᵀ)
        #    println("E: ", E)
        #    println("t: ", t)
        #    println("p: ", p)
        #    println("F: ", F)
        #    println("u: ", u)
        #end

    end

    return true 

end

function superelasticity_internal_force(superelastic_state::Union{SuperelasticStates, DFT_SEP_States, DFT_PSE_States},
    ind_G::Int,
    Be::AbstractMatrix{Float64},
    beam_ind::Int=0)

    sigma = Vec3(superelastic_state.superelasticⁿ⁺¹.sigma[:,ind_G])

    #if ind_G == 1 && beam_ind == 1
    #    println("Internal force: -----------")
    #    println("t: ", t)
    #    println("p: ", p)
    #end
    
    T̄ⁱⁿᵗ_gauss = Be' * sigma

    return T̄ⁱⁿᵗ_gauss

end

function superelasticity_matrix_tangent(Eᴬ::Float64, Eᴹ::Float64, Gᴬ::Float64, Gᴹ::Float64,
    εL::Float64, α::Float64, 
    R_f_AS::Float64, R_f_SA::Float64, 
    superelastic_state::Union{SuperelasticStates, DFT_SEP_States, DFT_PSE_States},
    ind_G::Int,
    Be::AbstractMatrix{Float64},
    beam_ind::Int=0)
    
    # Extract the variables from the structures for the matrix tangent computation
    xi_S = superelastic_state.superelasticⁿ⁺¹.xi_S[ind_G]
    xi_prev = superelastic_state.superelasticⁿ.xi_S[ind_G]
    eps_total = Vec3(superelastic_state.superelasticⁿ⁺¹.eps[:,ind_G])
    H_AS = superelastic_state.superelasticⁿ⁺¹.H_AS[ind_G]
    H_SA = superelastic_state.superelasticⁿ⁺¹.H_SA[ind_G]
    F = superelastic_state.superelasticⁿ⁺¹.F[ind_G]
    F_prev = superelastic_state.superelasticⁿ.F[ind_G]

    Cᴱᴬ = matrix_constitutive(Eᴬ, Gᴬ)
    Cᴱᴹ = matrix_constitutive(Eᴹ, Gᴹ)

    # Voigt model for A-M
    Cᴱ = (1 - xi_S) * Cᴱᴬ + xi_S * Cᴱᴹ
    E = (1 - xi_S) * Eᴬ + xi_S * Eᴹ

    # Reuss model for A-M
    #Cᴱ = inv((1-xi_S)*inv(Cᴱᴬ) + xi_S*inv(Cᴱᴹ))
    #E = (Eᴬ)/(1 + (Eᴬ/Eᴹ -1) * xi_S)

    # Compute the values of the direction of the transformation and the transformation strain
    e = deviatoric_tensor(eps_total)
    e_norm = norm_strain_deviatoric(e)
    if e_norm < 1e-10
        n = @SVector [0.0, 0.0, 0.0]
    else
        n = e / e_norm
    end
    u = xi_S * (n + α * I)


    #if ind_G == 1 && beam_ind == 1
    #    println("Matrix tangent: ----------")
    #    println("xi_S: ", xi_S)
    #    println("F: ", F)
    #    println("u: ", tensors_state.Tensorsⁿ⁺¹.u[:,ind_G])
    #end

    # Compute the current total elastic strain
    eps_total_elastic = eps_total - εL * u

    # Compute lambda_AS and lambda_SA
    if H_AS == 1.0
        lambda_AS = xi_S - xi_prev
    else
        lambda_AS = 0.0
    end

    if H_SA == 1.0
        lambda_SA = xi_S - xi_prev
    else
        lambda_SA = 0.0
    end
    
    # Compute the F_f_AS and F_f_SA in the current step
    F_f_AS = F - R_f_AS
    F_f_SA = F - R_f_SA

    # Compute the Jacobian matrix for the linear model in the current step
    a, b, c, d = jacobian_linear(F, F_prev, 
        F_f_AS, lambda_AS, H_AS,
        F_f_SA, lambda_SA, H_SA, 
        n, xi_S, Eᴬ, Eᴹ, Gᴬ, Gᴹ, εL, α, eps_total_elastic, beam_ind, ind_G)

    # Compute the inverse of the Jacobian matrix = [a_inv b_inv; c_inv d_inv]
    det = a * d - c * b
        
    a_inv = d / det
    b_inv = -b / det
    c_inv = -c / det
    d_inv = a / det

    # Compute A_AS = -dR_AS/dF_AS and A_SA = -dR_SA/dF_SA
    
    A_AS = - lambda_AS - H_AS * (1-xi_S)
    A_SA = - lambda_SA + H_SA * xi_S

    # Compute components of the consistent tangent matrix
    # Compute K11 in the term K1 = εL *  xi * dn 
    # Compute N_d = 1/||e|| * (I - n * n')

    n_t =  Diagonal(Vec3((3.0/2.0), (1.0/2.0), (1.0/2.0))) * n

    if e_norm != 0.0
        N_d = 1.0/e_norm * (Mat33{Float64}(ID3) - n * n_t')
    else
        N_d = Mat33{Float64}(ID3)
    end

    K11 = εL * xi_S * N_d * Mat33{Float64}(Idev)

    # Compute K21 K22 K31 K32
        
    # Compute T_1_AS, T_2_AS, T_1_SA, T_2_SA
    T_1_AS = a_inv * A_AS + b_inv * A_SA
    T_2_AS = E * α * (a_inv * A_AS + b_inv * A_SA)
    T_1_SA = c_inv * A_AS + d_inv * A_SA
    T_2_SA = E * α * (c_inv * A_AS + d_inv * A_SA)

    # Compute K21 and K22
    K21 = εL * (T_1_AS + T_1_SA ) * n * n' * Cᴱ
    K22 = εL * (T_2_AS + T_2_SA ) * n * I'

    # Compute K31 and K32
    K31 = εL * α * (T_1_AS + T_1_SA ) * I * n' * Cᴱ
    K32 = εL * α * (T_2_AS + T_2_SA ) * I * I'

    # Compute K41

    K411 = (T_1_AS + T_1_SA) * n' * Cᴱ + (T_2_AS + T_2_SA) * I'
    # Voigt Model
    K41 = (Cᴱᴹ - Cᴱᴬ) * (eps_total_elastic * K411) 
    # Reuss Model
    #dC = -Cᴱ * (inv(Cᴱᴹ) - inv(Cᴱᴬ)) * Cᴱ
    #K41 = dC * eps_total_elastic * K411

    # Compute the matrix consistent tangent matrix
    C_const = K41 + (Cᴱ - Cᴱ * (K11 + K21 + K22 + K31 + K32))

    #if ind_G == 1 && beam_ind == 1
    #    println("--------------------------------")
    #    println("K41: ", K41)
    #    println("K11: ", K11)
    #    println("K21: ", K21)
    #    println("K22: ", K22)
    #    println("K31: ", K31)
    #    println("K32: ", K32)
    #    println("C_const: ", C_const)
    #end

    K̄ⁱⁿᵗ_gauss = Be' * C_const * Be

    return K̄ⁱⁿᵗ_gauss
    
end








