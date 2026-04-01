function compute_DFT_PSE_contributions( 
        K̄ⁱⁿᵗ_inner, ū, Θ̅₁, Θ̅₂, B¹, B¹¹, B²¹, B¹², B²², B¹⁴, B²⁴, exact, Tₛ⁻¹Θ̅₁, Tₛ⁻¹Θ̅₂,
        P¹¹, P²¹, P¹², P²², P¹⁴, P²⁴, η, lₙ, D₃, Gᵀ¹, Gᵀ², Gᵀ⁴, Rₑ, r¹,
        B̄⁺¹¹, B̄⁺²¹, B̄⁺¹², B̄⁺²², B̄⁺¹⁴, B̄⁺²⁴,
        radius, radius_inner,
        E_outer, ν_outer, G_outer, y0, K,
        Eᴬ, Eᴹ, ν_inner, Gᴬ, Gᴹ, εL, α, R_s_AS, R_s_SA, R_f_AS, R_f_SA,
        l₀, nᴳ, zᴳ, ωᴳ, nˢ, qˢ, rˢ, tˢ, wˢ,
        zθ, wθ, zr, wr,
        DFT_PSE_state, beam_ind
    
        )
    # Initialize force and stiffness contributions
    strain_energy = 0.0
    activation = false
    local_converged = true  # Track convergence of local Newton iterations

    zero_row = SMatrix{1,3,Float64}(0.0, 0.0, 0.0)

    # ============================================================
    # Compute the plastic contributions of the outer layer
    # ============================================================

    # Initialize outer layer contributions
    T̄ⁱⁿᵗ_outer = zeros(Vec7{Float64})
    K̄ⁱⁿᵗ_outer = zeros(Mat77{Float64})

    id_G = 1
    for i in 1:nᴳ
       x1 = (l₀/2)*(zᴳ[i]+1)
       for j in eachindex(zr)
            r = zr[j]
            for k in eachindex(zθ)
                θ = zθ[k]
                x2 = r * cos(θ)
                x3 = r * sin(θ)

                sigma, C_const, converged = plasticity_return_mapping(E_outer, G_outer, ν_outer, y0, K,
                    DFT_PSE_state,
                    id_G, x1, x2, x3, l₀,
                    ū, Θ̅₁, Θ̅₂, beam_ind)

                if !converged
                    local_converged = false
                end

                Be   = matrix_Be(x1, x2, x3, l₀)
                pref = (l₀/2)
                wtot = ωᴳ[i]*wr[j]*wθ[k] * r

                T̄ⁱⁿᵗ_outer += (pref*wtot) * (Be' * sigma)
                K̄ⁱⁿᵗ_outer += (pref*wtot) * (Be' * C_const * Be)

                id_G += 1
            end
        end
    end

    T̄ⁱⁿᵗū_outer = T̄ⁱⁿᵗ_outer[1]
    T̄ⁱⁿᵗΘ̅₁_outer = Vec3(T̄ⁱⁿᵗ_outer[2], T̄ⁱⁿᵗ_outer[3], T̄ⁱⁿᵗ_outer[4])
    T̄ⁱⁿᵗΘ̅₂_outer = Vec3(T̄ⁱⁿᵗ_outer[5], T̄ⁱⁿᵗ_outer[6], T̄ⁱⁿᵗ_outer[7])


    # ============================================================
    # Compute the superelastic contributions of the inner layer
    # ============================================================
    # Loop over Gauss points to check if there is tranformation at any Gauss point
    id_G = 1
    for i in 1:nᴳ
       x1 = (l₀/2)*(zᴳ[i]+1)
       for j in 1:nˢ
            for k in 1:nˢ
                x2 = radius_inner * qˢ[j] * rˢ[k]
                x3 = radius_inner * tˢ[j] * rˢ[k]
                # Compute trial state at the current Gauss point
                trial_state(Eᴬ, Eᴹ, Gᴬ, Gᴹ, εL, α, R_s_AS, R_s_SA,
                   DFT_PSE_state,
                   id_G, x1, x2, x3, l₀,
                   ū, Θ̅₁, Θ̅₂, beam_ind)
                
                # Change Gauss ID
                id_G +=1
            end
        end
    end

    # Check if there is any activation ( H_AS = 1 or H_SA = 1 or H_SS = 1)
    activation = any(DFT_PSE_state.superelasticⁿ⁺¹.H_AS .== 1.0) || any(DFT_PSE_state.superelasticⁿ⁺¹.H_SA .== 1.0) || any(DFT_PSE_state.superelasticⁿ⁺¹.H_SS .== 1.0)


    # If there is activation or if xi_S currently not equal to 0
    # Compute the superlastic internal force and stiffness matrix
    if activation
        # Initialize the local vector and matrix for the superelasticity computataion
        T̄ⁱⁿᵗ_inner = zeros(Vec7{Float64})
        K̄ⁱⁿᵗ_inner = zeros(Mat77{Float64})
        id_G = 1
        for i in 1:nᴳ
            x1 = (l₀/2)*(zᴳ[i]+1)
            for j in 1:nˢ
                for k in 1:nˢ
                    x2 = radius_inner * qˢ[j] * rˢ[k]
                    x3 = radius_inner * tˢ[j] * rˢ[k]

                    # Compute the Newton local at each Gauss point to find Martensite fraction
                    converged = superelasticity_newton_local(Eᴬ, Eᴹ, Gᴬ, Gᴹ,
                        εL, α, 
                        R_f_AS, R_f_SA, 
                        DFT_PSE_state,
                        id_G, beam_ind)
                    
                    # If local Newton did not converge, signal non-convergence
                    if !converged
                        local_converged = false
                    end
                    
                    # Compute Matrix of shape function derivatives for this Gauss point
                    Be = matrix_Be(x1, x2, x3, l₀)

                    # Compute the variables used for Intergration Gauss method
                    pref = (l₀/2)*((pi*radius_inner^2)/nˢ)
                    wtot = ωᴳ[i]*wˢ[k]

                    # Compute the internal force in this Gauss point 
                    T̄ⁱⁿᵗ_gauss = superelasticity_internal_force(DFT_PSE_state, id_G, Be)
                    T̄ⁱⁿᵗ_inner += (pref*wtot)*T̄ⁱⁿᵗ_gauss

                    # Compute the stiffness matrix in this Gauss point
                    K̄ⁱⁿᵗ_gauss = superelasticity_matrix_tangent(Eᴬ, Eᴹ, Gᴬ, Gᴹ, εL, α, R_f_AS, R_f_SA,
                                                    DFT_PSE_state, id_G, Be, beam_ind)

                    K̄ⁱⁿᵗ_inner += (pref*wtot)*K̄ⁱⁿᵗ_gauss

                    # Change Gauss ID
                    id_G +=1
                end
            end
        end

        T̄ⁱⁿᵗū_inner  = T̄ⁱⁿᵗ_inner[1]
        T̄ⁱⁿᵗΘ̅₁_inner = Vec3(T̄ⁱⁿᵗ_inner[2], T̄ⁱⁿᵗ_inner[3], T̄ⁱⁿᵗ_inner[4])
        T̄ⁱⁿᵗΘ̅₂_inner = Vec3(T̄ⁱⁿᵗ_inner[5], T̄ⁱⁿᵗ_inner[6], T̄ⁱⁿᵗ_inner[7])
            
    else
        # If there is no activation, compute the elastic contributions of the outer layer
        K̄ⁱⁿᵗū_inner, K̄ⁱⁿᵗΘ̅_inner, K̄ⁱⁿᵗΘ̅Θ̅_inner = K̄ⁱⁿᵗ_inner

        # Compute elastic forces: T̄ⁱⁿᵗ = K̄ⁱⁿᵗ * D̄
        T̄ⁱⁿᵗū_inner  = K̄ⁱⁿᵗū_inner  * ū
        T̄ⁱⁿᵗΘ̅₁_inner = K̄ⁱⁿᵗΘ̅_inner  * Θ̅₁ + K̄ⁱⁿᵗΘ̅Θ̅_inner * Θ̅₂
        T̄ⁱⁿᵗΘ̅₂_inner = K̄ⁱⁿᵗΘ̅Θ̅_inner * Θ̅₁ + K̄ⁱⁿᵗΘ̅_inner  * Θ̅₂

        K̄ⁱⁿᵗ_inner = vcat(
            hcat(K̄ⁱⁿᵗū_inner, zero_row, zero_row),
            hcat(zero_row', K̄ⁱⁿᵗΘ̅_inner, K̄ⁱⁿᵗΘ̅Θ̅_inner),
            hcat(zero_row', K̄ⁱⁿᵗΘ̅Θ̅_inner, K̄ⁱⁿᵗΘ̅_inner)
        )
    end


    # ============================================================
    # Compute the global internal contribution
    # ============================================================
    
    # Assemble the full tangent stiffness matrix in the local coordinates
    K̄ⁱⁿᵗ = K̄ⁱⁿᵗ_outer + K̄ⁱⁿᵗ_inner

    # Transform from local coordinates (T̄ⁱⁿᵗ, K̄ⁱⁿᵗ) to global coordinates (Tⁱⁿᵗ, Kⁱⁿᵗ)
    # Extract components from T̄ⁱⁿᵗ (Vec7): [ū, Θ̅₁[1], Θ̅₁[2], Θ̅₁[3], Θ̅₂[1], Θ̅₂[2], Θ̅₂[3]] +  the inner layer contributions
    T̄ⁱⁿᵗū  = T̄ⁱⁿᵗū_outer + T̄ⁱⁿᵗū_inner
    T̄ⁱⁿᵗΘ̅₁ = T̄ⁱⁿᵗΘ̅₁_outer + T̄ⁱⁿᵗΘ̅₁_inner
    T̄ⁱⁿᵗΘ̅₂ = T̄ⁱⁿᵗΘ̅₂_outer + T̄ⁱⁿᵗΘ̅₂_inner

    # Transform forces: Tⁱⁿᵗ = Bᵀ T̄ⁱⁿᵗ
    Tⁱⁿᵗ¹ = B¹'*T̄ⁱⁿᵗū + B¹¹'*T̄ⁱⁿᵗΘ̅₁ + B²¹'*T̄ⁱⁿᵗΘ̅₂
    Tⁱⁿᵗ² =             B¹²'*T̄ⁱⁿᵗΘ̅₁ + B²²'*T̄ⁱⁿᵗΘ̅₂
    Tⁱⁿᵗ³ = -Tⁱⁿᵗ¹
    Tⁱⁿᵗ⁴ =             B¹⁴'*T̄ⁱⁿᵗΘ̅₁ + B²⁴'*T̄ⁱⁿᵗΘ̅₂

    Tⁱⁿᵗ = [Tⁱⁿᵗ¹; 
            Tⁱⁿᵗ²; 
            Tⁱⁿᵗ³; 
            Tⁱⁿᵗ⁴]

    # Compute [N̄ M̄⁺₁ M̄⁺₂] = B̄ᵀ T̄ⁱⁿᵗ (for geometric stiffness)
    N̄   = T̄ⁱⁿᵗū
    M̄⁺₁ = exact ? Tₛ⁻¹Θ̅₁' * T̄ⁱⁿᵗΘ̅₁  : T̄ⁱⁿᵗΘ̅₁
    M̄⁺₂ = exact ? Tₛ⁻¹Θ̅₂' * T̄ⁱⁿᵗΘ̅₂  : T̄ⁱⁿᵗΘ̅₂

    # Qₛ = Pᵀ [M̄⁺₁ M̄⁺₂]
    Qₛ¹ = P¹¹' * M̄⁺₁ + P²¹' * M̄⁺₂
    Qₛ² = P¹²' * M̄⁺₁ + P²²' * M̄⁺₂
    Qₛ⁴ = P¹⁴' * M̄⁺₁ + P²⁴' * M̄⁺₂
        
    # Q = S(Qₛ)
    Q¹ = skew(Qₛ¹)
    Q² = skew(Qₛ²)
    Q⁴ = skew(Qₛ⁴)

    a = @SVector [0, η*(M̄⁺₁[1] + M̄⁺₂[1])/lₙ + (M̄⁺₁[2] + M̄⁺₂[2])/lₙ, (M̄⁺₁[3] + M̄⁺₂[3])/lₙ]

    # DN̄ (DN̄¹¹ = DN̄³³ = -DN̄¹³ = -DN̄³¹)
    DN̄¹¹ = D₃*N̄

    # QGᵀ
    QGᵀ¹¹ = Q¹*Gᵀ¹
    QGᵀ¹² = Q¹*Gᵀ²
    QGᵀ¹⁴ = Q¹*Gᵀ⁴
    QGᵀ²² = Q²*Gᵀ²
    QGᵀ²⁴ = Q²*Gᵀ⁴
    QGᵀ⁴⁴ = Q⁴*Gᵀ⁴

    # EGa
    EGa¹ = Rₑ*Gᵀ¹'*a

    # EGar
    EGar¹¹ = EGa¹*r¹

    # Kₘ = DN̄ - EQGᵀEᵀ + EGar (geometric stiffness)
    Kₘ¹¹ = DN̄¹¹ - Rₑ*QGᵀ¹¹*Rₑ' + EGar¹¹
    Kₘ¹² =      - Rₑ*QGᵀ¹²*Rₑ'
    Kₘ¹⁴ =      - Rₑ*QGᵀ¹⁴*Rₑ'
    Kₘ²² =      - Rₑ*QGᵀ²²*Rₑ'
    Kₘ²⁴ =      - Rₑ*QGᵀ²⁴*Rₑ'
    Kₘ⁴⁴ =      - Rₑ*QGᵀ⁴⁴*Rₑ'

    # K̃ (geometric stiffness with exact rotation)
    if exact
            η₁, μ₁ = compute_η_μ(Θ̅₁)
            η₂, μ₂ = compute_η_μ(Θ̅₂)

            M̄₁ = T̄ⁱⁿᵗΘ̅₁
            M̄₂ = T̄ⁱⁿᵗΘ̅₂

            K̄ₕ₁ = compute_K̄ₕ(Θ̅₁, M̄₁, Tₛ⁻¹Θ̅₁, η₁, μ₁)
            K̄ₕ₂ = compute_K̄ₕ(Θ̅₂, M̄₂, Tₛ⁻¹Θ̅₂, η₂, μ₂)
    end

    K̃¹¹ = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹¹  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²¹  +  Kₘ¹¹   :   Kₘ¹¹ 
    K̃¹² = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹²  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²²  +  Kₘ¹²   :   Kₘ¹² 
    K̃¹⁴ = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ¹⁴   :   Kₘ¹⁴ 
    K̃²² = exact ?  B̄⁺¹²' * K̄ₕ₁ * B̄⁺¹²  +  B̄⁺²²' * K̄ₕ₂ * B̄⁺²²  +  Kₘ²²   :   Kₘ²² 
    K̃²⁴ = exact ?  B̄⁺¹²' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²²' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ²⁴   :   Kₘ²⁴ 
    K̃⁴⁴ = exact ?  B̄⁺¹⁴' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²⁴' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ⁴⁴   :   Kₘ⁴⁴ 


    B = vcat(
            hcat(B¹, zero_row, -B¹, zero_row),              
            hcat(B¹¹, B¹², -B¹¹, B¹⁴),                      
            hcat(B²¹, B²², -B²¹, B²⁴)                       
    )

    # Compute B' * K̄ⁱⁿᵗ * B directly (12×12 matrix)
    BᵀKⁱⁿᵗB = B' * K̄ⁱⁿᵗ * B
        
    # Construct full K̃ matrix (12×12): geometric stiffness matrix
    # Structure: K̃ = [K̃¹¹    K̃¹²    -K̃¹¹    K̃¹⁴  ]
    #                [K̃¹²'    K̃²²    -K̃¹²'    K̃²⁴  ]
    #                [-K̃¹¹'    -K̃¹²    K̃¹¹    -K̃¹⁴  ]
    #                [K̃¹⁴'    K̃²⁴'    -K̃¹⁴'    K̃⁴⁴  ]
    K̃ = hcat(vcat(K̃¹¹, K̃¹²', -K̃¹¹', K̃¹⁴'),
            vcat(K̃¹², K̃²², -K̃¹², K̃²⁴'),
            vcat(-K̃¹¹, -K̃¹²', K̃¹¹, -K̃¹⁴'),
            vcat(K̃¹⁴, K̃²⁴, -K̃¹⁴, K̃⁴⁴))

    # Compute the global stiffness matrix
    Kⁱⁿᵗ = BᵀKⁱⁿᵗB + K̃

    return Tⁱⁿᵗ, Kⁱⁿᵗ, strain_energy, local_converged
end
