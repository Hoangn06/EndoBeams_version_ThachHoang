function beam2beam_compute_Kc_Tc(beaminfo₁, beaminfo₂, ξᶜ₁, ξᶜ₂, gₙ, first_contact, ξᶜ₁_n, ξᶜ₂_n, xᵖ₁, xᵖ₂, G1P, G2P, gᵀ₁_n, gᵀ₂_n, gᵀᴾ₁_n, gᵀᴾ₂_n)

    u₁₁, u₂₁, R₁₁, R₂₁, init₁ = beaminfo₁
    u₁₂, u₂₂, R₁₂, R₂₂, init₂ = beaminfo₂
    
    # Compute the derivatives of the beam positions
    xᶜ₁, dxᶜ₁, ddxᶜ₁, dddxᶜ₁  = beam2beam_get_derivatives(ξᶜ₁, u₁₁, u₂₁, R₁₁, R₂₁, init₁)
    xᶜ₂, dxᶜ₂, ddxᶜ₂, dddxᶜ₂  = beam2beam_get_derivatives(ξᶜ₂, u₁₂, u₂₂, R₁₂, R₂₂, init₂)
    
    # Compute the RₑH₁Rₑᵀ and RₑdH₁Rₑᵀ matrices
    RₑH₁¹Rₑᵀ₁, RₑH₁²Rₑᵀ₁, RₑH₁³Rₑᵀ₁, RₑH₁⁴Rₑᵀ₁, RₑdH₁¹Rₑᵀ₁, RₑdH₁²Rₑᵀ₁, RₑdH₁³Rₑᵀ₁, RₑdH₁⁴Rₑᵀ₁, RₑddH₁¹Rₑᵀ₁, RₑddH₁²Rₑᵀ₁, RₑddH₁³Rₑᵀ₁, RₑddH₁⁴Rₑᵀ₁ = beam2beam_get_RₑH₁Rₑᵀ_and_its_derivatives(ξᶜ₁, u₁₁, u₂₁, R₁₁, R₂₁, init₁)
    RₑH₁¹Rₑᵀ₂, RₑH₁²Rₑᵀ₂, RₑH₁³Rₑᵀ₂, RₑH₁⁴Rₑᵀ₂, RₑdH₁¹Rₑᵀ₂, RₑdH₁²Rₑᵀ₂, RₑdH₁³Rₑᵀ₂, RₑdH₁⁴Rₑᵀ₂, RₑddH₁¹Rₑᵀ₂, RₑddH₁²Rₑᵀ₂, RₑddH₁³Rₑᵀ₂, RₑddH₁⁴Rₑᵀ₂ = beam2beam_get_RₑH₁Rₑᵀ_and_its_derivatives(ξᶜ₂, u₁₂, u₂₂, R₁₂, R₂₂, init₂)

    G1 = hcat(RₑH₁¹Rₑᵀ₁, RₑH₁²Rₑᵀ₁, RₑH₁³Rₑᵀ₁, RₑH₁⁴Rₑᵀ₁)
    G2 = hcat(RₑH₁¹Rₑᵀ₂, RₑH₁²Rₑᵀ₂, RₑH₁³Rₑᵀ₂, RₑH₁⁴Rₑᵀ₂)

    H1 = hcat(RₑdH₁¹Rₑᵀ₁, RₑdH₁²Rₑᵀ₁, RₑdH₁³Rₑᵀ₁, RₑdH₁⁴Rₑᵀ₁)
    H2 = hcat(RₑdH₁¹Rₑᵀ₂, RₑdH₁²Rₑᵀ₂, RₑdH₁³Rₑᵀ₂, RₑdH₁⁴Rₑᵀ₂)

    M1 = hcat(RₑddH₁¹Rₑᵀ₁, RₑddH₁²Rₑᵀ₁, RₑddH₁³Rₑᵀ₁, RₑddH₁⁴Rₑᵀ₁)
    M2 = hcat(RₑddH₁¹Rₑᵀ₂, RₑddH₁²Rₑᵀ₂, RₑddH₁³Rₑᵀ₂, RₑddH₁⁴Rₑᵀ₂)

    # Compute the normal contact vector
    xᶜ₁₂ = xᶜ₂-xᶜ₁
    dᵇˡ = norm(xᶜ₁₂)
    n =  xᶜ₁₂/dᵇˡ
    
    # Compute the G matrix 
    G = hcat(G2, -G1)
    
    # Compute the normal nodal contact force
    εₙ = 500 # Normal stiffness coefficient
    Tᶜ = εₙ * gₙ * G' * n
    
    # Compute the D1 and D2 matrices
    A = Mat22(dot(dxᶜ₂, dxᶜ₂)+dot((xᶜ₂-xᶜ₁),ddxᶜ₂), dot(dxᶜ₂,dxᶜ₁), -dot(dxᶜ₂,dxᶜ₁), -dot(dxᶜ₁,dxᶜ₁) + dot((xᶜ₂-xᶜ₁),ddxᶜ₁))    
    aux1 = - dxᶜ₂'*G2 - (xᶜ₂-xᶜ₁)'*H2
    aux2 = dxᶜ₂'*G1
    aux3 = -dxᶜ₁'*G2
    aux4 = -(xᶜ₂-xᶜ₁)'*H1 + dxᶜ₁'*G1
    B = Mat224(aux1[1], aux3[1], aux1[2], aux3[2], aux1[3], aux3[3], aux1[4], aux3[4], aux1[5], aux3[5], aux1[6], aux3[6],aux1[7], aux3[7], aux1[8], aux3[8], aux1[9], aux3[9], aux1[10], aux3[10],  aux1[11], aux3[11], aux1[12], aux3[12], aux2[1], aux4[1], aux2[2], aux4[2], aux2[3], aux4[3],  aux2[4], aux4[4], aux2[5], aux4[5], aux2[6], aux4[6],aux2[7], aux4[7], aux2[8], aux4[8], aux2[9], aux4[9], aux2[10], aux4[10], aux2[11], aux4[11], aux2[12], aux4[12])
    
    aux = A\B

    D2 = Vec24(aux[1,1], aux[1,2], aux[1,3], aux[1,4], aux[1,5], aux[1,6], aux[1,7], aux[1,8], aux[1,9], aux[1,10], aux[1,11], aux[1,12], aux[1,13], aux[1,14], aux[1,15], aux[1,16], aux[1,17], aux[1,18], aux[1,19], aux[1,20],aux[1,21],aux[1,22],aux[1,23],aux[1,24])
    D1 = Vec24(aux[2,1], aux[2,2], aux[2,3], aux[2,4], aux[2,5], aux[2,6], aux[2,7], aux[2,8], aux[2,9], aux[2,10], aux[2,11], aux[2,12], aux[2,13], aux[2,14], aux[2,15], aux[2,16], aux[2,17], aux[2,18], aux[2,19], aux[2,20],aux[2,21],aux[2,22],aux[2,23],aux[2,24])
    
    # Compute the matrix E for tangent stiffness matrix
    E = vcat(H2' * n * D2', -H1'* n * D1')
    
    # Compute the matrix F for tangent stiffness matrix
    F = D2 * ddxᶜ₂' * n * D2' - D1 * ddxᶜ₁' * n * D1'
    
    # Compute the matrix I for tangent stiffness matrix
    I = (1/dᵇˡ) * (G' + D2 * dxᶜ₂' - D1 * dxᶜ₁')*(ID3 - n * n') * (G + dxᶜ₂ * D2' - dxᶜ₁ * D1')

    Kcn = E + E' + F + I
    
    # Compute the tangent stiffness contribution for normal contact
    Kᶜ =  εₙ*(G' * n * n' * G + gₙ * Kcn )

    # Friction contribution to the contact force and matrix
    if first_contact
        Tᶜ = Tᶜ 
        Kᶜ = Kᶜ 
        gᵀᴾ₁ = 0
        gᵀᴾ₂ = 0
        gᵀ₁ = 0
        gᵀ₂ = 0
    else
        # Compute the total tangential gap 
        # Compute the tangential gap for beam 1
        s₁ = sign(ξᶜ₁-ξᶜ₁_n)
        dgᵀ₁ = norm(xᶜ₁-xᵖ₁)
        gᵀ₁ = gᵀ₁_n + s₁ * dgᵀ₁
        # Compute the tangential gap for beam 2
        s₂ = sign(ξᶜ₂-ξᶜ₂_n)
        dgᵀ₂ = norm(xᶜ₂-xᵖ₂)
        gᵀ₂ = gᵀ₂_n + s₂ * dgᵀ₂

        # Choose the penalty parameter for the Friction
        εₜ = 500 # Tangential stiffness coefficient
        μ = 0.3 # Friction coefficient

        # Compute the normal contact force
        Fᴺ = εₙ * gₙ

        #Check for the stick/slip conditions of beams
        # Check the condition of beam 1 / Return map for beam 1
        # Trial Phase
        gᵀᴱ₁_trial = gᵀ₁ - gᵀᴾ₁_n

        Fᵀ₁_trial = εₜ * gᵀᴱ₁_trial

        f₁_trial = abs(Fᵀ₁_trial) - μ * abs(Fᴺ)
        
        # Check Sttick or Slip for beam 1
        if f₁_trial ≤ 0
            # Stick condition for beam 1
            gᵀᴾ₁ = gᵀᴾ₁_n
            Fᵀ₁ = Fᵀ₁_trial
            gᵀᴱ₁ = gᵀᴱ₁_trial
            println("stick beam 1")

            stick_beam1 = true
        else
            # Slip condition for beam 1
            Fᵀ₁ = μ * abs(Fᴺ) * sign(gᵀᴱ₁_trial)
            gᵀᴾ₁ = gᵀᴾ₁_n + (1/εₜ) * (Fᵀ₁_trial - Fᵀ₁)
            gᵀᴱ₁ = gᵀ₁ - gᵀᴾ₁

            println("slip beam 1")
            stick_beam1 = false
        end

        # Check the condition of beam 2 / Return map for beam 2
        # Trial Phase
        gᵀᴱ₂_trial = gᵀ₂ - gᵀᴾ₂_n

        Fᵀ₂_trial = εₜ * gᵀᴱ₂_trial

        f₂_trial = abs(Fᵀ₂_trial) - μ * abs(Fᴺ)
        
        # Check Sttick or Slip for beam 2
        if f₂_trial ≤ 0
            # Stick condition for beam 2
            gᵀᴾ₂ = gᵀᴾ₂_n
            Fᵀ₂ = Fᵀ₂_trial
            gᵀᴱ₂ = gᵀᴱ₂_trial

            println("stick beam 2")

            stick_beam2 = true
        else
            # Slip condition for beam 2
            Fᵀ₂ = μ * abs(Fᴺ) * sign(gᵀᴱ₂_trial)
            gᵀᴾ₂ = gᵀᴾ₂_n + (1/εₜ) * (Fᵀ₂_trial - Fᵀ₂)
            gᵀᴱ₂ = gᵀ₂ - gᵀᴾ₂

            println("slip beam 2")

            stick_beam2 = false
        end
        
        # Compute the matrices/vectors used for the friction force and matrix
        # Compute the tangential vectors
        t₁ = (xᶜ₁-xᵖ₁)/norm(xᶜ₁-xᵖ₁)
        t₂ = (xᶜ₂-xᵖ₂)/norm(xᶜ₂-xᵖ₂)
        
        # Compute the S matrices
        S1 = dxᶜ₁ * D1' + hcat(zeros(SMatrix{3,12,Float64}), G1 - G1P)
        S2 = dxᶜ₂ * D2' + hcat(G2 - G2P, zeros(SMatrix{3,12,Float64}))

        # Compute the Y matrices
        Y11 = vcat(- H2' * dxᶜ₁, - H1' * dxᶜ₂) # Compute Y11

        Y12 = vcat(M2' * xᶜ₁₂ + 2 * H2' * dxᶜ₂ + G2' * ddxᶜ₂, -G1' * ddxᶜ₂) # Compute Y12

        Y21 = vcat(G2' * ddxᶜ₁ , M1' * xᶜ₁₂ - 2 * H1' * dxᶜ₁ - G1' * ddxᶜ₁ ) # Compute Y21

        # Compute p values
        p₁ = dot(dddxᶜ₁ , xᶜ₁₂) - 3 * dot(ddxᶜ₁, dxᶜ₁)
        p₂ = dot(ddxᶜ₁ , dxᶜ₂)
        p₃ = dot(dxᶜ₁ , ddxᶜ₂)
        p₄ = dot(dddxᶜ₂ , xᶜ₁₂) + 3 * dot(ddxᶜ₂, dxᶜ₂)

        # Compute R1
        R11 = -p₂ * D1 * D1' - p₃ * (D2 * D1' + D1 * D2') + p₄ * D2 * D2'

        R12 = D1 * Y11' + D2 * Y12' + Y11 * D1' - Y12 * D2'

        R131 = H2' * G2 + G2' * H2
        R132 = -H2' * G1
        R133 = -G1' * H2
        R134 = zeros(SMatrix{12,12,Float64})
        
        R13 = vcat(hcat(R131,R132), hcat(R133,R134))

        R1 = R11 + R12 + R13

        # Compute R2
        R21 = p₁ * D1 * D1' + p₂ * (D2 * D1' + D1 * D2') + p₃ * D2 * D2'

        R22 = D1 * Y21' - D2 * Y11' + Y21 * D1' - Y11 * D2'

        R231 = zeros(SMatrix{12,12,Float64})
        R232 = G2' * H1
        R233 = H1' * G2
        R234 = -H1' * G1 - G1' * H1
        
        R23 = vcat(hcat(R231, R232), hcat(R233, R234))

        R2 = R21 + R22 + R23

        # Compute a21 a22 a11 a21
        # Compute the inversed matrix of A
        Ainv = inv(A)

        a21 = Ainv[1,1]  # Top-left
        a22 = Ainv[1,2]  # Top-right  
        a11 = Ainv[2,1]  # Bottom-left
        a12 = Ainv[2,2]  # Bottom-right

        # Compute Z matrices
        Z1 = vcat(zeros(SMatrix{12,3,Float64}) , H1') * t₁ * D1'
        
        Z2 = vcat(H2' , zeros(SMatrix{12,3,Float64})) * t₂ * D2'

        # Compute Kct 
        Kct1 = (1/dgᵀ₁) * S1' * (ID3 - t₁ * t₁') * S1 + t₁' * ddxᶜ₁ * D1 * D1' + Z1 + Z1' + t₁' * dxᶜ₁ * (-a11 * R1 - a12 * R2)
        Kct2 = (1/dgᵀ₂) * S2' * (ID3 - t₂ * t₂') * S2 + t₂' * ddxᶜ₂ * D2 * D2' + Z2 + Z2' + t₂' *  dxᶜ₂ * (-a21 * R1 - a22 * R2)

        if stick_beam1
            Tᵀ₁ = εₜ * gᵀᴱ₁ * S1' * t₁
            Kᵀ₁ = εₜ * (S1' * t₁ * t₁' * S1 + gᵀᴱ₁ * Kct1)

            Tᶜ = Tᶜ + Tᵀ₁
            Kᶜ = Kᶜ + Kᵀ₁
        else
            Tᵀ₁ = μ * abs(Fᴺ) * sign(gᵀᴱ₁_trial) * S1' * t₁
            Kᵀ₁ = μ * εₙ *  sign(gᵀᴱ₁_trial) * (sign(gₙ) * S1' * t₁ * n' * G + abs(gₙ) * Kct1 )

            Tᶜ = Tᶜ + Tᵀ₁
            Kᶜ = Kᶜ + Kᵀ₁
        end

        if stick_beam2
            Tᵀ₂ = εₜ * gᵀᴱ₂ * S2' * t₂
            Kᵀ₂ = εₜ * (S2' * t₂ * t₂' * S2 + gᵀᴱ₂ * Kct2)

            Tᶜ = Tᶜ + Tᵀ₂
            Kᶜ = Kᶜ + Kᵀ₂
        else
            Tᵀ₂ = μ * abs(Fᴺ) * sign(gᵀᴱ₂_trial) * S2' * t₂
            Kᵀ₂ = μ * εₙ *  sign(gᵀᴱ₂_trial) * (sign(gₙ) * S2' * t₂ * n' * G + abs(gₙ) * Kct2 )

            Tᶜ = Tᶜ + Tᵀ₂
            Kᶜ = Kᶜ + Kᵀ₂

        end

    end

    return Tᶜ, Kᶜ, gᵀ₁, gᵀ₂, gᵀᴾ₁, gᵀᴾ₂
end


function beam2beam_compute!(candidate_elements, conf, state)

    @unpack nodes, beams = conf

    # Clear all contacts for the current time step before computing new ones
    Clear_All_Contacts!(state.beam2beam_contactsⁿ⁺¹)
    
    # Check each candidate element to see if it is in contact
    # If it is in contact, compute the contact force and matrix
    for idx = eachindex(candidate_elements)

        beam₁ = beams[candidate_elements[idx][1]]
        beam₂ = beams[candidate_elements[idx][2]]

        beaminfo₁ = get_beam_info(beam₁, nodes)
        beaminfo₂ = get_beam_info(beam₂, nodes)

        @timeit_debug "Compute distance candidates" dᵇˡ, ξᶜ₁, ξᶜ₂ =  beam2beam_compute_dist_candidate(beaminfo₁, beaminfo₂)

        if dᵇˡ < 1e9
            gₙ = regularize_gap_penalty_beams2beams(dᵇˡ - (beam₁.properties.radius + beam₂.properties.radius))

            if gₙ !=0
                # Initialize contact parameters for the previous time step

                # Get previous contact parameters if they exist
                if Check_BeamPair(state.beam2beam_contactsⁿ, beam₁.ind, beam₂.ind)
                    prev_contact = Get_BeamPair(state.beam2beam_contactsⁿ, beam₁.ind, beam₂.ind)
                    ξᶜ₁_n = prev_contact.ξᶜ₁
                    ξᶜ₂_n = prev_contact.ξᶜ₂
                    xᵖ₁ = prev_contact.xᵖ₁
                    xᵖ₂ = prev_contact.xᵖ₂
                    G1P = prev_contact.G1P
                    G2P = prev_contact.G2P
                    gᵀ₁_n = prev_contact.gᵀ₁
                    gᵀ₂_n = prev_contact.gᵀ₂
                    gᵀᴾ₁_n = prev_contact.gᵀᴾ₁
                    gᵀᴾ₂_n = prev_contact.gᵀᴾ₂
                    first_contact = false
                else
                    # First time contact
                    ξᶜ₁_n = 0
                    ξᶜ₂_n = 0
                    xᵖ₁ = zeros(Vec3{Float64})
                    xᵖ₂ = zeros(Vec3{Float64})
                    G1P = zeros(Mat312{Float64})
                    G2P = zeros(Mat312{Float64})
                    gᵀ₁_n = 0
                    gᵀ₂_n = 0
                    gᵀᴾ₁_n = 0
                    gᵀᴾ₂_n = 0
                    first_contact = true
                end

                @timeit_debug "Compute force and matrix"  Tᶜ, Kᶜ, gᵀ₁, gᵀ₂, gᵀᴾ₁, gᵀᴾ₂ =  beam2beam_compute_Kc_Tc(beaminfo₁, beaminfo₂, ξᶜ₁, ξᶜ₂, gₙ, 
                                                                    first_contact, ξᶜ₁_n, ξᶜ₂_n,
                                                                    xᵖ₁, xᵖ₂, 
                                                                    G1P, G2P,
                                                                    gᵀ₁_n, gᵀ₂_n,
                                                                    gᵀᴾ₁_n, gᵀᴾ₂_n)

                # Create or update contact structure
                new_contact = BeamPairInContact(beam₁.ind, beam₂.ind, ξᶜ₁, ξᶜ₂, xᵖ₁, xᵖ₂, G1P, G2P, gᵀ₁, gᵀ₂, gᵀᴾ₁, gᵀᴾ₂)
                Update_BeamPair!(state.beam2beam_contactsⁿ⁺¹, new_contact)
                
                @timeit_debug "Assemble" begin
                    dof₁₁ = nodes.global_dofs[beam₁.node1]
                    dof₁₂ = nodes.global_dofs[beam₁.node2]
                    dof₁  = vcat(dof₁₁, dof₁₂)

                    dof₂₁ = nodes.global_dofs[beam₂.node1]
                    dof₂₂ = nodes.global_dofs[beam₂.node2]
                    dof₂  = vcat(dof₂₁, dof₂₂)

                    # Force vector: DOF ordering is [dof₂ (12), dof₁ (12)]
                    contact_dofs = vcat(dof₂, dof₁)
                    state.forcesⁿ⁺¹.Tᶜ[contact_dofs] .-= Tᶜ

                    # Stiffness matrix: use pre-computed sparsity map for fast nzval assembly.
                    # candidate_elements guarantees beam₁.ind < beam₂.ind, so the key is direct.
                    spmap = state.b2b_sparsity.maps[(beam₁.ind, beam₂.ind)]
                    state.matricesⁿ⁺¹.K.nzval[spmap] .+= vec(Kᶜ)
                end
                
            end 
        end 
    end
end 


