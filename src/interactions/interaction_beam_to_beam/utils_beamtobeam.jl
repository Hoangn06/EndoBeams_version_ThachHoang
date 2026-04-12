#------------------------------------------------
# CONTACT DETECTION FUNCTION TO SEARCH CANDIDATES
#------------------------------------------------

function beams2beam_search_canditates(conf)

    @unpack nodes, beams = conf

    nbeams = length(beams)
    nnodes = length(nodes)
    dist_elements = 1E6*ones(Int, nbeams, nbeams)
    dist_elements_ini = 1E6*ones(Int, nbeams)
    dist_comp = 1E6*ones(Int, nbeams, nbeams)
    fact_safe = 1.1
    TNOneighbour = ones(Int, nbeams, nbeams)
    Xlist = zeros(Int, nnodes, 2)

    for (idx, b) in enumerate(beams)

        n₁ = b.node1
        n₂ = b.node2
        
        if Xlist[n₁, 1] == 0
            Xlist[n₁, 1] = idx
        else 
            Xlist[n₁, 2] = idx
        end 
        if Xlist[n₂, 1] == 0
            Xlist[n₂, 1] = idx
        else 
            Xlist[n₂, 2] = idx
        end 

    end 

    ind = findall(!iszero, Xlist[:,2])
 
    for i = eachindex(ind)
        e_neigh = Xlist[ind[i],:]
        beamₐ = e_neigh[1]
        beamᵦ = e_neigh[2]
        TNOneighbour[beamₐ,beamᵦ] = 0
        TNOneighbour[beamᵦ,beamₐ] = 0
    end

    for i = 1:nbeams
        TNOneighbour[i,i] = 0
    end

    for beamₐ in beams
        
        n₁ = beamₐ.node1
        n₂ = beamₐ.node2
        X₁, X₂ = nodes.X₀[n₁], nodes.X₀[n₂]
        u₁, u₂ = nodes.u[n₁], nodes.u[n₂]
        x₁ =  X₁ + u₁
        x₂ =  X₂ + u₂
        lₙb₁ = norm(x₂ - x₁)

        dist_elements_ini[beamₐ.ind] =  beamₐ.l₀

        xb₁ = 0.5*(x₁+x₂)

        for beamᵦ in beams[beamₐ.ind+1:end]
            
            n₁ = beamᵦ.node1
            n₂ = beamᵦ.node2
            X₁, X₂ = nodes.X₀[n₁], nodes.X₀[n₂]
            u₁, u₂ = nodes.u[n₁], nodes.u[n₂]
            x₁ =  X₁ + u₁
            x₂ =  X₂ + u₂
            lₙb₂ = norm(x₂ - x₁)
        
            xb₂ = 0.5*(x₁+x₂)
    
            dist_b₁_b₂= norm(xb₁-xb₂)
            dist_elements[beamₐ.ind,beamᵦ.ind] = dist_b₁_b₂
            dist_elements[beamᵦ.ind,beamₐ.ind] = dist_b₁_b₂
            
            dist_comp[beamₐ.ind,beamᵦ.ind] = 0.5*fact_safe*(lₙb₁+lₙb₂)
            dist_comp[beamᵦ.ind,beamₐ.ind] = 0.5*fact_safe*(lₙb₁+lₙb₂)

        end 

    end 

    check_2 = dist_elements .< dist_comp
    check_final = TNOneighbour .* check_2

    aux = findall(x->x==1, check_final) 

    candidate_elements = Vec2[]
    for i in eachindex(aux)
        push!(candidate_elements, sort([aux[i][1],aux[i][2]]))
    end 

    candidate_elements = unique(candidate_elements, dims=1)

    return candidate_elements

end 


#------------------------------------------------
# FUNCTION TO GET BEAM INFORMATION
#------------------------------------------------

function get_beam_info(beam, nodes)
     
    Nₐ = beam.node1
    Nᵦ = beam.node2

    # information from node 1 and 2
    X₁, X₂ = nodes.X₀[Nₐ], nodes.X₀[Nᵦ]
    u₁, u₂ = nodes.u[Nₐ], nodes.u[Nᵦ]
    R₁, R₂ = nodes.R[Nₐ], nodes.R[Nᵦ]

    init = (X₁, X₂, beam.l₀, beam.Rₑ⁰)      
    beaminfo = (u₁, u₂, R₁, R₂, init)

    return beaminfo

end 


#------------------------------------------------
# FUNCTION TO COMPUTE RₑH₁Rₑᵀ AND ITS DERIVATIVES
#------------------------------------------------
function beam2beam_get_RₑH₁Rₑᵀ_and_its_derivatives(ξ, u₁, u₂, R₁, R₂, init)

    X₁, X₂, l₀, Rₑ⁰ = init
    
    x₁ =  X₁ + u₁
    x₂ =  X₂ + u₂
    
    lₙ = norm(x₂-x₁)
    
    Rₑ, _, _, _, Gᵀ¹, Gᵀ², Gᵀ³, Gᵀ⁴, _ = local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰[:,2], lₙ)
    R̅₁ = Rₑ' * R₁ * Rₑ⁰
    R̅₂ = Rₑ' * R₂ * Rₑ⁰


    Θ̅₁ = toangle(R̅₁)
    Θ̅₂ = toangle(R̅₂)
    
    # Shape functions and its derivatives and second derivatives
    N₁ = 1-ξ/l₀
    N₂ = 1-N₁
    N₃ = ξ*(1-ξ/l₀)^2
    N₄ = -(1-ξ/l₀)*((ξ^2)/l₀)

    dN₁ = -1/l₀
    dN₂ = 1/l₀
    dN₃ = (1-ξ/l₀) * (1 - 3*ξ/l₀)
    dN₄ = (ξ/l₀) * (3*ξ/l₀ - 2)

    ddN₁ = 0
    ddN₂ = 0
    ddN₃ = (-1/l₀)*(1 - 3*ξ/l₀) + (1 - ξ/l₀) * (-3/l₀)
    ddN₄ = 1/l₀ * (3*ξ/l₀ - 2) + (ξ/l₀) * (3/l₀)

    # Displacement vector and its derivatives and second derivatives
    uᵗ = @SVector [0, N₃*Θ̅₁[3] + N₄*Θ̅₂[3], -N₃*Θ̅₁[2] + -N₄*Θ̅₂[2]]
    duᵗ = @SVector [0, dN₃*Θ̅₁[3] + dN₄*Θ̅₂[3], -dN₃*Θ̅₁[2] + -dN₄*Θ̅₂[2]]
    dduᵗ = @SVector [0, ddN₃*Θ̅₁[3] + ddN₄*Θ̅₂[3], -ddN₃*Θ̅₁[2] + -ddN₄*Θ̅₂[2]]

    Suᵗ = skew(uᵗ)
    Sduᵗ = skew(duᵗ)
    Sdduᵗ = skew(dduᵗ)
   
    # Compute terms related to RₑH₁Rₑᵀ
    P₁P¹ = @SMatrix [0 0 0; 0 (N₃+N₄)/lₙ 0;0 0 (N₃+N₄)/lₙ]
    P₁P² = @SMatrix [0 0 0; 0 0 N₃;0 -N₃ 0]
    P₁P³ = -P₁P¹
    P₁P⁴ = @SMatrix [0 0 0; 0 0 N₄;0 -N₄ 0]

    H₁¹ = N₁*ID3 + P₁P¹ - Suᵗ*Gᵀ¹
    H₁² =          P₁P² - Suᵗ*Gᵀ²
    H₁³ = N₂*ID3 + P₁P³ - Suᵗ*Gᵀ³
    H₁⁴ =          P₁P⁴ - Suᵗ*Gᵀ⁴
    
    # Compute terms related to RₑdH₁Rₑᵀ
    dP₁P¹ = @SMatrix [0 0 0; 0 (dN₃+dN₄)/lₙ 0;0 0 (dN₃+dN₄)/lₙ]
    dP₁P² = @SMatrix [0 0 0; 0 0 dN₃;0 -dN₃ 0]
    dP₁P³ = -dP₁P¹
    dP₁P⁴ = @SMatrix [0 0 0; 0 0 dN₄;0 -dN₄ 0]

    dH₁¹ = dN₁*ID3 + dP₁P¹ - Sduᵗ*Gᵀ¹
    dH₁² =           dP₁P² - Sduᵗ*Gᵀ²
    dH₁³ = dN₂*ID3 + dP₁P³ - Sduᵗ*Gᵀ³
    dH₁⁴ =           dP₁P⁴ - Sduᵗ*Gᵀ⁴

    # Compute terms related to RₑddH₁Rₑᵀ
    ddP₁P¹ = @SMatrix [0 0 0; 0 (ddN₃+ddN₄)/lₙ 0;0 0 (ddN₃+ddN₄)/lₙ]
    ddP₁P² = @SMatrix [0 0 0; 0 0 ddN₃;0 -ddN₃ 0]
    ddP₁P³ = -ddP₁P¹
    ddP₁P⁴ = @SMatrix [0 0 0; 0 0 ddN₄;0 -ddN₄ 0]

    ddH₁¹ = ddN₁*ID3 + ddP₁P¹ - Sdduᵗ*Gᵀ¹
    ddH₁² =           ddP₁P² - Sdduᵗ*Gᵀ²
    ddH₁³ = ddN₂*ID3 + ddP₁P³ - Sdduᵗ*Gᵀ³
    ddH₁⁴ =           ddP₁P⁴ - Sdduᵗ*Gᵀ⁴
    
    # Compute the RₑH₁Rₑᵀ components
    RₑH₁¹Rₑᵀ = Rₑ * H₁¹ * Rₑ'
    RₑH₁²Rₑᵀ = Rₑ * H₁² * Rₑ'
    RₑH₁³Rₑᵀ = Rₑ * H₁³ * Rₑ'
    RₑH₁⁴Rₑᵀ = Rₑ * H₁⁴ * Rₑ'
    
    # Compute the RₑdH₁Rₑᵀ components
    RₑdH₁¹Rₑᵀ = Rₑ * dH₁¹ * Rₑ'
    RₑdH₁²Rₑᵀ = Rₑ * dH₁² * Rₑ'
    RₑdH₁³Rₑᵀ = Rₑ * dH₁³ * Rₑ'
    RₑdH₁⁴Rₑᵀ = Rₑ * dH₁⁴ * Rₑ'

    # Compute the RₑddH₁Rₑᵀ components
    RₑddH₁¹Rₑᵀ = Rₑ * ddH₁¹ * Rₑ'
    RₑddH₁²Rₑᵀ = Rₑ * ddH₁² * Rₑ'
    RₑddH₁³Rₑᵀ = Rₑ * ddH₁³ * Rₑ'
    RₑddH₁⁴Rₑᵀ = Rₑ * ddH₁⁴ * Rₑ'

    return RₑH₁¹Rₑᵀ, RₑH₁²Rₑᵀ, RₑH₁³Rₑᵀ, RₑH₁⁴Rₑᵀ, RₑdH₁¹Rₑᵀ, RₑdH₁²Rₑᵀ, RₑdH₁³Rₑᵀ, RₑdH₁⁴Rₑᵀ, RₑddH₁¹Rₑᵀ, RₑddH₁²Rₑᵀ, RₑddH₁³Rₑᵀ, RₑddH₁⁴Rₑᵀ

end 

#------------------------------------------------
# FUNCTION TO COMPUTE DERIVATIVES OF THE BEAM POSITION
#------------------------------------------------
function beam2beam_get_derivatives(ξ, u₁, u₂, R₁, R₂, init)

    X₁, X₂, l₀, Rₑ⁰ = init
    
    x₁ =  X₁ + u₁
    x₂ =  X₂ + u₂
    
    lₙ = norm(x₂ - x₁)
    Rₑ, _, _, _, _, _, _, _, _ = local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰[:,2], lₙ)

    R̅₁ = Rₑ' * R₁ * Rₑ⁰
    R̅₂ = Rₑ' * R₂ * Rₑ⁰

    Θ̅₁ = toangle(R̅₁)
    Θ̅₂ = toangle(R̅₂)
    
    # Shape functions and its derivatives 
    N₁ = 1-ξ/l₀
    N₂ = 1-N₁
    N₃ = ξ*(1-ξ/l₀)^2
    N₄ = -(1-ξ/l₀)*((ξ^2)/l₀)

    dN₃ = (1-ξ/l₀) * (1 - 3*ξ/l₀)
    dN₄ = (ξ/l₀) * (3*ξ/l₀ - 2)

    ddN₃ = (-1/l₀)*(1 - 3*ξ/l₀) + (1 - ξ/l₀) * (-3/l₀)
    ddN₄ = 1/l₀ * (3*ξ/l₀ - 2) + (ξ/l₀) * (3/l₀)

    dddN₃ = 6/(l₀^2)
    dddN₄ = 6/(l₀^2) 
    
    # Displacement vector and its derivatives 
    uᵗ = @SVector [0, N₃*Θ̅₁[3] + N₄*Θ̅₂[3], -N₃*Θ̅₁[2] + -N₄*Θ̅₂[2]]
    duᵗ = @SVector [0, dN₃*Θ̅₁[3] + dN₄*Θ̅₂[3], -dN₃*Θ̅₁[2] + -dN₄*Θ̅₂[2]]
    dduᵗ = @SVector [0, ddN₃*Θ̅₁[3] + ddN₄*Θ̅₂[3], -ddN₃*Θ̅₁[2] + -ddN₄*Θ̅₂[2]]
    ddduᵗ = @SVector [0, dddN₃*Θ̅₁[3] + dddN₄*Θ̅₂[3], -dddN₃*Θ̅₁[2] + -dddN₄*Θ̅₂[2]]
    
    # Compute the position and its derivatives 
    xᴳ = N₁*x₁ + N₂*x₂ + Rₑ*uᵗ
    dxᴳ =  (1/l₀)*(x₂-x₁) + Rₑ*duᵗ
    ddxᴳ = Rₑ*dduᵗ
    dddxᴳ = Rₑ*ddduᵗ

    return xᴳ, dxᴳ, ddxᴳ, dddxᴳ

end 


#------------------------------------------------
# FUNCTION TO COMPUTE DISTANCE BETWEEN BEAMS
#------------------------------------------------

function beam2beam_compute_dist_candidate(beaminfo₁, beaminfo₂)

    ξᶜ₁ = 0
    ξᶜ₂ = 0

    (u₁₁, u₂₁, R₁₁, R₂₁, init₁) = beaminfo₁
    (u₁₂, u₂₂, R₁₂, R₂₂, init₂) = beaminfo₂

    _, _, l₀₁, _ = init₁
    _, _, l₀₂, _ = init₂

    xᶜ₁, dxᶜ₁, ddxᶜ₁, _  = beam2beam_get_derivatives(ξᶜ₁, u₁₁, u₂₁, R₁₁, R₂₁, init₁)
    xᶜ₂, dxᶜ₂, ddxᶜ₂, _  = beam2beam_get_derivatives(ξᶜ₂, u₁₂, u₂₂, R₁₂, R₂₂, init₂)

    A = Mat22(dot(dxᶜ₂, dxᶜ₂), dot(dxᶜ₂,dxᶜ₁), -dot(dxᶜ₂,dxᶜ₁), -dot(dxᶜ₁,dxᶜ₁))
    b = Vec2(-dot((xᶜ₂-xᶜ₁),dxᶜ₂), -dot((xᶜ₂-xᶜ₁),dxᶜ₁))

    ξᶜ = A\b

    ξᶜ₂ = ξᶜ[1]
    ξᶜ₁ = ξᶜ[2]

    tol = 1E-6
    norm_res = 1E6
    max_it = 50
    k = 0

    while norm_res>tol && k<max_it

        xᶜ₁, dxᶜ₁, ddxᶜ₁, _  = beam2beam_get_derivatives(ξᶜ₁, u₁₁, u₂₁, R₁₁, R₂₁, init₁)
        xᶜ₂, dxᶜ₂, ddxᶜ₂, _  = beam2beam_get_derivatives(ξᶜ₂, u₁₂, u₂₂, R₁₂, R₂₂, init₂)
        
        A = Mat22(dot(dxᶜ₂, dxᶜ₂)+dot((xᶜ₂-xᶜ₁),ddxᶜ₂), dot(dxᶜ₂,dxᶜ₁), -dot(dxᶜ₂,dxᶜ₁), -dot(dxᶜ₁,dxᶜ₁) + dot((xᶜ₂-xᶜ₁),ddxᶜ₁))
        b = Vec2(-dot((xᶜ₂-xᶜ₁),dxᶜ₂), -dot((xᶜ₂-xᶜ₁),dxᶜ₁))
                
        Δξᶜ = A\b

        Δξᶜ₂ = Δξᶜ[1]
        Δξᶜ₁ = Δξᶜ[2]

        ξᶜ₁ = ξᶜ₁ + Δξᶜ₁
        ξᶜ₂ = ξᶜ₂ + Δξᶜ₂
        
        norm_res = norm(b)
        k = k+1

    end 
    
    # Recompute the position at the candidate point
    xᶜ₁, _, _, _  = beam2beam_get_derivatives(ξᶜ₁, u₁₁, u₂₁, R₁₁, R₂₁, init₁)
    xᶜ₂, _, _, _  = beam2beam_get_derivatives(ξᶜ₂, u₁₂, u₂₂, R₁₂, R₂₂, init₂)

    # Newton distance local checking phase
    check_1 = ξᶜ₁  >= 0 && ξᶜ₁  <= l₀₁
    check_2 = ξᶜ₂ >=0 && ξᶜ₂ <= l₀₂
    
    # Check If the contact points are out of bounds
    if check_1==1 && check_2 == 1
        dᵇˡ = norm(xᶜ₂-xᶜ₁)
    else
        dᵇˡ = 1E9
    end
    
    # Check if the Newton local converged
    if norm_res>tol
        dᵇˡ = 1E9
    end

    return dᵇˡ, ξᶜ₁, ξᶜ₂

end 

#------------------------------------------------
# FUNCTION TO COMPUTE THE POSITION AND MATRIX G AT THE PREVIOUS CONTACT POINT
#------------------------------------------------
function beam2beam_get_position_and_matrix_G_at_previous_contact_point(state, conf)
    @unpack beams, nodes = conf

    # Loop over all active contacts in the previous timestep
    for (key, contact) in state.beam2beam_contactsⁿ.contacts
        # Extract beam indices from the key
        beam1_ind, beam2_ind = key
        
        # Get contact parameters
        ξᶜ₁_n = contact.ξᶜ₁
        ξᶜ₂_n = contact.ξᶜ₂
        
        # Get beam information
        beam₁ = beams[beam1_ind]
        beam₂ = beams[beam2_ind]
        beaminfo₁ = get_beam_info(beam₁, nodes)
        beaminfo₂ = get_beam_info(beam₂, nodes)
        
        # Compute the beam information at the current time step
        u₁₁, u₂₁, R₁₁, R₂₁, init₁ = beaminfo₁
        u₁₂, u₂₂, R₁₂, R₂₂, init₂ = beaminfo₂

        # Compute the derivatives of the beam positions with the previous contact point
        xᵖ₁, _, _, _  = beam2beam_get_derivatives(ξᶜ₁_n, u₁₁, u₂₁, R₁₁, R₂₁, init₁)
        xᵖ₂, _, _, _  = beam2beam_get_derivatives(ξᶜ₂_n, u₁₂, u₂₂, R₁₂, R₂₂, init₂)

        # Compute the matrix G at the previous contact point
        # Compute the RₑH₁Rₑᵀ and RₑdH₁Rₑᵀ matrices
        RₑH₁¹Rₑᵀ₁, RₑH₁²Rₑᵀ₁, RₑH₁³Rₑᵀ₁, RₑH₁⁴Rₑᵀ₁, _, _, _, _, _, _, _, _ = beam2beam_get_RₑH₁Rₑᵀ_and_its_derivatives(ξᶜ₁_n, u₁₁, u₂₁, R₁₁, R₂₁, init₁)
        RₑH₁¹Rₑᵀ₂, RₑH₁²Rₑᵀ₂, RₑH₁³Rₑᵀ₂, RₑH₁⁴Rₑᵀ₂, _, _, _, _, _, _, _, _ = beam2beam_get_RₑH₁Rₑᵀ_and_its_derivatives(ξᶜ₂_n, u₁₂, u₂₂, R₁₂, R₂₂, init₂)

        G1P = hcat(RₑH₁¹Rₑᵀ₁, RₑH₁²Rₑᵀ₁, RₑH₁³Rₑᵀ₁, RₑH₁⁴Rₑᵀ₁)
        G2P = hcat(RₑH₁¹Rₑᵀ₂, RₑH₁²Rₑᵀ₂, RₑH₁³Rₑᵀ₂, RₑH₁⁴Rₑᵀ₂)

        # Update the contact structure with computed positions and matrices
        updated_contact = BeamPairInContact(beam1_ind, beam2_ind, 
                                           ξᶜ₁_n, ξᶜ₂_n, 
                                           xᵖ₁, xᵖ₂, 
                                           G1P, G2P, 
                                           contact.gᵀ₁, contact.gᵀ₂, 
                                           contact.gᵀᴾ₁, contact.gᵀᴾ₂)

        Update_BeamPair!(state.beam2beam_contactsⁿ, updated_contact)
        
    end
    
end

function regularize_gap_penalty_beams2beams(gap::Float64)
  
    # If the gap is non-positive, use a simple quadratic penalty formulation
    if gap ≤ 0
        regularized_penalty = gap  # Regularized penalty
    else
        regularized_penalty = 0
    end
    
    return regularized_penalty  # Return the penalty
end