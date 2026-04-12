#--------------------------------------------------------------------------------------------------
# FUNCTIONS TO COMPUTE SPARSITY PATTERN FOR BEAM ELEMENTS
#------------------------------------------------------------------------------------------------------

# Function to compute sparsity pattern for beam elements
function sparsity_beams(nodes, beams, fixed_dofs)

    beam_ndofs = 12  # Number of dofs per beam element
    N = length(beams) * beam_ndofs^2  # Total non-zero entries in sparsity pattern for beams

    # # Add constraints to the sparsity pattern if any
    # if !isnothing(constraints)
    #     N += length(constraints) * beam_ndofs^2
    # end

    # Initialize sparsity pattern
    I = Vector{Int}(undef, N)    # Row indices for the sparse matrix
    J = Vector{Int}(undef, N)    # Column indices for the sparse matrix
    infos = Vector{Vector{Vec4{Int}}}(undef, N)  # Store information about non-zero entries

    k = 1  # Index for sparsity pattern
    # Loop through all beams to calculate the sparsity pattern
    for bi in eachindex(beams)
        n1, n2 = beams.node1[bi], beams.node2[bi]  # Get node indices for the beam
        dofs = vcat(nodes.local_dofs[n1], nodes.local_dofs[n2])  # Get dofs for both nodes
        local_dofs = 1  # Local DOF index for the beam element

        # Loop through all combinations of dofs for the beam
        for j in dofs, i in dofs
            I[k] = i
            J[k] = j
            infos[k] = [Vec4(1, bi, local_dofs, i in fixed_dofs || j in fixed_dofs)]  # Store info (type=1 for beam)
            k += 1
            local_dofs += 1
        end
    end

    # Create the sparse matrix using the computed sparsity pattern
    K = sparse(I, J, infos, maximum(I), maximum(J), vcat)

    # Initialize sparsity free vector to store free dofs
    sparsity_free = Vector{Int}(undef, length(K.nzval))
    ksf = 0  # Counter for free dofs

    # Create maps for beams and constraints sparsity
    beams_spmap = [MVector{144, Int}(undef) for _ in 1:length(beams)]
    # constraints_spmap = isnothing(constraints) ? [] : [MVector{144, Int}(undef) for _ in 1:length(constraints)]

    # Loop through non-zero entries in the sparse matrix
    for (i, infos) in enumerate(K.nzval)
        for (type, idx, local_dofs) in infos
            if type == 1
                beams_spmap[idx][local_dofs] = i  # For beams
            end
        end
        # Identify free dofs
        if infos[1][4] == 0
            ksf += 1
            sparsity_free[ksf] = i
        end
    end

    # Resize sparsity_free vector to remove unused entries
    resize!(sparsity_free, ksf)

    # Assign the sparsity map to each beam and constraint (if present)
    max_dofs = 0
    for i in eachindex(beams)
        beams.local_sparsity_map[i] = beams_spmap[i]
        beams.global_sparsity_map[i] = beams_spmap[i]
        if maximum(beams.global_sparsity_map[i])>max_dofs
            max_dofs = maximum(beams.global_sparsity_map[i])
        end
    end

    return I, J, sparsity_free, max_dofs
end

# Function to construct sparse matrices for beam configuration
function sparse_matrices_beams!(conf::BeamsConfiguration)

    @unpack beams, nodes, bcs, ndofs = conf

    free_dofs = bcs.free_dofs  # Free dofs from boundary conditions
    fixed_dofs = bcs.fixed_dofs  # Fixed dofs from boundary conditions
    nfreedofs = length(free_dofs)  # Number of free dofs

    # Compute sparsity pattern for beams and constraints (if present)
    I, J, sparsity_free = sparsity_beams(nodes, beams, fixed_dofs)

    # Create the tangent stiffness matrix (initially empty)
    Ktan = sparse(I, J, 0.0)

    # Extract the sub-matrix for free dofs
    Ktan_free = Ktan[free_dofs, free_dofs]

    # Create matrices structure for further use
    matrices = Matrices(I, J, sparsity_free)

    # Create nodal solution to hold the stiffness matrix and other info
    sol = Solution(Ktan, Ktan_free, ndofs, nfreedofs)

    return matrices, sol
end

# --------------------------------------------------------------------------------------------------
# SPARSITY PATTERN FOR BEAM-TO-BEAM CONTACT
# --------------------------------------------------------------------------------------------------

# Stores one 24×24 sparsity map per ordered beam pair (b1, b2).
struct Beam2BeamSparsityMaps
    maps::Dict{Tuple{Int,Int}, MVector{576, Int}}
end

# Build the combined sparsity pattern that includes:
#   - all single-beam entries         (type 1, 12×12 = 144 per beam)
#   - all canonical beam pair entries  (type 2, 24×24 = 576 per pair, b1 < b2 only)
function sparsity_beams_b2b(nodes, beams, fixed_dofs)

    beam_ndofs    = 12
    contact_ndofs = 24
    n_beams       = length(beams)
    n_pairs       = n_beams * (n_beams - 1) ÷ 2      # canonical pairs: b1 < b2 only

    N_beams   = n_beams * beam_ndofs^2                # 144 per beam
    N_contact = n_pairs * contact_ndofs^2             # 576 per canonical pair
    N         = N_beams + N_contact

    I     = Vector{Int}(undef, N)
    J     = Vector{Int}(undef, N)
    infos = Vector{Vector{Vec4{Int}}}(undef, N)

    k = 1

    # --- Beam loop (type = 1, unchanged) ---
    for bi in eachindex(beams)
        n1, n2 = beams.node1[bi], beams.node2[bi]
        dofs = vcat(nodes.local_dofs[n1], nodes.local_dofs[n2])
        local_dofs = 1
        for j in dofs, i in dofs
            I[k] = i
            J[k] = j
            infos[k] = [Vec4(1, bi, local_dofs, i in fixed_dofs || j in fixed_dofs)]
            k += 1
            local_dofs += 1
        end
    end

    # --- Contact pair loop (type = 2, b1 < b2 only) ---
    # DOF ordering: [local_dofs of b2 (12), local_dofs of b1 (12)]
    # This matches vcat(dof₂, dof₁) used in the assembly.
    ci = 1
    pair_ci = Dict{Tuple{Int,Int}, Int}()   # (b1, b2) with b1 < b2 → contact index ci

    for b1 in 1:n_beams, b2 in (b1+1):n_beams
        nb1_1, nb1_2 = beams.node1[b1], beams.node2[b1]   # nodes of b1 (smaller index)
        nb2_1, nb2_2 = beams.node1[b2], beams.node2[b2]   # nodes of b2 (larger index)
        # DOF order: b2 first (= dof₂ in assembly), b1 second (= dof₁ in assembly)
        dofs = vcat(nodes.local_dofs[nb2_1], nodes.local_dofs[nb2_2],
                    nodes.local_dofs[nb1_1], nodes.local_dofs[nb1_2])
        local_dofs = 1
        for j in dofs, i in dofs
            I[k] = i
            J[k] = j
            infos[k] = [Vec4(2, ci, local_dofs, i in fixed_dofs || j in fixed_dofs)]
            k += 1
            local_dofs += 1
        end
        pair_ci[(b1, b2)] = ci
        ci += 1
    end

    # Build combined sparse matrix
    K = sparse(I, J, infos, maximum(I), maximum(J), vcat)

    sparsity_free  = Vector{Int}(undef, length(K.nzval))
    ksf = 0

    beams_spmap   = [MVector{144, Int}(undef) for _ in 1:n_beams]
    contact_spmap = [MVector{576, Int}(undef) for _ in 1:n_pairs]

    for (i, infos_i) in enumerate(K.nzval)
        for info in infos_i
            type, idx, local_dofs = info[1], info[2], info[3]
            if type == 1
                beams_spmap[idx][local_dofs] = i
            elseif type == 2
                contact_spmap[idx][local_dofs] = i
            end
        end
        if infos_i[1][4] == 0
            ksf += 1
            sparsity_free[ksf] = i
        end
    end

    resize!(sparsity_free, ksf)

    # Store beam sparsity maps on beam structs
    max_dofs = 0
    for i in eachindex(beams)
        beams.local_sparsity_map[i] = beams_spmap[i]
        beams.global_sparsity_map[i] = beams_spmap[i]
        if maximum(beams.global_sparsity_map[i]) > max_dofs
            max_dofs = maximum(beams.global_sparsity_map[i])
        end
    end

    # Build Beam2BeamSparsityMaps
    b2b_maps = Dict{Tuple{Int,Int}, MVector{576, Int}}()
    for (pair, ci_val) in pair_ci
        b2b_maps[pair] = contact_spmap[ci_val]
    end
    b2b_sparsity = Beam2BeamSparsityMaps(b2b_maps)

    return I, J, sparsity_free, max_dofs, b2b_sparsity
end

# Build sparse matrices including beam-to-beam contact sparsity.
function sparse_matrices_beams_b2b!(conf::BeamsConfiguration)

    @unpack beams, nodes, bcs, ndofs = conf

    free_dofs  = bcs.free_dofs
    fixed_dofs = bcs.fixed_dofs
    nfreedofs  = length(free_dofs)

    I, J, sparsity_free, _, b2b_sparsity = sparsity_beams_b2b(nodes, beams, fixed_dofs)

    Ktan      = sparse(I, J, 0.0)
    Ktan_free = Ktan[free_dofs, free_dofs]

    matrices = Matrices(I, J, sparsity_free)
    sol      = Solution(Ktan, Ktan_free, ndofs, nfreedofs)

    return matrices, sol, b2b_sparsity
end

