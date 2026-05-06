
# -------------------------------------------------------------------------------------------
# Convert variables
# -------------------------------------------------------------------------------------------

"""
Convert a vector of `Vec3{Float64}` to a matrix of shape `n × 3`.

# Arguments
- `arr::Vector{Vec3{Float64}}`: Vector of 3D points.

# Returns
- `Matrix{Float64}`: Matrix with 3 columns (x, y, z coordinates), where each row 
  corresponds to a `Vec3` in the input vector.
"""
function vec3_array_to_matrix(arr)
    n = length(arr)
    mat = Matrix{Float64}(undef, n, 3)
    for i in 1:n
        v = arr[i]
        mat[i, 1] = v.x
        mat[i, 2] = v.y
        mat[i, 3] = v.z
    end
    return mat
end

"""
Convert a matrix of shape `n × 3` to a vector of `Vec3{Float64}`.

# Arguments
- `mat::Matrix{Float64}`: Matrix with 3 columns (x, y, z coordinates).

# Returns
- `Vector{Vec3{Float64}}`: Each row of the matrix is converted to a `Vec3`.
"""
function matrix_to_vec3_array(mat)
    @assert size(mat, 2) == 3 "Input matrix must have exactly 3 columns."

    return [Vec3(mat[i, 1], mat[i, 2], mat[i, 3]) for i in 1:size(mat, 1)]
end

"""
Convert a matrix of shape `n × 2` to a vector of `Vec2{Int}`.

# Arguments
- `mat::Matrix{Int}`: Matrix with 2 columns.

# Returns
- `Vector{Vec2{Int}}`: Each row of the matrix is converted to a `Vec2`.
"""
function matrix_to_vec2_array(mat)
    @assert size(mat, 2) == 2 "Input matrix must have exactly 2 columns."

    return [Vec2(mat[i, 1], mat[i, 2]) for i in 1:size(mat, 1)]
end

"""
Convert a matrix of shape `n × 9` to a vector of `Mat33{Float64}` matrices.

# Arguments
- `mat::Matrix{Float64}`: Matrix with 9 columns representing 3×3 matrices row-wise.

# Returns
- `Vector{Mat33{Float64}}`: Each row of the matrix converted to a `Mat33`.
"""
function matrix_to_mat33_array(mat::Matrix{Float64})
    @assert size(mat, 2) == 9 "Input matrix must have exactly 9 columns."

    vec = Vector{Mat33{Float64}}(undef, size(mat, 1))

    for i in 1:size(mat, 1)
        # # reshape row vector into 3x3 matrix in column-major order
        # # note: Julia is column-major, so transpose if needed
        # # here, we assume row-major in input, so reshape and transpose
        # row = mat[i, :]
        # mat33 = reshape(row, 3, 3)'  # transpose to get proper orientation
        vec[i] = Mat33(mat[i, :])
    end

    return vec
end


# -------------------------------------------------------------------------------------------
# Write
# -------------------------------------------------------------------------------------------
"""
Write node and beam solution data (`u`, `R`, `Re0`, and non-linear material states) to text files.

Arguments:
- `nodes`: Object containing displacement `u` and rotation matrices `R` (per node).
- `beams`: Object containing reference element rotations `Rₑ⁰` (per beam).
- `num_nodes`: Number of nodes.
- `num_beams`: Number of beams.
- `output_dir`: Directory where files will be written.
- `CLEAN_FOLDER`: If true, clears existing `.vtk`, `.vtu`, `.txt` files in the folder.

# Keyword arguments
- `material_states`: If not `nothing`, superelastic beams append one row of `xi_S` (martensite fraction
  per Gauss section) to `xi_S.txt`. Values are taken from `superelasticⁿ`, i.e. the state after
  `update_converged!` for the last completed time step (not the trial `superelasticⁿ⁺¹` during iteration).
"""
function export_stent_solution_to_txt(nodes, beams, num_nodes, num_beams, output_dir::String, flag_clean_folder::Bool=true;
    material_states::Union{Nothing, BeamMaterialStates}=nothing)
    
    # Clean existing output files if requested
    if flag_clean_folder && output_dir != ""
        savedir = pwd()
        cd(output_dir)
        for ext in [".vtk", ".vtu", ".txt"]
            for file in filter(f -> endswith(f, ext), readdir())
                rm(file, force=true)
            end
        end
        cd(savedir)
    end

    # Ensure directory path ends without double slash
    if !endswith(output_dir, "/")
        output_dir *= "/"
    end

    # Write node displacements and rotations
    for i in 1:num_nodes
        # Write displacement vector u (3 values)
        open(output_dir * "u.txt", "a") do io
            writedlm(io, [nodes.u[i]'])  # Transpose to row
        end

        # Write rotation matrix R (flattened 3x3)
        R = nodes.R[i]
        open(output_dir * "R.txt", "a") do io
            writedlm(io, [vec(R)'])  # Flatten and transpose to row
        end
    end

    # Write beam reference rotations
    for i in 1:num_beams
        Re0 = beams.Rₑ⁰[i]
        open(output_dir * "Re0.txt", "a") do io
            writedlm(io, [vec(Re0)'])  # Flatten 3x3 matrix to row
        end
    end

    # Write data from plastic beams
    if any(beams.material[i] == :plastic for i in 1:num_beams)
        pl_sigma_cols      = Vector{Float64}[]
        pl_epsilon_cols    = Vector{Float64}[]
        pl_epsilon_pl_cols = Vector{Float64}[]
        pl_e_pl_cols       = Vector{Float64}[]
        pl_F_cols          = Vector{Float64}[]
        pl_beam_indices    = Int[]
        for i in 1:num_beams
            if beams.material[i] == :plastic
                pl_state = material_states.by_beam[beams.ind[i]]::PlasticStates
                push!(pl_sigma_cols,      vec(pl_state.Plasticⁿ.sigma))
                push!(pl_epsilon_cols,    vec(pl_state.Plasticⁿ.epsilon))
                push!(pl_epsilon_pl_cols, vec(pl_state.Plasticⁿ.epsilon_pl))
                push!(pl_e_pl_cols,       pl_state.Plasticⁿ.e_pl)
                push!(pl_F_cols,          pl_state.Plasticⁿ.F)
                push!(pl_beam_indices, beams.ind[i])
            end
        end
        open(output_dir * "pl_sigma.txt", "w") do io
            writedlm(io, pl_beam_indices')            # first row: beam indices
            writedlm(io, hcat(pl_sigma_cols...))      # (3*nG*nS) × num_plastic_beams
        end
        open(output_dir * "pl_epsilon.txt", "w") do io
            writedlm(io, pl_beam_indices')            # first row: beam indices
            writedlm(io, hcat(pl_epsilon_cols...))    # (3*nG*nS) × num_plastic_beams
        end
        open(output_dir * "pl_epsilon_pl.txt", "w") do io
            writedlm(io, pl_beam_indices')             # first row: beam indices
            writedlm(io, hcat(pl_epsilon_pl_cols...))  # (3*nG*nS) × num_plastic_beams
        end
        open(output_dir * "pl_e_pl.txt", "w") do io
            writedlm(io, pl_beam_indices')           # first row: beam indices
            writedlm(io, hcat(pl_e_pl_cols...))      # (nG*nS) × num_plastic_beams
        end
        open(output_dir * "pl_F.txt", "w") do io
            writedlm(io, pl_beam_indices')        # first row: beam indices
            writedlm(io, hcat(pl_F_cols...))      # (nG*nS) × num_plastic_beams
        end
    end

    # Write data from DFT_EP beams
    if any(beams.material[i] == :DFT_EP for i in 1:num_beams)
        ep_sigma_cols      = Vector{Float64}[]
        ep_epsilon_cols    = Vector{Float64}[]
        ep_epsilon_pl_cols = Vector{Float64}[]
        ep_e_pl_cols       = Vector{Float64}[]
        ep_F_cols          = Vector{Float64}[]
        ep_beam_indices    = Int[]
        for i in 1:num_beams
            if beams.material[i] == :DFT_EP
                st = material_states.by_beam[beams.ind[i]]::DFT_EP_States
                push!(ep_sigma_cols,      vec(st.Plasticⁿ.sigma))
                push!(ep_epsilon_cols,    vec(st.Plasticⁿ.epsilon))
                push!(ep_epsilon_pl_cols, vec(st.Plasticⁿ.epsilon_pl))
                push!(ep_e_pl_cols,       st.Plasticⁿ.e_pl)
                push!(ep_F_cols,          st.Plasticⁿ.F)
                push!(ep_beam_indices, beams.ind[i])
            end
        end
        open(output_dir * "ep_sigma.txt", "w") do io
            writedlm(io, ep_beam_indices')
            writedlm(io, hcat(ep_sigma_cols...))      # (3*nG*nS) × num_DFT_EP_beams
        end
        open(output_dir * "ep_epsilon.txt", "w") do io
            writedlm(io, ep_beam_indices')
            writedlm(io, hcat(ep_epsilon_cols...))    # (3*nG*nS) × num_DFT_EP_beams
        end
        open(output_dir * "ep_epsilon_pl.txt", "w") do io
            writedlm(io, ep_beam_indices')
            writedlm(io, hcat(ep_epsilon_pl_cols...)) # (3*nG*nS) × num_DFT_EP_beams
        end
        open(output_dir * "ep_e_pl.txt", "w") do io
            writedlm(io, ep_beam_indices')
            writedlm(io, hcat(ep_e_pl_cols...))       # (nG*nS) × num_DFT_EP_beams
        end
        open(output_dir * "ep_F.txt", "w") do io
            writedlm(io, ep_beam_indices')
            writedlm(io, hcat(ep_F_cols...))          # (nG*nS) × num_DFT_EP_beams
        end
    end

    # Write data from DFT_PE beams
    if any(beams.material[i] == :DFT_PE for i in 1:num_beams)
        pe_sigma_cols      = Vector{Float64}[]
        pe_epsilon_cols    = Vector{Float64}[]
        pe_epsilon_pl_cols = Vector{Float64}[]
        pe_e_pl_cols       = Vector{Float64}[]
        pe_F_cols          = Vector{Float64}[]
        pe_beam_indices    = Int[]
        for i in 1:num_beams
            if beams.material[i] == :DFT_PE
                st = material_states.by_beam[beams.ind[i]]::DFT_PE_States
                push!(pe_sigma_cols,      vec(st.Plasticⁿ.sigma))
                push!(pe_epsilon_cols,    vec(st.Plasticⁿ.epsilon))
                push!(pe_epsilon_pl_cols, vec(st.Plasticⁿ.epsilon_pl))
                push!(pe_e_pl_cols,       st.Plasticⁿ.e_pl)
                push!(pe_F_cols,          st.Plasticⁿ.F)
                push!(pe_beam_indices, beams.ind[i])
            end
        end
        open(output_dir * "pe_sigma.txt", "w") do io
            writedlm(io, pe_beam_indices')
            writedlm(io, hcat(pe_sigma_cols...))      # (3*nG*nS) × num_DFT_PE_beams
        end
        open(output_dir * "pe_epsilon.txt", "w") do io
            writedlm(io, pe_beam_indices')
            writedlm(io, hcat(pe_epsilon_cols...))    # (3*nG*nS) × num_DFT_PE_beams
        end
        open(output_dir * "pe_epsilon_pl.txt", "w") do io
            writedlm(io, pe_beam_indices')
            writedlm(io, hcat(pe_epsilon_pl_cols...)) # (3*nG*nS) × num_DFT_PE_beams
        end
        open(output_dir * "pe_e_pl.txt", "w") do io
            writedlm(io, pe_beam_indices')
            writedlm(io, hcat(pe_e_pl_cols...))       # (nG*nS) × num_DFT_PE_beams
        end
        open(output_dir * "pe_F.txt", "w") do io
            writedlm(io, pe_beam_indices')
            writedlm(io, hcat(pe_F_cols...))          # (nG*nS) × num_DFT_PE_beams
        end
    end

    # Write data from DFT_SEP beams
    if any(beams.material[i] == :DFT_SEP for i in 1:num_beams)
        sep_se_xi_S_cols      = Vector{Float64}[]
        sep_se_F_cols         = Vector{Float64}[]
        sep_se_eps_cols       = Vector{Float64}[]
        sep_se_sigma_cols     = Vector{Float64}[]
        sep_pl_sigma_cols     = Vector{Float64}[]
        sep_pl_epsilon_cols   = Vector{Float64}[]
        sep_pl_epsilon_pl_cols= Vector{Float64}[]
        sep_pl_e_pl_cols      = Vector{Float64}[]
        sep_pl_F_cols         = Vector{Float64}[]
        sep_beam_indices      = Int[]
        for i in 1:num_beams
            if beams.material[i] == :DFT_SEP
                st = material_states.by_beam[beams.ind[i]]::DFT_SEP_States
                push!(sep_se_xi_S_cols,       st.superelasticⁿ.xi_S)
                push!(sep_se_F_cols,          st.superelasticⁿ.F)
                push!(sep_se_eps_cols,        vec(st.superelasticⁿ.eps))
                push!(sep_se_sigma_cols,      vec(st.superelasticⁿ.sigma))
                push!(sep_pl_sigma_cols,      vec(st.Plasticⁿ.sigma))
                push!(sep_pl_epsilon_cols,    vec(st.Plasticⁿ.epsilon))
                push!(sep_pl_epsilon_pl_cols, vec(st.Plasticⁿ.epsilon_pl))
                push!(sep_pl_e_pl_cols,       st.Plasticⁿ.e_pl)
                push!(sep_pl_F_cols,          st.Plasticⁿ.F)
                push!(sep_beam_indices, beams.ind[i])
            end
        end
        open(output_dir * "sep_se_xi_S.txt", "w") do io
            writedlm(io, sep_beam_indices')           # first row: beam indices
            writedlm(io, hcat(sep_se_xi_S_cols...))  # (nG*nS) × num_DFT_SEP_beams
        end
        open(output_dir * "sep_se_F.txt", "w") do io
            writedlm(io, sep_beam_indices')
            writedlm(io, hcat(sep_se_F_cols...))     # (nG*nS) × num_DFT_SEP_beams
        end
        open(output_dir * "sep_se_eps.txt", "w") do io
            writedlm(io, sep_beam_indices')
            writedlm(io, hcat(sep_se_eps_cols...))   # (3*nG*nS) × num_DFT_SEP_beams
        end
        open(output_dir * "sep_se_sigma.txt", "w") do io
            writedlm(io, sep_beam_indices')
            writedlm(io, hcat(sep_se_sigma_cols...)) # (3*nG*nS) × num_DFT_SEP_beams
        end
        open(output_dir * "sep_pl_sigma.txt", "w") do io
            writedlm(io, sep_beam_indices')
            writedlm(io, hcat(sep_pl_sigma_cols...)) # (3*nG*nS) × num_DFT_SEP_beams
        end
        open(output_dir * "sep_pl_epsilon.txt", "w") do io
            writedlm(io, sep_beam_indices')
            writedlm(io, hcat(sep_pl_epsilon_cols...)) # (3*nG*nS) × num_DFT_SEP_beams
        end
        open(output_dir * "sep_pl_epsilon_pl.txt", "w") do io
            writedlm(io, sep_beam_indices')
            writedlm(io, hcat(sep_pl_epsilon_pl_cols...)) # (3*nG*nS) × num_DFT_SEP_beams
        end
        open(output_dir * "sep_pl_e_pl.txt", "w") do io
            writedlm(io, sep_beam_indices')
            writedlm(io, hcat(sep_pl_e_pl_cols...))  # (nG*nS) × num_DFT_SEP_beams
        end
        open(output_dir * "sep_pl_F.txt", "w") do io
            writedlm(io, sep_beam_indices')
            writedlm(io, hcat(sep_pl_F_cols...))     # (nG*nS) × num_DFT_SEP_beams
        end
    end

    # Write data from DFT_PSE beams
    if any(beams.material[i] == :DFT_PSE for i in 1:num_beams)
        pse_se_xi_S_cols      = Vector{Float64}[]
        pse_se_F_cols         = Vector{Float64}[]
        pse_se_eps_cols       = Vector{Float64}[]
        pse_se_sigma_cols     = Vector{Float64}[]
        pse_pl_sigma_cols     = Vector{Float64}[]
        pse_pl_epsilon_cols   = Vector{Float64}[]
        pse_pl_epsilon_pl_cols= Vector{Float64}[]
        pse_pl_e_pl_cols      = Vector{Float64}[]
        pse_pl_F_cols         = Vector{Float64}[]
        pse_beam_indices      = Int[]
        for i in 1:num_beams
            if beams.material[i] == :DFT_PSE
                st = material_states.by_beam[beams.ind[i]]::DFT_PSE_States
                push!(pse_se_xi_S_cols,       st.superelasticⁿ.xi_S)
                push!(pse_se_F_cols,          st.superelasticⁿ.F)
                push!(pse_se_eps_cols,        vec(st.superelasticⁿ.eps))
                push!(pse_se_sigma_cols,      vec(st.superelasticⁿ.sigma))
                push!(pse_pl_sigma_cols,      vec(st.Plasticⁿ.sigma))
                push!(pse_pl_epsilon_cols,    vec(st.Plasticⁿ.epsilon))
                push!(pse_pl_epsilon_pl_cols, vec(st.Plasticⁿ.epsilon_pl))
                push!(pse_pl_e_pl_cols,       st.Plasticⁿ.e_pl)
                push!(pse_pl_F_cols,          st.Plasticⁿ.F)
                push!(pse_beam_indices, beams.ind[i])
            end
        end
        open(output_dir * "pse_se_xi_S.txt", "w") do io
            writedlm(io, pse_beam_indices')
            writedlm(io, hcat(pse_se_xi_S_cols...))  # (nG*nS) × num_DFT_PSE_beams
        end
        open(output_dir * "pse_se_F.txt", "w") do io
            writedlm(io, pse_beam_indices')
            writedlm(io, hcat(pse_se_F_cols...))     # (nG*nS) × num_DFT_PSE_beams
        end
        open(output_dir * "pse_se_eps.txt", "w") do io
            writedlm(io, pse_beam_indices')
            writedlm(io, hcat(pse_se_eps_cols...))   # (3*nG*nS) × num_DFT_PSE_beams
        end
        open(output_dir * "pse_se_sigma.txt", "w") do io
            writedlm(io, pse_beam_indices')
            writedlm(io, hcat(pse_se_sigma_cols...)) # (3*nG*nS) × num_DFT_PSE_beams
        end
        open(output_dir * "pse_pl_sigma.txt", "w") do io
            writedlm(io, pse_beam_indices')
            writedlm(io, hcat(pse_pl_sigma_cols...)) # (3*nG*nS) × num_DFT_PSE_beams
        end
        open(output_dir * "pse_pl_epsilon.txt", "w") do io
            writedlm(io, pse_beam_indices')
            writedlm(io, hcat(pse_pl_epsilon_cols...)) # (3*nG*nS) × num_DFT_PSE_beams
        end
        open(output_dir * "pse_pl_epsilon_pl.txt", "w") do io
            writedlm(io, pse_beam_indices')
            writedlm(io, hcat(pse_pl_epsilon_pl_cols...)) # (3*nG*nS) × num_DFT_PSE_beams
        end
        open(output_dir * "pse_pl_e_pl.txt", "w") do io
            writedlm(io, pse_beam_indices')
            writedlm(io, hcat(pse_pl_e_pl_cols...))  # (nG*nS) × num_DFT_PSE_beams
        end
        open(output_dir * "pse_pl_F.txt", "w") do io
            writedlm(io, pse_beam_indices')
            writedlm(io, hcat(pse_pl_F_cols...))     # (nG*nS) × num_DFT_PSE_beams
        end
    end

    # Write data from DFT_PP beams
    if any(beams.material[i] == :DFT_PP for i in 1:num_beams)
        pp_outer_sigma_cols     = Vector{Float64}[]
        pp_outer_epsilon_cols   = Vector{Float64}[]
        pp_outer_epsilon_pl_cols= Vector{Float64}[]
        pp_outer_e_pl_cols      = Vector{Float64}[]
        pp_outer_F_cols         = Vector{Float64}[]
        pp_inner_sigma_cols     = Vector{Float64}[]
        pp_inner_epsilon_cols   = Vector{Float64}[]
        pp_inner_epsilon_pl_cols= Vector{Float64}[]
        pp_inner_e_pl_cols      = Vector{Float64}[]
        pp_inner_F_cols         = Vector{Float64}[]
        pp_beam_indices         = Int[]
        for i in 1:num_beams
            if beams.material[i] == :DFT_PP
                st = material_states.by_beam[beams.ind[i]]::DFT_PP_States
                push!(pp_outer_sigma_cols,      vec(st.Plasticⁿ_outer.sigma))
                push!(pp_outer_epsilon_cols,    vec(st.Plasticⁿ_outer.epsilon))
                push!(pp_outer_epsilon_pl_cols, vec(st.Plasticⁿ_outer.epsilon_pl))
                push!(pp_outer_e_pl_cols,       st.Plasticⁿ_outer.e_pl)
                push!(pp_outer_F_cols,          st.Plasticⁿ_outer.F)
                push!(pp_inner_sigma_cols,      vec(st.Plasticⁿ_inner.sigma))
                push!(pp_inner_epsilon_cols,    vec(st.Plasticⁿ_inner.epsilon))
                push!(pp_inner_epsilon_pl_cols, vec(st.Plasticⁿ_inner.epsilon_pl))
                push!(pp_inner_e_pl_cols,       st.Plasticⁿ_inner.e_pl)
                push!(pp_inner_F_cols,          st.Plasticⁿ_inner.F)
                push!(pp_beam_indices, beams.ind[i])
            end
        end
        open(output_dir * "pp_outer_sigma.txt", "w") do io
            writedlm(io, pp_beam_indices')
            writedlm(io, hcat(pp_outer_sigma_cols...))      # (3*nG*nS) × num_DFT_PP_beams
        end
        open(output_dir * "pp_outer_epsilon.txt", "w") do io
            writedlm(io, pp_beam_indices')
            writedlm(io, hcat(pp_outer_epsilon_cols...))    # (3*nG*nS) × num_DFT_PP_beams
        end
        open(output_dir * "pp_outer_epsilon_pl.txt", "w") do io
            writedlm(io, pp_beam_indices')
            writedlm(io, hcat(pp_outer_epsilon_pl_cols...)) # (3*nG*nS) × num_DFT_PP_beams
        end
        open(output_dir * "pp_outer_e_pl.txt", "w") do io
            writedlm(io, pp_beam_indices')
            writedlm(io, hcat(pp_outer_e_pl_cols...))       # (nG*nS) × num_DFT_PP_beams
        end
        open(output_dir * "pp_outer_F.txt", "w") do io
            writedlm(io, pp_beam_indices')
            writedlm(io, hcat(pp_outer_F_cols...))          # (nG*nS) × num_DFT_PP_beams
        end
        open(output_dir * "pp_inner_sigma.txt", "w") do io
            writedlm(io, pp_beam_indices')
            writedlm(io, hcat(pp_inner_sigma_cols...))      # (3*nG*nS) × num_DFT_PP_beams
        end
        open(output_dir * "pp_inner_epsilon.txt", "w") do io
            writedlm(io, pp_beam_indices')
            writedlm(io, hcat(pp_inner_epsilon_cols...))    # (3*nG*nS) × num_DFT_PP_beams
        end
        open(output_dir * "pp_inner_epsilon_pl.txt", "w") do io
            writedlm(io, pp_beam_indices')
            writedlm(io, hcat(pp_inner_epsilon_pl_cols...)) # (3*nG*nS) × num_DFT_PP_beams
        end
        open(output_dir * "pp_inner_e_pl.txt", "w") do io
            writedlm(io, pp_beam_indices')
            writedlm(io, hcat(pp_inner_e_pl_cols...))       # (nG*nS) × num_DFT_PP_beams
        end
        open(output_dir * "pp_inner_F.txt", "w") do io
            writedlm(io, pp_beam_indices')
            writedlm(io, hcat(pp_inner_F_cols...))          # (nG*nS) × num_DFT_PP_beams
        end
    end

    # Write data from superelastic beams
    if any(beams.material[i] == :superelastic for i in 1:num_beams)
        xi_S_cols   = Vector{Float64}[]
        F_cols      = Vector{Float64}[]
        eps_cols    = Vector{Float64}[]
        sigma_cols  = Vector{Float64}[]
        beam_indices = Int[]
        for i in 1:num_beams
            if beams.material[i] == :superelastic
                se_state = material_states.by_beam[beams.ind[i]]::SuperelasticStates
                push!(xi_S_cols,  se_state.superelasticⁿ.xi_S)
                push!(F_cols,     se_state.superelasticⁿ.F)
                push!(eps_cols,   vec(se_state.superelasticⁿ.eps))
                push!(sigma_cols, vec(se_state.superelasticⁿ.sigma))
                push!(beam_indices, beams.ind[i])
            end
        end
        open(output_dir * "se_xi_S.txt", "w") do io
            writedlm(io, beam_indices')        # first row: beam indices
            writedlm(io, hcat(xi_S_cols...))  # (nG*nS) × num_superelastic_beams
        end
        open(output_dir * "se_F.txt", "w") do io
            writedlm(io, beam_indices')       # first row: beam indices
            writedlm(io, hcat(F_cols...))     # (nG*nS) × num_superelastic_beams
        end
        open(output_dir * "se_eps.txt", "w") do io
            writedlm(io, beam_indices')       # first row: beam indices
            writedlm(io, hcat(eps_cols...))   # (3*nG*nS) × num_superelastic_beams
        end
        open(output_dir * "se_sigma.txt", "w") do io
            writedlm(io, beam_indices')         # first row: beam indices
            writedlm(io, hcat(sigma_cols...))   # (3*nG*nS) × num_superelastic_beams
        end
    end
end

