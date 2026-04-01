#----------------------------------
# BEAM STRUCTURES DEFINITION
#----------------------------------

struct Beam{Tp}
    
    ind::Int # index of this beam 
    node1::Int # index of node 1   
    node2::Int # index of node 2
    material::Symbol # material type
    l₀::Float64 # initial beam length
    Rₑ⁰::Mat33{Float64} # initial beam rotation matrix
    properties::Tp # beam internal matrix
    local_sparsity_map::SVector{144, Int} # DOF mapping for local matrix assembly                  
    global_sparsity_map::SVector{144, Int} # DOF mapping for global matrix assembly  

end

#----------------------------------
# FUNCTIONS TO BUILD ELASTIC BEAM
#----------------------------------
function ElasticBeams(nodes, connectivity::AbstractMatrix,E, ν, ρ, radius, damping, Rₑ⁰=nothing)

    material = :elastic

    beams = StructArray((
        construct_beam_elastic_matconn(i, nodes, connectivity, material, E, ν, ρ, radius, damping, Rₑ⁰)
        for i in 1:size(connectivity, 1)
    ))
    return beams
end

function ElasticBeams(nodes, connectivity::AbstractVector, E, ν, ρ, radius, damping, Rₑ⁰=nothing)

    material = :elastic

    beams = StructArray((
            construct_beam_elastic_vecconn(i, nodes, connectivity, material, E, ν, ρ, radius, damping, Rₑ⁰)
            for i in 1:length(connectivity)
        ))
    return beams
end

#----------------------------------
# FUNCTIONS TO BUILD SUPERELASTIC BEAM
#----------------------------------
function SuperElasticBeams(nodes, connectivity::AbstractMatrix, 
    E, ν, ρ, εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ,
    radius, damping,Rₑ⁰=nothing)

    material = :superelastic

    beams = StructArray((
            construct_beam_superelastic_matconn(i, nodes, connectivity, material,
                                    E, ν, ρ, radius, damping,
                                    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, Rₑ⁰)
            for i in 1:size(connectivity, 1)
        ))
    return beams
end

function SuperElasticBeams(nodes, connectivity::AbstractVector,
    E, ν, ρ, εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ,
    radius, damping, Rₑ⁰=nothing)

    material = :superelastic

    beams = StructArray((
            construct_beam_superelastic_vecconn(i, nodes, connectivity, material,
                                    E, ν, ρ, radius, damping,
                                    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, Rₑ⁰)
            for i in 1:length(connectivity)
        ))
    return beams
end

#----------------------------------
# FUNCTIONS TO BUILD PLASTIC BEAM
#----------------------------------
function PlasticBeams(nodes, connectivity::AbstractMatrix,
    E, ν, ρ, 
    y0, K,
    radius, damping, Rₑ⁰=nothing)

    material = :plastic
    beams = StructArray((
        construct_beam_plastic_matconn(i, nodes, connectivity, material, E, ν, ρ, radius, damping, y0, K, Rₑ⁰)
        for i in 1:size(connectivity, 1)
    ))

    return beams
end

function PlasticBeams(nodes, connectivity::AbstractVector,
    E, ν, ρ, 
    y0, K,
    radius, damping, Rₑ⁰=nothing)

    material = :plastic
    beams = StructArray((
        construct_beam_plastic_vecconn(i, nodes, connectivity, material, E, ν, ρ, radius, damping, y0, K, Rₑ⁰)
        for i in 1:length(connectivity)
    ))
    return beams
end

#----------------------------------
# FUNCTIONS TO BUILD DFT SUPERELASTIC - PLASTIC BEAM
#----------------------------------
function DFT_SEPBeams(nodes, connectivity::AbstractMatrix,
    E, ν, ρ, εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ,
    E_inner, ν_inner, ρ_inner, radius_inner, y0, K,
    radius, damping, Nr::Int=4, Nθ::Int=12, Rₑ⁰=nothing)

    beams = StructArray([
            construct_beam_DFT_SEP_matconn(i, nodes, connectivity, :DFT_SEP,
                                    E, ν, ρ, radius, damping,
                                    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, E_inner, ν_inner, ρ_inner, radius_inner, y0, K, Nr, Nθ, Rₑ⁰)
            for i in 1:size(connectivity, 1)
        ])
    return beams
end

function DFT_SEPBeams(nodes, connectivity::AbstractVector,
    E, ν, ρ, εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ,
    E_inner, ν_inner, ρ_inner, radius_inner, y0, K,
    radius, damping, Nr::Int=4, Nθ::Int=12, Rₑ⁰=nothing)

    beams = StructArray([
            construct_beam_DFT_SEP_vecconn(i, nodes, connectivity, :DFT_SEP,
                                    E, ν, ρ, radius, damping,
                                    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, E_inner, ν_inner, ρ_inner, radius_inner, y0, K, Nr, Nθ, Rₑ⁰)
            for i in 1:length(connectivity)
        ])
    return beams
end

#----------------------------------
# FUNCTIONS TO BUILD DFT PLASTIC-SUPERELASTIC BEAM
#----------------------------------
function DFT_PSEBeams(nodes, connectivity::AbstractMatrix,
    E, ν, ρ, y0, K,
    E_inner, ν_inner, ρ_inner, radius_inner,
    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ,
    radius, damping, Nr::Int=4, Nθ::Int=12, Rₑ⁰=nothing)

    beams = StructArray([
            construct_beam_DFT_PSE_matconn(i, nodes, connectivity, :DFT_PSE,
                                    E, ν, ρ, radius, y0, K, damping,
                                    E_inner, ν_inner, ρ_inner, radius_inner,
                                    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, Nr, Nθ, Rₑ⁰)
            for i in 1:size(connectivity, 1)
        ])
    return beams
end

function DFT_PSEBeams(nodes, connectivity::AbstractVector,
    E, ν, ρ, y0, K,
    E_inner, ν_inner, ρ_inner, radius_inner,
    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ,
    radius, damping, Nr::Int=4, Nθ::Int=12, Rₑ⁰=nothing)

    beams = StructArray([
            construct_beam_DFT_PSE_vecconn(i, nodes, connectivity, :DFT_PSE,
                                    E, ν, ρ, radius, y0, K, damping,
                                    E_inner, ν_inner, ρ_inner, radius_inner,
                                    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, Nr, Nθ, Rₑ⁰)
            for i in 1:length(connectivity)
        ])
    return beams
end

#----------------------------------
# FUNCTIONS TO BUILD DFT PLASTIC-PLASTIC BEAM
#----------------------------------
function DFT_PPBeams(nodes, connectivity::AbstractMatrix,
    E, ν, ρ, y0_outer, K_outer,
    E_inner, ν_inner, ρ_inner, radius_inner,
    y0_inner, K_inner,
    radius, damping, Nr::Int=4, Nθ::Int=12, Rₑ⁰=nothing)

    beams = StructArray([
            construct_beam_DFT_PP_matconn(i, nodes, connectivity, :DFT_PP,
                                    E, ν, ρ, radius, y0_outer, K_outer, damping,
                                    E_inner, ν_inner, ρ_inner, radius_inner,
                                    y0_inner, K_inner, Nr, Nθ, Rₑ⁰)
            for i in 1:size(connectivity, 1)
        ])
    return beams
end

function DFT_PPBeams(nodes, connectivity::AbstractVector,
    E, ν, ρ, y0_outer, K_outer,
    E_inner, ν_inner, ρ_inner, radius_inner,
    y0_inner, K_inner,
    radius, damping, Nr::Int=4, Nθ::Int=12, Rₑ⁰=nothing)

    beams = StructArray([
            construct_beam_DFT_PP_vecconn(i, nodes, connectivity, :DFT_PP,
                                    E, ν, ρ, radius, y0_outer, K_outer, damping,
                                    E_inner, ν_inner, ρ_inner, radius_inner,
                                    y0_inner, K_inner, Nr, Nθ, Rₑ⁰)
            for i in 1:length(connectivity)
        ])
    return beams
end
#----------------------------------
# FUNCTIONS TO BUILD DFT ELASTIC-PLASTIC BEAM
#----------------------------------
function DFT_EPBeams(nodes, connectivity::AbstractVector,
    E, ν, ρ,
    E_inner, ν_inner, ρ_inner, radius_inner, y0, K,
    radius, damping, Rₑ⁰=nothing)

    material = :DFT_EP

    beams = StructArray((
            construct_beam_DFT_EP_vecconn(i, nodes, connectivity, material,
                                    E, ν, ρ, radius, damping,
                                    E_inner, ν_inner, ρ_inner, radius_inner, y0, K, Rₑ⁰)
            for i in 1:length(connectivity)
        ))
    return beams
end

function DFT_EPBeams(nodes, connectivity::AbstractMatrix,
    E, ν, ρ,
    E_inner, ν_inner, ρ_inner, radius_inner, y0, K,
    radius, damping, Rₑ⁰=nothing)

    material = :DFT_EP
    beams = StructArray((
        construct_beam_DFT_EP_matconn(i, nodes, connectivity, material, 
                                      E, ν, ρ, radius, damping, 
                                      E_inner, ν_inner, ρ_inner, radius_inner, y0, K, Rₑ⁰)
        for i in 1:size(connectivity, 1)
    ))
    return beams
end


#----------------------------------
# FUNCTIONS TO BUILD DFT PLASTIC-ELASTIC BEAM
#----------------------------------
function DFT_PEBeams(nodes, connectivity::AbstractVector,
    E, ν, ρ,
    radius, y0, K, damping,
    E_inner, ν_inner, ρ_inner, radius_inner,
    Nr::Int=4, Nθ::Int=12, Rₑ⁰=nothing)

    beams = StructArray([
            construct_beam_DFT_PE_vecconn(i, nodes, connectivity, :DFT_PE,
                                    E, ν, ρ, radius, y0, K, damping,
                                    E_inner, ν_inner, ρ_inner, radius_inner, Nr, Nθ, Rₑ⁰)
            for i in 1:length(connectivity)
        ])
    return beams
end

function DFT_PEBeams(nodes, connectivity::AbstractMatrix,
    E, ν, ρ,
    radius, y0, K, damping,
    E_inner, ν_inner, ρ_inner, radius_inner,
    Nr::Int=4, Nθ::Int=12, Rₑ⁰=nothing)

    beams = StructArray([
        construct_beam_DFT_PE_matconn(i, nodes, connectivity, :DFT_PE,
                                      E, ν, ρ, radius, y0, K, damping,
                                      E_inner, ν_inner, ρ_inner, radius_inner, Nr, Nθ, Rₑ⁰)
        for i in 1:size(connectivity, 1)
    ])
    return beams
end


#----------------------------------
# UTILS TO BUILD BEAMS
#----------------------------------

# Function to compute the rotation matrix Rₑ⁰ between two vectors X1 and X2.
# It returns the rotation matrix from X1 to X2, given the direction cosine between them.
function get_Rₑ⁰(X1, X2, Float64=Float64) 
    
    l0 = norm(X2 - X1)  # Length of the vector difference between X2 and X1
    E1 = Vec3(1, 0, 0)  # Unit vector along the x-axis
    E1_0 = (X2 - X1) / l0  # Normalize the difference vector between X1 and X2
    
    v = cross(E1, E1_0)  # Cross product of E1 and the normalized vector
    
    s = norm(v)  # Magnitude of the cross product (sin(theta))
    c = dot(E1, E1_0)  # Dot product between E1 and the normalized vector (cos(theta))
    
    # Handle the special cases where the cosine of the angle is very close to -1 or 1
    if c < -(1 - 2 * eps(Float64))
        # If the cosine is approximately -1, return a 180-degree rotation matrix
        return Mat33{Float64}(-1, 0, 0, 0, 1, 0, 0, 0, -1)
    elseif c > (1 - 2 * eps(Float64))
        # If the cosine is approximately 1, return the identity matrix (no rotation)
        return ID3
    else 
        # Otherwise, compute the Rodrigues' rotation formula for rotation matrix
        Sv = skew(v)  # Skew-symmetric matrix of v
        return ID3 + Sv + ((1 - c) / s^2) * Sv * Sv
    end      
end

# Computes the stiffness matrix components for a beam element.
# This function calculates the axial stiffness, rotational stiffness, and additional rotational stiffness 
# components for a beam element using its material properties and geometric properties.
#
# Inputs:
# - `E`   : Young's Modulus (Pa) of the beam material
# - `G`   : Shear Modulus (Pa) of the beam material
# - `Iₒ`  : Polar moment of inertia (m⁴) of the beam cross-section
# - `A`   : Cross-sectional area (m²) of the beam
# - `I₂₂` : Second moment of area (m⁴) about the second axis (bending stiffness in the y-axis)
# - `I₃₃` : Second moment of area (m⁴) about the third axis (bending stiffness in the z-axis)
# - `l₀`  : Length of the beam element (m)
#
# Returns:
# - `K̄ⁱⁿᵗū`  : Axial stiffness matrix component (a scalar)
# - `K̄ⁱⁿᵗΘ̅`  : Rotational stiffness matrix (a diagonal matrix with rotational stiffness components)
# - `K̄ⁱⁿᵗΘ̅Θ̅`: Additional rotational stiffness matrix (a diagonal matrix with additional rotational stiffness components)
@inline function K̄ⁱⁿᵗ_beam(E, G, Iₒ, A, I₂₂, I₃₃, l₀)
    
    # Compute the axial stiffness matrix component (K̄ⁱⁿᵗū)
    K̄ⁱⁿᵗū = A * E / l₀  # Formula for axial stiffness (A: cross-sectional area, E: Young's Modulus, l₀: beam length)
    
    # Compute the rotational stiffness matrix components (K̄ⁱⁿᵗΘ̅ and K̄ⁱⁿᵗΘ̅Θ̅)
    # K̄ⁱⁿᵗΘ̅ represents the torsion and bending stiffness
    K̄ⁱⁿᵗΘ̅ = Diagonal(@SVector [G * Iₒ / l₀, 4 * E * I₃₃ / l₀, 4 * E * I₂₂ / l₀])
    
    # K̄ⁱⁿᵗΘ̅Θ̅ represents additional rotational stiffness components
    K̄ⁱⁿᵗΘ̅Θ̅ = Diagonal(@SVector [-G * Iₒ / l₀, 2 * E * I₃₃ / l₀, 2 * E * I₂₂ / l₀])
    
    # Return the computed axial and rotational stiffness matrices
    return K̄ⁱⁿᵗū, K̄ⁱⁿᵗΘ̅, K̄ⁱⁿᵗΘ̅Θ̅
end

# Trapezoidal quadrature on [0,2π] and on [radius_inner, radius] for outer layer
function trapezoidal_quadrature_vectors(radius::Float64, radius_inner::Float64, Nr::Int, Nθ::Int)
    Δθ = 2π / Nθ
    zθ = collect(Float64, (0:Nθ) .* Δθ)
    wθ = vcat(0.5 * Δθ, fill(Δθ, Nθ - 1), 0.5 * Δθ)
    Δr = (radius - radius_inner) / Nr
    zr = collect(Float64, radius_inner .+ (0:Nr) .* Δr)
    wr = vcat(0.5 * Δr, fill(Δr, Nr - 1), 0.5 * Δr)
    return zθ, wθ, zr, wr
end