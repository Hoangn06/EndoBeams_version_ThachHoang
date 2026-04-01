#----------------------------------
# ELASTIC BEAM PROPERTIES DEFINITION
#----------------------------------

# Elastic beam properties structure
struct ElasticBeamProperties{TK, TJ}
    radius::Float64          # Beam radius
    E::Float64               # Young's modulus
    K̄ⁱⁿᵗ::TK                 # Internal stiffness matrix
    Jᵨ::TJ                   # Rotational inertia matrix
    Aᵨ::Float64              # Cross-sectional area
    damping::Float64         # Damping coefficient
end

# Function to create an elastic beam element
# With matrix connectivity
function construct_beam_elastic_matconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, damping, Rₑ⁰)

    Ei   = E       isa AbstractVector ? E[i]       : E
    νi   = ν       isa AbstractVector ? ν[i]       : ν
    ρi   = ρ       isa AbstractVector ? ρ[i]       : ρ
    ri   = radius  isa AbstractVector ? radius[i]  : radius
    ci   = damping isa AbstractVector ? damping[i] : damping
    mat = material isa AbstractVector ? material[i] : material

    return BeamElastic(i,
        nodes[connectivity[i, 1]],
        nodes[connectivity[i, 2]],
        mat, Ei, νi, ρi, ri, ci,
        Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i, :] : Rₑ⁰)
end


# Function to create an elastic beam element
# With vector connectivity
function construct_beam_elastic_vecconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, damping, Rₑ⁰)

    Ei   = E       isa AbstractVector ? E[i]       : E
    νi   = ν       isa AbstractVector ? ν[i]       : ν
    ρi   = ρ       isa AbstractVector ? ρ[i]       : ρ
    ri   = radius  isa AbstractVector ? radius[i]  : radius
    ci   = damping isa AbstractVector ? damping[i] : damping
    mat = material isa AbstractVector ? material[i] : material

    return BeamElastic(i,
        nodes[connectivity[i][1]],
        nodes[connectivity[i][2]],
        mat,Ei, νi, ρi, ri, ci,
        Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : Rₑ⁰)
end


# Function to create an elastic beam element's data structure
function BeamElastic(ind, node1::NodeBeam, node2::NodeBeam, material, E, ν, ρ, radius, damping, Rₑ⁰)

    i1 = node1.i
    i2 = node2.i
    l₀ = norm(node1.X₀ - node2.X₀)

    beamprops = ElasticBeamProperties(l₀, E, ν, ρ, radius, damping)
    
    return Beam{typeof(beamprops)}(ind, i1, i2, material, l₀, Rₑ⁰, beamprops, zeros(Int, 144), zeros(Int, 144))
end

# Default orientation for elastic when Rₑ⁰ is not provided
BeamElastic(ind, node1::NodeBeam, node2::NodeBeam, material, E, ν, ρ, radius, damping, Rₑ⁰::Nothing) =
    BeamElastic(ind, node1, node2, material, E, ν, ρ, radius, damping, get_Rₑ⁰(node1.X₀, node2.X₀))

# Function to build Elastic Beam properties structure
function ElasticBeamProperties(l₀, E, ν, ρ, radius, damping)
    G = E / (2 * (1 + ν))                        # Shear modulus
    A = pi * radius^2                            # Cross-sectional area
    I₂₂ = pi * radius^4 / 4                      # Second moment of inertia (axis 2)
    I₃₃ = I₂₂                                    # Second moment of inertia (axis 3)
    Iₒ = I₂₂ + I₃₃                               # Polar moment of inertia
    Jᵨ = ρ * Diagonal(Vec3(Iₒ, I₂₂, I₃₃))        # Rotational inertia matrix
    Aᵨ = ρ * A                                   # Area mass density
    K̄ⁱⁿᵗ = K̄ⁱⁿᵗ_beam(E, G, Iₒ, A, I₂₂, I₃₃, l₀)  # Beam stiffness matrix

    return ElasticBeamProperties{typeof(K̄ⁱⁿᵗ), typeof(Jᵨ)}(radius, E, K̄ⁱⁿᵗ, Jᵨ, Aᵨ, damping)
end


