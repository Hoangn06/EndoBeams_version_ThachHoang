#----------------------------------
# PLASTIC BEAM PROPERTIES DEFINITION
#----------------------------------

# Plastic beam properties structure
struct PlasticBeamProperties{TK,TJ}
    radius::Float64          # Beam radius
    E::Float64               # Elastic Young's modulus
    G::Float64               # Shear modulus
    ν::Float64               # Poisson's ratio

    y0::Float64
    K::Float64

    # Elastic and dynamic properties
    K̄ⁱⁿᵗ::TK                 # Internal elastic stiffness matrix for Austenite
    Jᵨ::TJ                   # Rotational inertia matrix
    Aᵨ::Float64              # Cross-sectional area
    damping::Float64         # Damping coefficient
end

#----------------------------------
# FUNCTIONS TO BUILD PLASTIC BEAM STRUCTURES
#----------------------------------

# Create a plastic beam element
# With matrix connectivity
function construct_beam_plastic_matconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, damping,
    y0, K, Rₑ⁰)

    Ei   = E       isa AbstractVector ? E[i]       : E
    νi   = ν       isa AbstractVector ? ν[i]       : ν
    ρi   = ρ       isa AbstractVector ? ρ[i]       : ρ
    ri   = radius  isa AbstractVector ? radius[i]  : radius
    ci   = damping isa AbstractVector ? damping[i] : damping

    y0i = y0 isa AbstractVector ? y0[i] : y0
    Ki = K isa AbstractVector ? K[i] : K

    mat = material isa AbstractVector ? material[i] : material

    return BeamPlastic(i,
        nodes[connectivity[i, 1]],
        nodes[connectivity[i, 2]],
        mat,Ei, νi, ρi, ri, ci,
        y0i, Ki,
        Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i, :] : Rₑ⁰)
end

# Function to create a plastic beam element
# With vector connectivity
function construct_beam_plastic_vecconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, damping,
    y0, K, Rₑ⁰)

    Ei   = E       isa AbstractVector ? E[i]       : E
    νi   = ν       isa AbstractVector ? ν[i]       : ν
    ρi   = ρ       isa AbstractVector ? ρ[i]       : ρ
    ri   = radius  isa AbstractVector ? radius[i]  : radius
    ci   = damping isa AbstractVector ? damping[i] : damping

    y0i = y0 isa AbstractVector ? y0[i] : y0
    Ki = K isa AbstractVector ? K[i] : K

    mat = material isa AbstractVector ? material[i] : material

    return BeamPlastic(i,
        nodes[connectivity[i][1]],
        nodes[connectivity[i][2]],
        mat,Ei, νi, ρi, ri, ci,
        y0i, Ki,
        Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : Rₑ⁰)
end

# Function to create a plastic beam element's data structure
function BeamPlastic(ind, node1::NodeBeam, node2::NodeBeam, material,
        E, ν, ρ, radius, damping,
        y0, K, Rₑ⁰)

    i1 = node1.i
    i2 = node2.i
    l₀ = norm(node1.X₀ - node2.X₀)

    beamprops = PlasticBeamProperties(l₀, E, ν, ρ, radius, damping, y0, K)
    
    return Beam{typeof(beamprops)}(ind, i1, i2, material, l₀, Rₑ⁰, beamprops, zeros(Int, 144), zeros(Int, 144))
end

# Default orientation for plastic when Rₑ⁰ is not provided
BeamPlastic(ind, node1::NodeBeam, node2::NodeBeam, material, E, ν, ρ, radius, damping,
                y0, K, Rₑ⁰::Nothing) =
    BeamPlastic(ind, node1, node2, material, E, ν, ρ, radius, damping,
                y0, K, get_Rₑ⁰(node1.X₀, node2.X₀))


# Function to build Super-elastic Beam properties structure
function PlasticBeamProperties(l₀, E, ν, ρ, radius, damping, y0, K)
    
    G = E / (2 * (1 + ν))                        # Shear modulus
    A = pi * radius^2                            # Cross-sectional area
    I₂₂ = pi * radius^4 / 4                      # Second moment of inertia (axis 2)
    I₃₃ = I₂₂                                    # Second moment of inertia (axis 3)
    Iₒ = I₂₂ + I₃₃                               # Polar moment of inertia
    Jᵨ = ρ * Diagonal(Vec3(Iₒ, I₂₂, I₃₃))        # Rotational inertia matrix
    Aᵨ = ρ * A
    K̄ⁱⁿᵗ = K̄ⁱⁿᵗ_beam(E, G, Iₒ, A, I₂₂, I₃₃, l₀)  # Beam stiffness matrix

    return PlasticBeamProperties{typeof(K̄ⁱⁿᵗ), typeof(Jᵨ)}(
        radius, E, G, ν,
        y0, K,
        K̄ⁱⁿᵗ, Jᵨ, Aᵨ, damping)
end

