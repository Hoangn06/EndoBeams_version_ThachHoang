#----------------------------------
# DFT ELASTIC-PLASTIC BEAM PROPERTIES DEFINITION
# (Composite between elastic outer layer + plastic inner layer)
#----------------------------------

# DFT Elastic-Plastic beam properties structure
struct DFT_EPBeamProperties{TO,TJ}
    # Elastic properties for the outer layer
    radius::Float64          # Beam radius
    E_outer::Float64         # Outer layer Young's modulus
    ν_outer::Float64         # Outer layer Poisson's ratio
    G_outer::Float64         # Outer layer Shear modulus

    # Plastic properties for the inner layer
    radius_inner::Float64    # Inner radius
    E_inner::Float64         # Inner layer Young's modulus
    ν_inner::Float64         # Inner layer Poisson's ratio
    G_inner::Float64         # Inner layer Shear modulus
    y0::Float64              # Inner layer yield stress
    K::Float64               # Inner layer plastic hardening modulus
    

    K̄ⁱⁿᵗ_outer::TO           # Internal elastic stiffness matrix for outer layer
    Jᵨ::TJ                   # Rotational inertia matrix
    Aᵨ::Float64              # Cross-sectional area
    damping::Float64         # Damping coefficient
end

#----------------------------------
# FUNCTIONS TO BUILD DFT Elastic-Plastic BEAM STRUCTURES
#----------------------------------

# Create a DFT Elastic-Plastic beam element
# With matrix connectivity
function construct_beam_DFT_EP_matconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, damping,
    E_inner, ν_inner, ρ_inner, radius_inner, y0, K, Rₑ⁰)

    Ei   = E       isa AbstractVector ? E[i]       : E
    νi   = ν       isa AbstractVector ? ν[i]       : ν
    ρi   = ρ       isa AbstractVector ? ρ[i]       : ρ
    ri   = radius  isa AbstractVector ? radius[i]  : radius
    ci   = damping isa AbstractVector ? damping[i] : damping

    E_inneri = E_inner    isa AbstractVector ? E_inner[i]    : E_inner
    ν_inneri = ν_inner    isa AbstractVector ? ν_inner[i]    : ν_inner
    ρ_inneri = ρ_inner    isa AbstractVector ? ρ_inner[i]    : ρ_inner
    radius_inneri = radius_inner isa AbstractVector ? radius_inner[i] : radius_inner
    y0i = y0 isa AbstractVector ? y0[i] : y0
    Ki = K isa AbstractVector ? K[i] : K

    mat = material isa AbstractVector ? material[i] : material

    return BeamDFT_EP(i,
        nodes[connectivity[i, 1]],
        nodes[connectivity[i, 2]],
        mat,Ei, νi, ρi, ri, ci,
        E_inneri, ν_inneri, ρ_inneri, radius_inneri, y0i, Ki,
        Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i, :] : Rₑ⁰)
end

# Function to create a DFT Elastic-Plastic beam element
# With vector connectivity
function construct_beam_DFT_EP_vecconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, damping,
    E_inner, ν_inner, ρ_inner, radius_inner, y0, K, Rₑ⁰)

    Ei   = E       isa AbstractVector ? E[i]       : E
    νi   = ν       isa AbstractVector ? ν[i]       : ν
    ρi   = ρ       isa AbstractVector ? ρ[i]       : ρ
    ri   = radius  isa AbstractVector ? radius[i]  : radius
    ci   = damping isa AbstractVector ? damping[i] : damping

    E_inneri = E_inner    isa AbstractVector ? E_inner[i]    : E_inner
    ν_inneri = ν_inner    isa AbstractVector ? ν_inner[i]    : ν_inner
    ρ_inneri = ρ_inner    isa AbstractVector ? ρ_inner[i]    : ρ_inner
    radius_inneri = radius_inner isa AbstractVector ? radius_inner[i] : radius_inner
    y0i = y0 isa AbstractVector ? y0[i] : y0
    Ki = K isa AbstractVector ? K[i] : K

    mat = material isa AbstractVector ? material[i] : material

    return BeamDFT_EP(i,
        nodes[connectivity[i][1]],
        nodes[connectivity[i][2]],
        mat,Ei, νi, ρi, ri, ci,
        E_inneri, ν_inneri, ρ_inneri, radius_inneri, y0i, Ki,
        Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : Rₑ⁰)
end

# Function to create a DFT Elastic-Plastic beam element's data structure
function BeamDFT_EP(ind, node1::NodeBeam, node2::NodeBeam, material,
        E, ν, ρ, radius, damping,
        E_inner, ν_inner, ρ_inner, radius_inner, y0, K, Rₑ⁰)

    i1 = node1.i
    i2 = node2.i
    l₀ = norm(node1.X₀ - node2.X₀)

    beamprops = DFT_EPBeamProperties(l₀, E, ν, ρ, radius, damping, E_inner, ν_inner, ρ_inner, radius_inner, y0, K)
    
    return Beam{typeof(beamprops)}(ind, i1, i2, material, l₀, Rₑ⁰, beamprops, zeros(Int, 144), zeros(Int, 144))
end

# Default orientation for DFT Elastic-Plastic when Rₑ⁰ is not provided
BeamDFT_EP(ind, node1::NodeBeam, node2::NodeBeam, material, E, ν, ρ, radius, damping,
                E_inner, ν_inner, ρ_inner, radius_inner, y0, K, Rₑ⁰::Nothing) =
    BeamDFT_EP(ind, node1, node2, material, E, ν, ρ, radius, damping,
                E_inner, ν_inner, ρ_inner, radius_inner, y0, K, get_Rₑ⁰(node1.X₀, node2.X₀))


# Function to build DFT Elastic-Plastic Beam properties structure
function DFT_EPBeamProperties(l₀, E, ν, ρ, radius, damping, E_inner, ν_inner, ρ_inner, radius_inner, y0, K)
    
    # Compute the elastic outer layer properties
    G_outer = E / (2 * (1 + ν))                        # Shear modulus
    # Compute the cross section area and moments of inertia for the outer layer
    A_outer = pi * ((radius)^2 - (radius_inner)^2)
    I₂₂_outer = pi * ((radius)^4 - (radius_inner)^4) / 4
    I₃₃_outer = I₂₂_outer
    Iₒ_outer = I₂₂_outer + I₃₃_outer
    Jᵨ_outer = ρ * Diagonal(Vec3(Iₒ_outer, I₂₂_outer, I₃₃_outer))
    Aᵨ_outer = ρ * A_outer
    # Compute the tangent stiffness matrix for the outer layer in case of Austenite
    K̄ⁱⁿᵗ_outer = K̄ⁱⁿᵗ_beam(E, G_outer, Iₒ_outer, A_outer, I₂₂_outer, I₃₃_outer, l₀)

    # Compute the elastic properties for the inner layer
    G_inner = E_inner / (2 * (1 + ν_inner))
    # Compute the cross section area and moments of inertia for the inner layer
    A_inner = pi * radius_inner^2
    I₂₂_inner = pi * radius_inner^4 / 4
    I₃₃_inner = I₂₂_inner
    Iₒ_inner = I₂₂_inner + I₃₃_inner
    Jᵨ_inner = ρ_inner * Diagonal(Vec3(Iₒ_inner, I₂₂_inner, I₃₃_inner))
    Aᵨ_inner = ρ_inner * A_inner

    # Total Rotational inertia matrix 
    Jᵨ = Jᵨ_outer + Jᵨ_inner
    # Total cross section area
    Aᵨ = Aᵨ_outer + Aᵨ_inner

    return DFT_EPBeamProperties{typeof(K̄ⁱⁿᵗ_outer), typeof(Jᵨ)}(
        radius, 
        E, ν, G_outer, 
        radius_inner,
        E_inner, ν_inner, G_inner, y0, K,
        K̄ⁱⁿᵗ_outer, Jᵨ, Aᵨ, damping)
end

