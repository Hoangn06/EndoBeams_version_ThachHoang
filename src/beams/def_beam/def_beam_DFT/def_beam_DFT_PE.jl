#----------------------------------
# DFT PLASTIC-ELASTIC BEAM PROPERTIES DEFINITION
# (Composite between plastic outer layer + elastic inner layer)
#----------------------------------

# DFT Plastic-Elastic beam properties structure
struct DFT_PEBeamProperties{TO,TJ}
    # Elastic properties for the outer layer
    radius::Float64          # Beam radius
    E_outer::Float64         # Outer layer Young's modulus
    ν_outer::Float64         # Outer layer Poisson's ratio
    G_outer::Float64         # Outer layer Shear modulus
    y0::Float64              # Outer layer yield stress
    K::Float64               # Outer layer plastic hardening modulus

    # Plastic properties for the inner layer
    radius_inner::Float64    # Inner radius
    E_inner::Float64         # Inner layer Young's modulus
    ν_inner::Float64         # Inner layer Poisson's ratio
    G_inner::Float64         # Inner layer Shear modulus
    
    # Elastic and dynamic properties
    K̄ⁱⁿᵗ_inner::TO           # Internal elastic stiffness matrix for inner layer
    Jᵨ::TJ                   # Rotational inertia matrix
    Aᵨ::Float64              # Cross-sectional area
    damping::Float64         # Damping coefficient

    # Quadrature nodes and weights for the outer layer cross-section
    zθ::Vector{Float64}
    wθ::Vector{Float64}
    zr::Vector{Float64}
    wr::Vector{Float64}
end

#----------------------------------
# FUNCTIONS TO BUILD DFT Elastic-Plastic BEAM STRUCTURES
#----------------------------------

# Create a DFT Plastic-Elastic beam element
# With matrix connectivity
function construct_beam_DFT_PE_matconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, y0, K, damping,
    E_inner, ν_inner, ρ_inner, radius_inner, Nr, Nθ, Rₑ⁰)

    Ei   = E       isa AbstractVector ? E[i]       : E
    νi   = ν       isa AbstractVector ? ν[i]       : ν
    ρi   = ρ       isa AbstractVector ? ρ[i]       : ρ
    ri   = radius  isa AbstractVector ? radius[i]  : radius
    y0i = y0 isa AbstractVector ? y0[i] : y0
    Ki = K isa AbstractVector ? K[i] : K
    ci   = damping isa AbstractVector ? damping[i] : damping

    E_inneri = E_inner    isa AbstractVector ? E_inner[i]    : E_inner
    ν_inneri = ν_inner    isa AbstractVector ? ν_inner[i]    : ν_inner
    ρ_inneri = ρ_inner    isa AbstractVector ? ρ_inner[i]    : ρ_inner
    radius_inneri = radius_inner isa AbstractVector ? radius_inner[i] : radius_inner

    mat = material isa AbstractVector ? material[i] : material

    return BeamDFT_PE(i,
        nodes[connectivity[i, 1]],
        nodes[connectivity[i, 2]],
        mat,Ei, νi, ρi, ri, y0i, Ki, ci,
        E_inneri, ν_inneri, ρ_inneri, radius_inneri,
        Nr, Nθ, Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i, :] : Rₑ⁰)
end

# Function to create a DFT Plastic-Elastic beam element
# With vector connectivity
function construct_beam_DFT_PE_vecconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, y0, K, damping,
    E_inner, ν_inner, ρ_inner, radius_inner, Nr, Nθ, Rₑ⁰)

    Ei   = E       isa AbstractVector ? E[i]       : E
    νi   = ν       isa AbstractVector ? ν[i]       : ν
    ρi   = ρ       isa AbstractVector ? ρ[i]       : ρ
    ri   = radius  isa AbstractVector ? radius[i]  : radius
    y0i = y0 isa AbstractVector ? y0[i] : y0
    Ki = K isa AbstractVector ? K[i] : K
    ci   = damping isa AbstractVector ? damping[i] : damping

    E_inneri = E_inner    isa AbstractVector ? E_inner[i]    : E_inner
    ν_inneri = ν_inner    isa AbstractVector ? ν_inner[i]    : ν_inner
    ρ_inneri = ρ_inner    isa AbstractVector ? ρ_inner[i]    : ρ_inner
    radius_inneri = radius_inner isa AbstractVector ? radius_inner[i] : radius_inner

    mat = material isa AbstractVector ? material[i] : material

    return BeamDFT_PE(i,
        nodes[connectivity[i][1]],
        nodes[connectivity[i][2]],
        mat,Ei, νi, ρi, ri, y0i, Ki, ci,
        E_inneri, ν_inneri, ρ_inneri, radius_inneri,
        Nr, Nθ, Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : Rₑ⁰)
end

# Function to create a DFT Plastic-Elastic beam element's data structure
function BeamDFT_PE(ind, node1::NodeBeam, node2::NodeBeam, material,
        E, ν, ρ, radius, y0, K, damping,
        E_inner, ν_inner, ρ_inner, radius_inner, Nr::Int, Nθ::Int, Rₑ⁰)

    i1 = node1.i
    i2 = node2.i
    l₀ = norm(node1.X₀ - node2.X₀)

    beamprops = DFT_PEBeamProperties(l₀, E, ν, ρ, radius, y0, K, damping, E_inner, ν_inner, ρ_inner, radius_inner, Nr, Nθ)
    
    return Beam{typeof(beamprops)}(ind, i1, i2, material, l₀, Rₑ⁰, beamprops, zeros(Int, 144), zeros(Int, 144))
end

# Default orientation for DFT Plastic-Elastic when Rₑ⁰ is not provided 
BeamDFT_PE(ind, node1::NodeBeam, node2::NodeBeam, material, E, ν, ρ, radius, y0, K, damping,
                E_inner, ν_inner, ρ_inner, radius_inner, Nr::Int, Nθ::Int, Rₑ⁰::Nothing) =
    BeamDFT_PE(ind, node1, node2, material, E, ν, ρ, radius, y0, K, damping,
                E_inner, ν_inner, ρ_inner, radius_inner, Nr, Nθ, get_Rₑ⁰(node1.X₀, node2.X₀))

# Function to build DFT Plastic-Elastic Beam properties structure
function DFT_PEBeamProperties(l₀, E, ν, ρ, radius, y0, K, damping, E_inner, ν_inner, ρ_inner, radius_inner, Nr::Int, Nθ::Int)
    
    # Compute the elastic outer layer properties
    G_outer = E / (2 * (1 + ν))                        # Shear modulus
    # Compute the cross section area and moments of inertia for the outer layer
    A_outer = pi * ((radius)^2 - (radius_inner)^2)
    I₂₂_outer = pi * ((radius)^4 - (radius_inner)^4) / 4
    I₃₃_outer = I₂₂_outer
    Iₒ_outer = I₂₂_outer + I₃₃_outer
    Jᵨ_outer = ρ * Diagonal(Vec3(Iₒ_outer, I₂₂_outer, I₃₃_outer))
    Aᵨ_outer = ρ * A_outer

    # Compute the elastic properties for the inner layer
    G_inner = E_inner / (2 * (1 + ν_inner))
    # Compute the cross section area and moments of inertia for the inner layer
    A_inner = pi * radius_inner^2
    I₂₂_inner = pi * radius_inner^4 / 4
    I₃₃_inner = I₂₂_inner
    Iₒ_inner = I₂₂_inner + I₃₃_inner
    Jᵨ_inner = ρ_inner * Diagonal(Vec3(Iₒ_inner, I₂₂_inner, I₃₃_inner))
    Aᵨ_inner = ρ_inner * A_inner
    # Compute the tangent stiffness matrix for the inner layer in case of Austenite
    K̄ⁱⁿᵗ_inner = K̄ⁱⁿᵗ_beam(E_inner, G_inner, Iₒ_inner, A_inner, I₂₂_inner, I₃₃_inner, l₀)

    # Total Rotational inertia matrix 
    Jᵨ = Jᵨ_outer + Jᵨ_inner
    # Total cross section area
    Aᵨ = Aᵨ_outer + Aᵨ_inner

    zθ, wθ, zr, wr = trapezoidal_quadrature_vectors(radius, radius_inner, Nr, Nθ)

    return DFT_PEBeamProperties{typeof(K̄ⁱⁿᵗ_inner), typeof(Jᵨ)}(
        radius, 
        E, ν, G_outer, y0, K,
        radius_inner,
        E_inner, ν_inner, G_inner,
        K̄ⁱⁿᵗ_inner, Jᵨ, Aᵨ, damping,
        zθ, wθ, zr, wr)
end


