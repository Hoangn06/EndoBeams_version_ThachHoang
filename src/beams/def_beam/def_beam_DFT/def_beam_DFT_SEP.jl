#----------------------------------
# DFT-SEP BEAM PROPERTIES DEFINITION
# (Composite between superelastic outer layer + plastic inner layer)
#----------------------------------

# DFT-SEP beam properties structure
struct DFT_SEPBeamProperties{TO,TJ}
    
    # Superelastic properties for the outer layer
    radius::Float64          # Beam radius
    Eᴬ::Float64               # Austenite Young's modulus
    Eᴹ::Float64               # Martensite Young's modulus
    ν::Float64               # Poisson's ratio
    Gᴬ::Float64               # Austenite Shear modulus
    Gᴹ::Float64               # Martensite Shear modulus
    εL::Float64      # max transformation strain (scalar)
    α::Float64       # Superelasticity parameter (Drucker–Prager coupling)
    R_s_SA::Float64 # start stress for S→A
    R_s_AS::Float64 # start stress for A→S
    R_f_SA::Float64 # finish stress for S→A
    R_f_AS::Float64 # finish stress for A→S

    # Plastic properties for the inner layer
    radius_inner::Float64           # Inner layer radius
    E_inner::Float64               # Inner layer Young's modulus
    ν_inner::Float64               # Inner layer Poisson's ratio
    G_inner::Float64               # Inner layer Shear modulus
    y0::Float64                    # Yield stress
    K::Float64                      # Plastic hardening modulus

    # Elastic and dynamic properties
    K̄ⁱⁿᵗ_outer::TO                 # Internal elastic stiffness matrix for outer layer
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
# FUNCTIONS TO BUILD DFT-SEP BEAM STRUCTURES
#----------------------------------

# Create a DFT-SEP beam element
# With matrix connectivity
function construct_beam_DFT_SEP_matconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, damping,
    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS,Eᴹ, E_inner, ν_inner, ρ_inner,radius_inner, y0, K, Nr, Nθ, Rₑ⁰)

    Ei   = E       isa AbstractVector ? E[i]       : E
    νi   = ν       isa AbstractVector ? ν[i]       : ν
    ρi   = ρ       isa AbstractVector ? ρ[i]       : ρ
    ri   = radius  isa AbstractVector ? radius[i]  : radius
    ci   = damping isa AbstractVector ? damping[i] : damping

    εLi     = εL            isa AbstractVector ? εL[i]            : εL
    σS_ASi  = sigma_S_AS    isa AbstractVector ? sigma_S_AS[i]    : sigma_S_AS
    σF_ASi  = sigma_F_AS    isa AbstractVector ? sigma_F_AS[i]    : sigma_F_AS
    σS_SAi  = sigma_S_SA    isa AbstractVector ? sigma_S_SA[i]    : sigma_S_SA
    σF_SAi  = sigma_F_SA    isa AbstractVector ? sigma_F_SA[i]    : sigma_F_SA
    σc_SASi = sigma_c_SAS   isa AbstractVector ? sigma_c_SAS[i]   : sigma_c_SAS
    Eᴹi     = Eᴹ            isa AbstractVector ? Eᴹ[i]            : Eᴹ

    E_inneri = E_inner    isa AbstractVector ? E_inner[i]    : E_inner
    ν_inneri = ν_inner    isa AbstractVector ? ν_inner[i]    : ν_inner
    ρ_inneri = ρ_inner    isa AbstractVector ? ρ_inner[i]    : ρ_inner
    radius_inneri = radius_inner isa AbstractVector ? radius_inner[i] : radius_inner
    y0i = y0 isa AbstractVector ? y0[i] : y0
    Ki = K isa AbstractVector ? K[i] : K
    

    mat = material isa AbstractVector ? material[i] : material

    return BeamDFT_SEP(i,
        nodes[connectivity[i, 1]],
        nodes[connectivity[i, 2]],
        mat,Ei, νi, ρi, ri, ci,
        εLi, σS_ASi, σF_ASi, σS_SAi, σF_SAi, σc_SASi,
        Eᴹi, E_inneri, ν_inneri, ρ_inneri,radius_inneri, y0i, Ki,
        Nr, Nθ, Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i, :] : Rₑ⁰)
end

# Function to create a DFT-SEP beam element
# With vector connectivity
function construct_beam_DFT_SEP_vecconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, damping,
    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, E_inner, ν_inner, ρ_inner,radius_inner, y0, K, Nr, Nθ, Rₑ⁰)

    Ei   = E       isa AbstractVector ? E[i]       : E
    νi   = ν       isa AbstractVector ? ν[i]       : ν
    ρi   = ρ       isa AbstractVector ? ρ[i]       : ρ
    ri   = radius  isa AbstractVector ? radius[i]  : radius
    ci   = damping isa AbstractVector ? damping[i] : damping

    εLi     = εL            isa AbstractVector ? εL[i]            : εL
    σS_ASi  = sigma_S_AS    isa AbstractVector ? sigma_S_AS[i]    : sigma_S_AS
    σF_ASi  = sigma_F_AS    isa AbstractVector ? sigma_F_AS[i]    : sigma_F_AS
    σS_SAi  = sigma_S_SA    isa AbstractVector ? sigma_S_SA[i]    : sigma_S_SA
    σF_SAi  = sigma_F_SA    isa AbstractVector ? sigma_F_SA[i]    : sigma_F_SA
    σc_SASi = sigma_c_SAS   isa AbstractVector ? sigma_c_SAS[i]   : sigma_c_SAS
    Eᴹi     = Eᴹ            isa AbstractVector ? Eᴹ[i]            : Eᴹ

    E_inneri = E_inner    isa AbstractVector ? E_inner[i]    : E_inner
    ν_inneri = ν_inner    isa AbstractVector ? ν_inner[i]    : ν_inner
    ρ_inneri = ρ_inner    isa AbstractVector ? ρ_inner[i]    : ρ_inner
    radius_inneri = radius_inner isa AbstractVector ? radius_inner[i] : radius_inner
    y0i = y0 isa AbstractVector ? y0[i] : y0
    Ki = K isa AbstractVector ? K[i] : K

    mat = material isa AbstractVector ? material[i] : material

    return BeamDFT_SEP(i,
        nodes[connectivity[i][1]],
        nodes[connectivity[i][2]],
        mat,Ei, νi, ρi, ri, ci,
        εLi, σS_ASi, σF_ASi, σS_SAi, σF_SAi, σc_SASi,
        Eᴹi, E_inneri, ν_inneri, ρ_inneri,radius_inneri, y0i, Ki,
        Nr, Nθ, Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : Rₑ⁰)
end

# Function to create a DFT-SEP beam element's data structure
function BeamDFT_SEP(ind, node1::NodeBeam, node2::NodeBeam, material,
        E, ν, ρ, radius, damping,
        εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, E_inner, ν_inner, ρ_inner,radius_inner, y0, K, Nr::Int, Nθ::Int, Rₑ⁰)

    i1 = node1.i
    i2 = node2.i
    l₀ = norm(node1.X₀ - node2.X₀)

    beamprops = DFT_SEPBeamProperties(l₀, E, ν, ρ, radius, damping, εL,
        sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, E_inner, ν_inner, ρ_inner,radius_inner, y0, K, Nr, Nθ)
    
    return Beam{typeof(beamprops)}(ind, i1, i2, material, l₀, Rₑ⁰, beamprops, zeros(Int, 144), zeros(Int, 144))
end

# Default orientation when Rₑ⁰ is not provided (Rₑ⁰ must be last)
BeamDFT_SEP(ind, node1::NodeBeam, node2::NodeBeam, material, E, ν, ρ, radius, damping,
                εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, E_inner, ν_inner, ρ_inner,radius_inner, y0, K, Nr::Int, Nθ::Int, Rₑ⁰::Nothing) =
    BeamDFT_SEP(ind, node1, node2, material, E, ν, ρ, radius, damping,
                εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, E_inner, ν_inner, ρ_inner,radius_inner, y0, K, Nr, Nθ, get_Rₑ⁰(node1.X₀, node2.X₀))


# Function to build DFT-SEP Beam properties structure
function DFT_SEPBeamProperties(l₀, E, ν, ρ, radius, damping, εL, sigma_S_AS, sigma_F_AS, 
        sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, E_inner, ν_inner, ρ_inner,radius_inner, y0, K, Nr::Int, Nθ::Int)
    
    # Compute the superelastic outer layer properties
    Gᴬ = E / (2 * (1 + ν))                        # Shear modulus
    Gᴹ = Eᴹ / (2 * (1 + ν))                        # Shear modulus
    # Compute the asymmetry parameter
    α = sqrt(2/3)*((sigma_c_SAS-sigma_S_AS)/(sigma_c_SAS+sigma_S_AS)) # Asymmetry parameter
    # Compute boundaries for Drucker-Prager functions
    rfac = α + sqrt(2.0/3.0)
    R_s_AS = rfac * sigma_S_AS
    R_f_AS = rfac * sigma_F_AS
    R_s_SA = rfac * sigma_S_SA
    R_f_SA = rfac * sigma_F_SA
    # Compute the cross section area and moments of inertia for the outer layer
    A_outer = pi * ((radius)^2 - (radius_inner)^2)
    I₂₂_outer = pi * ((radius)^4 - (radius_inner)^4) / 4
    I₃₃_outer = I₂₂_outer
    Iₒ_outer = I₂₂_outer + I₃₃_outer
    Jᵨ_outer = ρ * Diagonal(Vec3(Iₒ_outer, I₂₂_outer, I₃₃_outer))
    Aᵨ_outer = ρ * A_outer
    # Compute the tangent stiffness matrix for the outer layer in case of Austenite
    K̄ⁱⁿᵗ_outer = K̄ⁱⁿᵗ_beam(E, Gᴬ, Iₒ_outer, A_outer, I₂₂_outer, I₃₃_outer, l₀)

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

    zθ, wθ, zr, wr = trapezoidal_quadrature_vectors(radius, radius_inner, Nr, Nθ)

    return DFT_SEPBeamProperties{typeof(K̄ⁱⁿᵗ_outer), typeof(Jᵨ)}(
        radius,
        E, Eᴹ, ν, Gᴬ, Gᴹ,
        εL, α, 
        R_s_SA, R_s_AS, R_f_SA, R_f_AS,
        radius_inner,
        E_inner, ν_inner, G_inner, y0, K,
        K̄ⁱⁿᵗ_outer, Jᵨ, Aᵨ, damping,
        zθ, wθ, zr, wr)
end

