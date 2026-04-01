#----------------------------------
# SUPER-ELASTIC BEAM PROPERTIES DEFINITION
#----------------------------------

# Super-elastic beam properties structure
struct SuperElasticBeamProperties{TK,TJ}
    radius::Float64          # Beam radius
    Eᴬ::Float64               # Austenite Young's modulus
    Eᴹ::Float64               # Martensite Young's modulus
    Gᴬ::Float64               # Austenite Shear modulus
    Gᴹ::Float64               # Martensite Shear modulus
    ν::Float64               # Poisson's ratio

    # Transformation properties
    εL::Float64      # max transformation strain (scalar)
    α::Float64       # Superelasticity parameter (Drucker–Prager coupling)

    # Boundaries for transformation
    R_s_SA::Float64 # start stress for S→A
    R_s_AS::Float64 # start stress for A→S
    R_f_SA::Float64 # finish stress for S→A
    R_f_AS::Float64 # finish stress for A→S

    # Elastic and dynamic properties
    K̄ⁱⁿᵗ::TK                 # Internal elastic stiffness matrix for Austenite
    Jᵨ::TJ                   # Rotational inertia matrix
    Aᵨ::Float64              # Cross-sectional area
    damping::Float64         # Damping coefficient
end

#----------------------------------
# FUNCTIONS TO BUILD SUPER-ELASTIC BEAM STRUCTURES
#----------------------------------

# Create a superelastic beam element
# With matrix connectivity
function construct_beam_superelastic_matconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, damping,
    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS,Eᴹ, Rₑ⁰)

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

    mat = material isa AbstractVector ? material[i] : material

    return BeamSuperelastic(i,
        nodes[connectivity[i, 1]],
        nodes[connectivity[i, 2]],
        mat,Ei, νi, ρi, ri, ci,
        εLi, σS_ASi, σF_ASi, σS_SAi, σF_SAi, σc_SASi,
        Eᴹi,
        Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i, :] : Rₑ⁰)
end

# Function to create a superelastic beam element
# With vector connectivity
function construct_beam_superelastic_vecconn(i, nodes, connectivity, material,
    E, ν, ρ, radius, damping,
    εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, Rₑ⁰)

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

    mat = material isa AbstractVector ? material[i] : material

    return BeamSuperelastic(i,
        nodes[connectivity[i][1]],
        nodes[connectivity[i][2]],
        mat,Ei, νi, ρi, ri, ci,
        εLi, σS_ASi, σF_ASi, σS_SAi, σF_SAi, σc_SASi,
        Eᴹi,
        Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : Rₑ⁰)
end

# Function to create a superelastic beam element's data structure
function BeamSuperelastic(ind, node1::NodeBeam, node2::NodeBeam, material,
        E, ν, ρ, radius, damping,
        εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, Rₑ⁰)

    i1 = node1.i
    i2 = node2.i
    l₀ = norm(node1.X₀ - node2.X₀)

    beamprops = SuperElasticBeamProperties(l₀, E, ν, ρ, radius, damping, εL,
        sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ)
    
    return Beam{typeof(beamprops)}(ind, i1, i2, material, l₀, Rₑ⁰, beamprops, zeros(Int, 144), zeros(Int, 144))
end

# Default orientation for superelastic when Rₑ⁰ is not provided
BeamSuperelastic(ind, node1::NodeBeam, node2::NodeBeam, material, E, ν, ρ, radius, damping,
                εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, Rₑ⁰::Nothing) =
    BeamSuperelastic(ind, node1, node2, material, E, ν, ρ, radius, damping,
                εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, get_Rₑ⁰(node1.X₀, node2.X₀))


# Function to build Super-elastic Beam properties structure
function SuperElasticBeamProperties(l₀, E, ν, ρ, radius, damping, εL, sigma_S_AS, sigma_F_AS, 
        sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ)
    
    Gᴬ = E / (2 * (1 + ν))                        # Shear modulus
    Gᴹ = Eᴹ / (2 * (1 + ν))                        # Shear modulus
    A = pi * radius^2                            # Cross-sectional area
    I₂₂ = pi * radius^4 / 4                      # Second moment of inertia (axis 2)
    I₃₃ = I₂₂                                    # Second moment of inertia (axis 3)
    Iₒ = I₂₂ + I₃₃                               # Polar moment of inertia
    Jᵨ = ρ * Diagonal(Vec3(Iₒ, I₂₂, I₃₃))        # Rotational inertia matrix
    Aᵨ = ρ * A
    K̄ⁱⁿᵗ = K̄ⁱⁿᵗ_beam(E, Gᴬ, Iₒ, A, I₂₂, I₃₃, l₀)  # Beam stiffness matrix

    α = sqrt(2/3)*((sigma_c_SAS-sigma_S_AS)/(sigma_c_SAS+sigma_S_AS)) # Asymmetry parameter

    # Compute boundaries for Drucker-Prager functions
    rfac = α + sqrt(2.0/3.0)
    R_s_AS = rfac * sigma_S_AS
    R_f_AS = rfac * sigma_F_AS
    R_s_SA = rfac * sigma_S_SA
    R_f_SA = rfac * sigma_F_SA

    return SuperElasticBeamProperties{typeof(K̄ⁱⁿᵗ), typeof(Jᵨ)}(
        radius, E, Eᴹ, Gᴬ, Gᴹ, ν,
        εL, α, 
        R_s_SA, R_s_AS, R_f_SA, R_f_AS, 
        K̄ⁱⁿᵗ, Jᵨ, Aᵨ, damping)
end

