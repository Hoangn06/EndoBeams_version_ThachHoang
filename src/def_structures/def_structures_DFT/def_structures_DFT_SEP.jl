#----------------------------------
# Structures to store the state of the DFT beam
# (Composite between superelastic outer tube + plastic inner fibers)

#----------------------------------
struct DFT_SEP_States
    superelasticⁿ::SuperelasticState   # Superelastic state at the current time step
    superelasticⁿ⁺¹::SuperelasticState # Superelastic state at the next time step
    Plasticⁿ::PlasticState             # Plastic state at the current time step
    Plasticⁿ⁺¹::PlasticState           # Plastic state at the next time step
end

#----------------------------------
# Functions to initialize the DFT structures
#----------------------------------

function DFT_SEP_States(params::SimulationParams)
    nG = params.nᴳ_beams
    # Outer (superelastic) integration uses trapezoidal grid in (r, θ)
    n_outer = (params.Nr + 1) * (params.Nθ + 1)
    n_total = nG * n_outer

    # Inner (plastic) integration uses section quadrature (nˢ × nˢ)
    n_inner = params.nˢ_beams * params.nˢ_beams
    n_total_inner = nG * n_inner

    superelasticⁿ = SuperelasticState(zeros(n_total), zeros(n_total), zeros(n_total), zeros(n_total), zeros(n_total), zeros(3, n_total), zeros(3, n_total))
    superelasticⁿ⁺¹ = deepcopy(superelasticⁿ)

    Plasticⁿ = PlasticState(zeros(3, n_total_inner), zeros(3, n_total_inner), zeros(3, n_total_inner), zeros(n_total_inner), zeros(n_total_inner))
    Plasticⁿ⁺¹ = deepcopy(Plasticⁿ)

    return DFT_SEP_States(superelasticⁿ, superelasticⁿ⁺¹, Plasticⁿ, Plasticⁿ⁺¹)
end

