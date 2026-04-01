#----------------------------------
# Structures to store the state of the DFT beam
# (Composite between superelastic outer tube + plastic inner fibers)

#----------------------------------
struct DFT_EP_States
    Plasticⁿ::PlasticState             # Plastic state at the current time step
    Plasticⁿ⁺¹::PlasticState           # Plastic state at the next time step
end

#----------------------------------
# Functions to initialize the DFT structures
#----------------------------------

function DFT_EP_States(params::SimulationParams)
    nG = params.nᴳ_beams
    
    n_inner = params.nˢ_beams * params.nˢ_beams
    n_total_inner = nG * n_inner

    Plasticⁿ = PlasticState(zeros(3, n_total_inner), zeros(3, n_total_inner), zeros(3, n_total_inner), zeros(n_total_inner), zeros(n_total_inner))
    Plasticⁿ⁺¹ = deepcopy(Plasticⁿ)

    return DFT_EP_States(Plasticⁿ, Plasticⁿ⁺¹)
end

