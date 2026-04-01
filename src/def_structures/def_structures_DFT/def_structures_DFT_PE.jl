#----------------------------------
# Structures to store the state of the DFT beam
# (Composite between plastic outer layer + elastic inner layer)
#----------------------------------

struct DFT_PE_States
    Plasticⁿ::PlasticState             # Plastic state at the current time step
    Plasticⁿ⁺¹::PlasticState           # Plastic state at the next time step
end

#----------------------------------
# Functions to initialize the DFT structures
#----------------------------------

function DFT_PE_States(params::SimulationParams)
    nG = params.nᴳ_beams
    
    n_outer = (params.Nr+1) * (params.Nθ+1)
    n_total_outer = nG * n_outer

    Plasticⁿ   = PlasticState(zeros(3, n_total_outer), zeros(3, n_total_outer), zeros(3, n_total_outer), zeros(n_total_outer), zeros(n_total_outer))
    Plasticⁿ⁺¹ = deepcopy(Plasticⁿ)

    return DFT_PE_States(Plasticⁿ, Plasticⁿ⁺¹)
end

