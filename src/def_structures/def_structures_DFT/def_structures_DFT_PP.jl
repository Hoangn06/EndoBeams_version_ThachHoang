#----------------------------------
# Structures to store the state of the DFT beam
# (Composite between plastic outer fibers + plastic inner fibers)

#----------------------------------
struct DFT_PP_States
    Plasticⁿ_outer::PlasticState             # Plastic state at the current time step
    Plasticⁿ⁺¹_outer::PlasticState           # Plastic state at the next time step
    Plasticⁿ_inner::PlasticState             # Plastic state at the current time step
    Plasticⁿ⁺¹_inner::PlasticState           # Plastic state at the next time step
end

#----------------------------------
# Functions to initialize the DFT structures
#----------------------------------

function DFT_PP_States(params::SimulationParams)
    nG = params.nᴳ_beams
    
    n_inner = params.nˢ_beams * params.nˢ_beams
    n_total_inner = nG * n_inner

    n_outer = (params.Nr+1) * (params.Nθ+1)
    n_total_outer = nG * n_outer

    Plasticⁿ_outer = PlasticState(zeros(3, n_total_outer), zeros(3, n_total_outer), zeros(3, n_total_outer), zeros(n_total_outer), zeros(n_total_outer))
    Plasticⁿ⁺¹_outer = deepcopy(Plasticⁿ_outer)

    Plasticⁿ_inner = PlasticState(zeros(3, n_total_inner), zeros(3, n_total_inner), zeros(3, n_total_inner), zeros(n_total_inner), zeros(n_total_inner))
    Plasticⁿ⁺¹_inner = deepcopy(Plasticⁿ_inner)

    return DFT_PP_States(Plasticⁿ_outer, Plasticⁿ⁺¹_outer, Plasticⁿ_inner, Plasticⁿ⁺¹_inner)
end

