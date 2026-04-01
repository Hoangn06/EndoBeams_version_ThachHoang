#----------------------------------
# Create structures used for plasticity computation
#----------------------------------
# Structure to store tensors in Voigt form (6 components) at each Gauss point
struct PlasticState
    sigma::Matrix{Float64} 
    epsilon::Matrix{Float64}
    epsilon_pl::Matrix{Float64}
    e_pl::Vector{Float64}
    F::Vector{Float64}
end

# Structure to store the state of the Plastic parameters
struct PlasticStates
    Plasticⁿ::PlasticState # Plastic parameters at the current time step
    Plasticⁿ⁺¹::PlasticState # Plastic parameters at the next time step
end


#----------------------------------
# Functions to initialize tensors used in the plasticity computation
#----------------------------------
# Constructor for Tensorsplastic structure, initializing parameters to zero
function PlasticState(params::SimulationParams)
    nG = params.nᴳ_beams
    nS = params.nˢ_beams * params.nˢ_beams
    n_total = nG * nS
    return PlasticState(zeros(3, n_total), zeros(3, n_total), zeros(3, n_total), zeros(n_total), zeros(n_total))
end


# Constructor for PlasticStates structure using SimulationParams
function PlasticStates(params::SimulationParams)
    Plasticⁿ = PlasticState(params)
    Plasticⁿ⁺¹ = deepcopy(Plasticⁿ)

    return PlasticStates(Plasticⁿ, Plasticⁿ⁺¹)
end


