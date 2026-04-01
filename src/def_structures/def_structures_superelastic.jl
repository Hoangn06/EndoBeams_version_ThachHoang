#----------------------------------
# Structures to store the state of the superelasticit
#----------------------------------
struct SuperelasticState
    xi_S::Vector{Float64}
    F::Vector{Float64}
    H_AS::Vector{Float64}
    H_SA::Vector{Float64}
    H_SS::Vector{Float64}
    eps::Matrix{Float64}
    sigma::Matrix{Float64}
end

struct SuperelasticStates
    superelasticⁿ::SuperelasticState # Superelastic state at the current time step
    superelasticⁿ⁺¹::SuperelasticState # Superelastic state at the next time step
end


#----------------------------------
# Functions to initialize the Superelastic structures
#----------------------------------

function SuperelasticState(params::SimulationParams)
    nG = params.nᴳ_beams
    nS = params.nˢ_beams * params.nˢ_beams
    n_total = nG * nS

    return SuperelasticState(zeros(n_total), zeros(n_total), zeros(n_total), zeros(n_total), zeros(n_total), zeros(3, n_total), zeros(3, n_total))
end

function SuperelasticStates(params::SimulationParams)
    superelasticⁿ = SuperelasticState(params)
    superelasticⁿ⁺¹ = deepcopy(superelasticⁿ)
    return SuperelasticStates(superelasticⁿ, superelasticⁿ⁺¹)
end

