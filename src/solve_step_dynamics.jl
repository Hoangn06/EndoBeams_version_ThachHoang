function solve_step_dynamics!(conf::SimulationConfiguration, state::SimulationState, params::SimulationParams, inter::Union{Nothing, Interaction},beam2beam, tⁿ⁺¹, Δt, solver) 
    
    # Update external loads at the current time step
    update_loads!(conf, state, tⁿ⁺¹)
    update_boundary_conditions!(conf, tⁿ⁺¹)

    # Update moving surfaces (nothing if fixed surfaces or no interaction)
    if isa(inter, MultiRigidInteraction)
        for ri in inter.interactions
            update_surface!(ri.master, tⁿ⁺¹)
        end
    elseif isa(inter, Interaction)
        update_surface!(inter.master, tⁿ⁺¹)
    end

    # Predictor phase: Estimate the solution for the current step
    @timeit_debug "Predictor" predictor!(conf, state, params, Δt, solver)

    # Corrector phase: Iteratively refine the solution until convergence
    @timeit_debug "Corrector" k = corrector!(conf, state, params, inter, beam2beam, Δt, solver)

    return k
end
