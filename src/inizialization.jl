# Cleans the output folders from files of the precedent computation
function clean_folders(output_dir)

    if output_dir != ""
        dir = pwd()
        cd(output_dir)
        foreach(rm, filter(endswith(".vtp"), readdir()))
        foreach(rm, filter(endswith(".vtu"), readdir()))
        foreach(rm, filter(endswith(".vtk"), readdir()))
        foreach(rm, filter(endswith(".pvd"), readdir()))
        foreach(rm, filter(endswith(".txt"), readdir()))
        cd(dir)
    end   
end 

# Initializes simulation state, pre-allocates memory, and generates initial VTK files.
function setup_state_simulation(conf::BeamsConfiguration,params::SimulationParams, inter::Union{Nothing, Interaction}, beam2beam, output_dir)
 
    # Group beam state variables into a single structure
    state = SimulationState(conf, params, beam2beam)

    # Prepare VTK data for visualization
    vtkdata = VTKDataBeams(conf, output_dir)

    return state, VTKData(vtkdata)

end

# Initializes the simulation state for a beam-only configuration
function initialize_state_simulation!(conf::BeamsConfiguration, state::SimulationState, params::SimulationParams)
    
    # Assemble system matrices and initialize state variables
    assemble!(conf, state, params)
    state.matricesⁿ.K .= state.matricesⁿ⁺¹.K
    state.matricesⁿ.C .= state.matricesⁿ⁺¹.C 
    state.matricesⁿ.M .= state.matricesⁿ⁺¹.M

end
