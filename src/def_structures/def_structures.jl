# Structure for force vectors required at each time step in the solver
struct Forces
    fᵉˣᵗ::Vector{Float64}    # External forces vector
    Tⁱⁿᵗ::Vector{Float64}    # Internal forces vector
    Tᵏ::Vector{Float64}      # Kinetic forces vector (inertia effects)
    Tᶜ::Vector{Float64}      # Contact forces vector
end  

# Constructor for Solution structure, initializing force vectors based on beams configuration
function Forces(conf::BeamsConfiguration)

    ndofs = conf.ndofs  # Total dofs

    # Initialize external forces vector with zeros
    fᵉˣᵗ = zeros(ndofs)

    # Populate external forces vector based on loaded dofs in configuration
    if conf.loads !== nothing && conf.loads.concentrated_force !== nothing
        for i in conf.loads.concentrated_force.loaded_dofs
            fᵉˣᵗ[i] = conf.loads.concentrated_force.force_function(0,i)
        end
    end 
    
    # Initialize remaining force vectors with zeros
    Tⁱⁿᵗ = zeros(ndofs)  # Internal forces
    Tᵏ = zeros(ndofs)     # Kinetic forces
    Tᶜ = zeros(ndofs)     # Contact forces

    return Forces(fᵉˣᵗ, Tⁱⁿᵗ, Tᵏ, Tᶜ)
end

# Structure to hold nodal solutions with preallocated vectors for solver computations
struct Solution 
    D::Vector{Float64}              # Displacement vector
    Ḋ::Vector{Float64}             # Velocity vector
    D̈::Vector{Float64}             # Acceleration vector
    Ḋⁿ::Vector{Float64}             # Old velocity vector
    r::Vector{Float64}              # Residual forces vector
    Ktan::SparseMatrixCSC{Float64,Int} # Global stiffness matrix in sparse format
    ΔD::Vector{Float64}             # Displacement increment vector
    temp::Vector{Float64}           # Temporary vector for intermediate computations
    r_free::Vector{Float64}         # Free dofs residual vector
    Ktan_free::SparseMatrixCSC{Float64,Int} # Free DOF stiffness matrix
    ΔD_free::Vector{Float64}        # Free DOF displacement increment vector
end 

# Constructor for Solution structure, setting up required vector and matrix allocations
function Solution(Ktan, Ktan_free, ndofs, nfreedofs)
    return Solution(
        zeros(ndofs),          # Displacement vector
        zeros(ndofs),          # Velocity vector
        zeros(ndofs),          # Acceleration vector
        zeros(ndofs),          # Old velocity vector
        zeros(ndofs),          # Residual forces vector
        Ktan,                  # Global stiffness matrix (sparse)
        zeros(ndofs),          # Displacement increment
        zeros(ndofs),          # Temporary vector
        zeros(nfreedofs),      # Free DOF residual forces vector
        Ktan_free,             # Free DOF stiffness matrix
        zeros(nfreedofs)       # Free DOF displacement increment
    )
end

# Structure to hold global matrices (stiffness, damping, mass) for solver operations
struct Matrices
    K::SparseMatrixCSC{Float64,Int} # Global stiffness matrix (sparse)
    C::SparseMatrixCSC{Float64,Int} # Damping matrix (sparse)
    M::SparseMatrixCSC{Float64,Int} # Mass matrix (sparse)
    sparsity_free::Vector{Int}          # Free dofs sparsity pattern
end 

# Constructor for Matrices structure, initializing matrices and forces for solver computations
function Matrices(I, J, sparsity_free)

    # Create sparse matrices for stiffness, damping, and mass
    K = sparse(I, J, 0.)
    C = sparse(I, J, 0.)
    M = sparse(I, J, 0.)
    
    return Matrices(K, C, M, sparsity_free)
end

# Structure to store energy contributions in the simulation
mutable struct Energy
    strain_energy::Float64      # Strain energy contribution
    kinetic_energy::Float64     # Kinetic energy contribution
    contact_energy::Float64     # Contact energy contribution
end

# Constructor for Energy structure, initializing energy values to zero
function Energy()
    return Energy(0.0, 0.0, 0.0)
end

# Union of all per-beam constitutive state objects (linear elastic beams have no entry in the dict).
const NonlinearMaterialState = Union{
    SuperelasticStates, PlasticStates,
    DFT_SEP_States, DFT_EP_States, DFT_PSE_States, DFT_PP_States, DFT_PE_States,
}

struct BeamMaterialStates
    by_beam::Dict{Int, NonlinearMaterialState}
end


struct SimulationState
    forcesⁿ::Forces
    forcesⁿ⁺¹::Forces
    matricesⁿ::Matrices
    matricesⁿ⁺¹::Matrices
    solⁿ⁺¹::Solution
    energyⁿ⁺¹::Energy
    material_states::BeamMaterialStates
    beam2beam_contactsⁿ⁺¹::Union{Nothing, Any}  # Beam-to-beam contact tracking for the current time step
    beam2beam_contactsⁿ::Union{Nothing, Any}     # Beam-to-beam contact tracking for the previous time step
    b2b_sparsity::Union{Nothing, Any}            # Pre-computed sparsity maps for all beam-to-beam pairs
end

function SimulationState(conf::BeamsConfiguration, params::SimulationParams, beam2beam)
    
    ndofs = conf.ndofs
    free_dofs = conf.bcs.free_dofs

    # Allocate state variables for beams
    forcesⁿ = Forces(conf)
    forcesⁿ⁺¹ = deepcopy(forcesⁿ)
    energyⁿ⁺¹ = Energy()
    material_states = BeamMaterialStates(conf, params)

    if beam2beam
        # Build sparse matrices including all possible beam-to-beam contact entries
        matricesⁿ, solutionⁿ⁺¹, b2b_sparsity = sparse_matrices_beams_b2b!(conf)
        beam2beam_contactsⁿ⁺¹ = Beam2BeamContacts()
        beam2beam_contactsⁿ    = deepcopy(beam2beam_contactsⁿ⁺¹)

    else
        matricesⁿ, solutionⁿ⁺¹ = sparse_matrices_beams!(conf)
        b2b_sparsity       = nothing
        beam2beam_contactsⁿ⁺¹ = nothing
        beam2beam_contactsⁿ    = nothing
    end

    matricesⁿ⁺¹ = deepcopy(matricesⁿ)

    return SimulationState(forcesⁿ, forcesⁿ⁺¹, matricesⁿ, matricesⁿ⁺¹, solutionⁿ⁺¹, energyⁿ⁺¹, material_states, beam2beam_contactsⁿ⁺¹, beam2beam_contactsⁿ, b2b_sparsity)
end

# Function to create the appropriate material state for a given beam
function material_state_for_beam(material::Symbol, params::SimulationParams)
    material == :superelastic && return SuperelasticStates(params)
    material == :plastic && return PlasticStates(params)
    material == :DFT_SEP && return DFT_SEP_States(params)
    material == :DFT_EP && return DFT_EP_States(params)
    material == :DFT_PSE && return DFT_PSE_States(params)
    material == :DFT_PP && return DFT_PP_States(params)
    material == :DFT_PE && return DFT_PE_States(params)
    return nothing
end

# Function to create state for all beams
function BeamMaterialStates(conf::BeamsConfiguration, params::SimulationParams)
    d = Dict{Int, NonlinearMaterialState}()
    for beam in LazyRows(conf.beams)
        st = material_state_for_beam(beam.material, params)
        if st !== nothing
            d[beam.ind] = st
        end
    end
    return BeamMaterialStates(d)
end
