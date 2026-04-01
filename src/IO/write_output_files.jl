using DelimitedFiles

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ROTATION DATA STORAGE AND WRITING
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Structure to store rotation data during simulation for all nodes
mutable struct RotationStorage
    times::Vector{Float64}
    rotations::Vector{Vector{Float64}}  # Each entry contains 
    n_nodes::Int  # Node number
    
    function RotationStorage()
        new(Float64[], Vector{Float64}[], 0)
    end
end

# Function to accumulate rotation data when simulation converges
function accumulate_rotation_data!(storage::RotationStorage, conf::BeamsConfiguration, t::Float64)
    n_nodes = length(conf.nodes)
    if n_nodes >= 1
        # Set number of nodes on first call
        if storage.n_nodes == 0
            storage.n_nodes = n_nodes
        end
        
        push!(storage.times, t)
        
        all_rotations = Float64[]
        for node_idx in 1:n_nodes
            # Extract total rotation angle from rotation matrix R
            # toangle converts rotation matrix to axis-angle representation
            total_rotation = toangle(conf.nodes.R[node_idx])
            push!(all_rotations, total_rotation[1])  # X-component
            push!(all_rotations, total_rotation[2])  # Y-component
            push!(all_rotations, total_rotation[3])  # Z-component
        end
        push!(storage.rotations, all_rotations)
    end
end

# Function to write rotation data to CSV file (Excel-compatible)
function write_rotation_data(storage::RotationStorage, output_file::String="rotation_data.csv")
    n_timesteps = length(storage.times)
    n_nodes = storage.n_nodes
    
    if n_timesteps == 0
        return
    end
    
    # Create data matrix: 1 time column + 3*n_nodes rotation columns
    data = zeros(n_timesteps, 1 + 3 * n_nodes)
    
    for i in 1:n_timesteps
        data[i, 1] = storage.times[i]
        for j in 1:(3 * n_nodes)
            data[i, 1 + j] = storage.rotations[i][j]
        end
    end
    
    # Write to CSV file with headers
    open(output_file, "w") do io
        header = "Time"
        for node in 1:n_nodes
            header *= ",RotX$node,RotY$node,RotZ$node"
        end
        println(io, header)
        # Write data
        writedlm(io, data, ',')
    end
end

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# DISPLACEMENT DATA STORAGE AND WRITING
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Structure to store displacement data during simulation for all nodes
mutable struct DisplacementStorage
    times::Vector{Float64}
    displacements::Vector{Vector{Float64}}  # Each entry contains [ux_node1, uy_node1, uz_node1, ux_node2, uy_node2, uz_node2, ...]
    n_nodes::Int  # Number of nodes (set on first accumulation)
    
    function DisplacementStorage()
        new(Float64[], Vector{Float64}[], 0)
    end
end

# Function to accumulate displacement data when simulation converges
# Stores displacement (X, Y, Z) for all nodes
function accumulate_displacement_data!(storage::DisplacementStorage, conf::BeamsConfiguration, t::Float64)
    n_nodes = length(conf.nodes)
    if n_nodes >= 1
        # Set number of nodes on first call
        if storage.n_nodes == 0
            storage.n_nodes = n_nodes
        end
        
        push!(storage.times, t)
        
        all_displacements = Float64[]
        for node_idx in 1:n_nodes
            # Extract displacement from node
            u = conf.nodes.u[node_idx]
            push!(all_displacements, u[1])  # X-component
            push!(all_displacements, u[2])  # Y-component
            push!(all_displacements, u[3])  # Z-component
        end
        push!(storage.displacements, all_displacements)
    end
end

# Function to write displacement data to CSV file (Excel-compatible)
function write_displacement_data(storage::DisplacementStorage, output_file::String="displacement_data.csv")
    n_timesteps = length(storage.times)
    n_nodes = storage.n_nodes
    
    if n_timesteps == 0
        return
    end
    
    # Create data matrix: 1 time column + 3*n_nodes displacement columns
    data = zeros(n_timesteps, 1 + 3 * n_nodes)
    
    for i in 1:n_timesteps
        data[i, 1] = storage.times[i]
        for j in 1:(3 * n_nodes)
            data[i, 1 + j] = storage.displacements[i][j]
        end
    end
    
    # Write to CSV file with headers
    open(output_file, "w") do io
        header = "Time"
        for node in 1:n_nodes
            header *= ",DispX$node,DispY$node,DispZ$node"
        end
        println(io, header)
        # Write data
        writedlm(io, data, ',')
    end
end

# Structure to store stress and strain VM data during simulation for all Gauss points
mutable struct StressStrainStorage
    times::Vector{Float64}
    stress_VM_values::Vector{Vector{Float64}}  # Each entry is a vector of n_gauss_points values
    strain_VM_values::Vector{Vector{Float64}}  # Each entry is a vector of n_gauss_points values
    
    function StressStrainStorage()
        new(Float64[], Vector{Float64}[], Vector{Float64}[])
    end
end

# Function to accumulate stress and strain VM data when simulation converges
# Stores stress_VM and strain_VM for all Gauss points
function accumulate_stress_strain_data!(storage::StressStrainStorage, state::SimulationState, t::Float64, stress_VM::Vector{Float64}, strain_VM::Vector{Float64})
    push!(storage.times, t)
    push!(storage.stress_VM_values, copy(stress_VM))
    push!(storage.strain_VM_values, copy(strain_VM))
end

# Function to write stress and strain VM data to CSV file (Excel-compatible)
# Format: Time, stress_GP1, strain_GP1, stress_GP2, strain_GP2, ..., stress_GPn, strain_GPn
function write_stress_strain_data(storage::StressStrainStorage, output_file::String="stress_strain_data.csv")
    n_timesteps = length(storage.times)
    
    # Check if there's any data to write
    if n_timesteps == 0 || isempty(storage.stress_VM_values)
        @warn "No stress-strain data to write. Skipping stress_strain_data.csv output."
        return
    end
    
    n_gauss_points = length(storage.stress_VM_values[1])  # Get number of Gauss points from first entry
    
    # Create data matrix: 1 time column + n_gauss_points*2 stress/strain columns
    data = zeros(n_timesteps, 1 + 2 * n_gauss_points)
    
    for i in 1:n_timesteps
        data[i, 1] = storage.times[i]
        for gp in 1:n_gauss_points
            data[i, 2 + 2*(gp-1)] = storage.stress_VM_values[i][gp]      # stress column
            data[i, 3 + 2*(gp-1)] = storage.strain_VM_values[i][gp]      # strain column
        end
    end
    
    # Write to CSV file with headers
    open(output_file, "w") do io
        # Build header: Time, stress_GP1, strain_GP1, stress_GP2, strain_GP2, ...
        header = "Time"
        for gp in 1:n_gauss_points
            header *= ",stress_GP$gp,strain_GP$gp"
        end
        println(io, header)
        # Write data
        writedlm(io, data, ',')
    end
end

