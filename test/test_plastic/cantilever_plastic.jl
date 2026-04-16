
#----------------------------------------------------
# PACKAGE INITIALIZATION
#----------------------------------------------------

using Pkg
using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, LinearSolve, Test, Revise, EndoBeams, BenchmarkTools
using Plots

# BEAM DEFINITIONS
#----------------------------------------------------
# Validation: cantilever beam
# Diameter D = 1.49 mm
# Length L = 20 mm
# Clamped at node 1 (all 6 DOFs fixed)
# Transverse force F applied at tip (node nnodes) in X direction
# Properties:
# E = 60000 MPa, Eᴹ = 20000 MPa
# v = 0.3
# epsilon_L = 7.5%/sqrt(2/3)
# sigma_S_AS = 500 MPa, sigma_F_AS = 600 MPa
# sigma_S_SA = 300 MPa, sigma_F_SA = 200 MPa


# Define node positions and connectivity for beams
num_elements = 10
nnodes = num_elements + 1    # Total number of nodes
beam_length = 20            # Total length of the beam
positions = zeros(nnodes, 3)  # Initialize position matrix
for i in 1:nnodes
    positions[i, :] = [0.0, (i-1) * beam_length / num_elements, 0.0]
end

nbeams = nnodes - 1                    # Number of beam elements
# Generate connectivity for beam elements
connectivity = zeros(Int, nbeams, 2)
for i in 1:nbeams
    connectivity[i, :] = [i, i+1]
end

# Initial conditions for displacements, velocities, accelerations, and rotations
initial_displacements = zeros(size(positions))       # Zero initial displacements
initial_velocities = zeros(size(positions))       # Zero initial velocities
initial_accelerations = zeros(size(positions))       # Zero initial accelerations
initial_rotations = zeros(size(positions))        # Zero initial rotations
initial_angular_velocities = zeros(size(positions))   # Zero initial angular velocities
initial_angular_accelerations = zeros(size(positions))  # Zero initial angular accelerations
plane = "xy"                         # Plane of the problem

# Material and geometric properties for beams
E = 51000                         # Young's modulus (MPa)
ν = 0.3                              # Poisson's ratio
ρ = 6.5e-9                                # Density (tonnes/mm³)
radius = 0.035                         # Beam radius (mm)
damping = 0                          # Damping coefficient
y0 = 170.0                         # Yield stress (MPa)
K = 20000.0                         # Plastic modulus (MPa)

# Build nodes and beams
nodes = NodesBeams(
positions, 
initial_displacements, 
initial_velocities, 
initial_accelerations, 
initial_rotations, 
initial_angular_velocities, 
initial_angular_accelerations, 
plane
)
beams = PlasticBeams(nodes, connectivity, E, ν, ρ, y0, K, radius, damping)

#----------------------------------
# BEAMS CONFIGURATION DEFINITIONS
#----------------------------------

# External force: global DOF index = (node_index - 1) * 6 + component (1–3 = ux,uy,uz; 4–6 = rotations)
# Beam along +Y → apply force in Y at the tip node (uy DOF).
simulation_end_time = 30.0              # End time (s): loading 0→15s, unloading 15→30s
t_peak              = 15.0              # Time at peak force
loaded_dofs = [(nnodes-1)*6 + 1]        # Tip node: X displacement (transverse force)
force_max = 20.0                        # Maximum force (N)
force_function(t, i) = t <= t_peak ? force_max * t / t_peak : force_max * (simulation_end_time - t) / t_peak   # 0→60N→0N
concentrated_force = ConcentratedForce(force_function, loaded_dofs)  

# Degrees of freedom (DOFs) definition
ndofs = nnodes * 6                     # Total number of DOFs (6 per node for displacement and rotation)
blocked_dofs = 1:6
encastre = Encastre(blocked_dofs)

# Beam configuration struct initialization
conf = BeamsConfiguration(nodes, beams, Loads(concentrated_force), BoundaryConditions(encastre, ndofs))  # Store the configuration for beams

#----------------------------------------------------
# SOLVER DEFINITIONS
#----------------------------------------------------

# HHT (Houbolt-Hughes-Taylor) time stepping parameters
α = -0.05   # Typically between 0 and -1, used for numerical stability
β = 0.25 * (1 - α)^2  # Damping parameter based on α
γ = 0.5 * (1 - 2 * α)  # Time-stepping parameter

# General time stepping parameters
initial_timestep = 1e-2  # Initial time step size
min_timestep = 1e-8    # Minimum allowed time step
max_timestep = 1e-2   # Maximum allowed time step (could be adjusted based on system behavior)
output_timestep = 1E-2   # Time step for output plotting or visualization
simulation_end_time = 30.0 # End time — matches force_function definition above

# Convergence criteria for the solver
tolerance_residual = 1e-5   # Residual tolerance for convergence checks
tolerance_displacement = 1e-5    # Tolerance for changes in displacement (ΔD)
max_iterations = 20      # Maximum number of iterations for the solver

# Store solver parameters in a structured Params object
params = SimulationParams(;
α, β, γ, initial_timestep, min_timestep, max_timestep, output_timestep, simulation_end_time, tolerance_residual, tolerance_displacement, max_iterations, output_dir = "test/output3D"
, verbose = true)

#----------------------------------------------------
# START SIMULATION
#----------------------------------------------------

# Run the solver with the defined configurations and parameters
run_simulation!(conf, params)

println("Simulation finished")

#----------------------------------------------------
# PLOT FORCE VS X-DISPLACEMENT OF NODE 16
#----------------------------------------------------

# Read displacement data from CSV file
# Format: Time, DispX1, DispY1, DispZ1, DispX2, DispY2, DispZ2, ...
csv_file = joinpath("test", "output3D", "displacement_data.csv")
displacement_data, header = readdlm(csv_file, ',', header=true)

# Column for X displacement of tip node
tracked_node = nnodes
disp_x_col = 2 + 3 * (tracked_node - 1)
times         = Float64.(displacement_data[:, 1])
displacements = Float64.(displacement_data[:, disp_x_col])

# Reconstruct the applied force at each output time
forces = [t <= t_peak ? force_max * t / t_peak : force_max * (simulation_end_time - t) / t_peak for t in times]

# Force vs X-displacement plot
p = plot(displacements, forces,
    xlabel="X Displacement of tip node $tracked_node (mm)",
    ylabel="Force (N)",
    title="Cantilever — Force vs Tip X-Displacement",
    linewidth=2,
    grid=true,
    legend=false,
    dpi=300,
)

plot_file = joinpath("test", "output3D", "force_vs_displacement.png")
savefig(p, plot_file)
println("Plot saved to: $plot_file")
display(p)












