
#----------------------------------------------------
# PACKAGE INITIALIZATION
#----------------------------------------------------

using Pkg
using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, LinearSolve, Test, Revise, EndoBeams, BenchmarkTools
using Plots


# Spring geometry parameters
coil_diameter      = 5.5    # Outer diameter of the coil (mm)
coil_height        = 7.5    # Axial height (pitch) per coil (mm)
total_spring_height = 27.61  # Total free length of the spring (mm)

# Derived parameter: number of coils from total height and pitch per coil
num_coils = total_spring_height / coil_height

# Discretization
elements_per_coil  = 20                        # Beam elements per coil (controls helix smoothness)
num_elements       = round(Int, elements_per_coil * num_coils)
nnodes             = num_elements + 1

# Build helix node positions
R           = coil_diameter / 2               # Coil radius
total_angle = 2π * num_coils                  # Total swept angle (rad)
positions   = zeros(nnodes, 3)
for i in 1:nnodes
    θ = (i - 1) * total_angle / num_elements
    positions[i, :] = [R * cos(θ), R * sin(θ), coil_height * θ / (2π)]
end

nbeams       = num_elements
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
E = 32000                         # Young's modulus (MPa)
ν = 0.3                              # Poisson's ratio
ρ = 6.5e-9*1e6                                # Density (tonnes/mm³)
radius = 0.5                         # Beam radius (mm)
damping = 0.0                          # Damping coefficient
εL = 0.07/sqrt(2/3)                          # Max transformation strain
sigma_S_AS = 350                    # Start stress for A→S (MPa)
sigma_F_AS = 375                    # Finish stress for A→S (MPa)
sigma_S_SA = 200                    # Start stress for S→A (MPa)
sigma_F_SA = 10                    # Finish stress for S→A (MPa)
sigma_c_SAS = 350                     # Critical stress for SAS transformation (MPa)
Eᴹ = 12000                         # Martensite Young's modulus (MPa)

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
beams = SuperElasticBeams(nodes, connectivity, E, ν, ρ, εL, sigma_S_AS, sigma_F_AS, sigma_S_SA, sigma_F_SA, sigma_c_SAS, Eᴹ, radius, damping)
#----------------------------------
# BEAMS CONFIGURATION DEFINITIONS
#----------------------------------

# Degrees of freedom (DOFs) definition
ndofs = nnodes * 6                     # Total number of DOFs (6 per node for displacement and rotation)
blocked_dofs = [1:6;
(nnodes-1)*6 + 1;
(nnodes-1)*6 + 2
]                      # First 6 DOFs are fixed (e.g., at the boundary)
encastre = Encastre(blocked_dofs)       # Encastre structured

# Displacement imposed
imposed_dofs = Int[]
displacement_applied(t) = 0.0
displacementimposed = ImposedDisplacement(imposed_dofs,displacement_applied)

# External force
loaded_dofs = [(nnodes-1)*6 + 3]   # Last node: uz DOF
force_function(t,i) = t <= 30.0 ? 40.0 * t / 30.0 : 40.0 * (60.0 - t) / 30.0   # 0→30N (t=0→30s), 30→0N (t=30→60s)
concentrated_force = ConcentratedForce(force_function, loaded_dofs)  

# Beam configuration struct initialization
conf = BeamsConfiguration(nodes, beams, Loads(concentrated_force), BoundaryConditions(encastre,displacementimposed, ndofs))  # Store the configuration for beams

#----------------------------------------------------
# SOLVER DEFINITIONS
#----------------------------------------------------

# HHT (Houbolt-Hughes-Taylor) time stepping parameters
α = -0.05   # Typically between 0 and -1, used for numerical stability
β = 0.25 * (1 - α)^2  # Damping parameter based on α
γ = 0.5 * (1 - 2 * α)  # Time-stepping parameter

# General time stepping parameters
initial_timestep = 1e-2  # Initial time step size
min_timestep = 1e-10    # Minimum allowed time step
max_timestep = 1e-2   # Maximum allowed time step (could be adjusted based on system behavior)
output_timestep = 0.5   # Time step for output plotting or visualization
simulation_end_time = 0 # End time for the simulation (loading 0→30s, unloading 30→60s)

# Convergence criteria for the solver
tolerance_residual = 1e-4   # Residual tolerance for convergence checks
tolerance_displacement = 1e-4    # Tolerance for changes in displacement (ΔD)
max_iterations = 20      # Maximum number of iterations for the solver

# Store solver parameters in a structured Params object
params = SimulationParams(;
α, β, γ, initial_timestep, min_timestep, max_timestep, output_timestep, simulation_end_time, tolerance_residual, tolerance_displacement, max_iterations, output_dir = "validation_spring/output3D"
, verbose = true)

#----------------------------------------------------
# START SIMULATION
#----------------------------------------------------

# Run the solver with the defined configurations and parameters
run_simulation!(conf, params)

println("Simulation finished")

#----------------------------------------------------
# PLOT FORCE VS DISPLACEMENT
#----------------------------------------------------

# Read displacement data from CSV file
csv_file = joinpath("validation_spring", "output3D", "displacement_data.csv")
displacement_data, header = readdlm(csv_file, ',', header=true)

times = Float64.(displacement_data[:, 1])

# Top node uz column: 1 (time) + (nnodes-1)*3 columns for previous nodes + 3 (uz = 3rd component)
uz_col = 1 + (nnodes - 1) * 3 + 3
top_uz = Float64.(displacement_data[:, uz_col])

# Reconstruct applied force at each output time
forces = [t <= 30.0 ? 30.0 * t / 30.0 : 30.0 * (60.0 - t) / 30.0 for t in times]

# Plot Force vs. axial displacement of the top node
p = plot(top_uz, forces,
    xlabel    = "Axial displacement uz of top node (mm)",
    ylabel    = "Applied force (N)",
    title     = "Spring tensile test — Force vs. Displacement",
    label     = "FEM",
    linewidth=2,
    grid=true,
    legend=false,
    dpi=300,
)

plot_file = joinpath("validation_spring", "output3D", "force_vs_displacement.png")
savefig(p, plot_file)
println("Plot saved to: $plot_file")
display(p)











