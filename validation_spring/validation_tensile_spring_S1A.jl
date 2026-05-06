
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
total_spring_height = 22.87  # Total free length of the spring (mm)

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
ρ = 6.5e-9                                # Density (tonnes/mm³)
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
#beams = ElasticBeams(nodes, connectivity, E, ν, ρ, radius, damping)
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

# Displacement-controlled: axial uz at last node, 0 → u_max mm over 0 → 30 s, then hold
tip_uz_dof = (nnodes - 1) * 6 + 3
u_max = 30.0   # mm
displacement_function(t) = t <= 30.0 ? u_max * t / 30.0 : u_max
displacementimposed = ImposedDisplacement([tip_uz_dof], displacement_function)

# No external concentrated force (displacement-controlled)
force_function(t, i) = 0.0
concentrated_force = ConcentratedForce(force_function, Int[])

# Beam configuration struct initialization
conf = BeamsConfiguration(nodes, beams, Loads(concentrated_force), BoundaryConditions(encastre, displacementimposed, ndofs))

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
simulation_end_time = 0 # End time (loading ramp 0→u_max mm over 0→30 s)

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
# PLOT DISPLACEMENT VS TIME (displacement-controlled tensile)
#----------------------------------------------------

# Read displacement data from CSV file
csv_file = joinpath("validation_spring", "output3D", "displacement_data.csv")
displacement_data, header = readdlm(csv_file, ',', header=true)

times = Float64.(displacement_data[:, 1])

# Top node uz column: 1 (time) + (nnodes-1)*3 columns for previous nodes + 3 (uz = 3rd component)
uz_col = 1 + (nnodes - 1) * 3 + 3
top_uz = Float64.(displacement_data[:, uz_col])

# Imposed uz history (same law as displacement_function)
imposed_uz = [t <= 30.0 ? u_max * t / 30.0 : u_max for t in times]

# Axial displacement vs time: FEM vs imposed BC
p = plot(times, top_uz,
    xlabel    = "Time (s)",
    ylabel    = "Axial displacement uz of top node (mm)",
    title     = "Spring tensile test — displacement vs time",
    label     = "FEM",
    linewidth=2,
    grid=true,
    legend=:topleft,
    dpi=300,
)
plot!(p, times, imposed_uz, label="Imposed", linestyle=:dash, linewidth=2)

plot_file = joinpath("validation_spring", "output3D", "displacement_vs_time.png")
savefig(p, plot_file)
println("Plot saved to: $plot_file")
display(p)











