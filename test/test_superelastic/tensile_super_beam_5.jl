
#----------------------------------------------------
# PACKAGE INITIALIZATION
#----------------------------------------------------

using Pkg
using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, LinearSolve, Test, Revise, EndoBeams, BenchmarkTools
using Plots

#----------------------------------------------------# BEAM DEFINITIONS
#----------------------------------------------------

# Define node positions and connectivity for beams
positions = [0.0   0.0   0.0;
             0.0   2.0   0.0;
             0.0   4.0   0.0;
             0.0   6.0   0.0;
             0.0   8.0   0.0;
             0.0  10.0   0.0;
             0.0  12.0   0.0;
             0.0  14.0   0.0;
             0.0  16.0   0.0;
             0.0  18.0   0.0;
             0.0  20.0   0.0]           
nnodes = size(positions, 1)    # Total number of nodes

nbeams = nnodes - 1                    # Number of beam elements
connectivity = [1 2;
        2 3;
        3 4;
        4 5;
        5 6;
        6 7;
        7 8;
        8 9;
        9 10;
        10 11]  # Connectivity for beam elements

 # Initial conditions for displacements, velocities, accelerations, and rotations
initial_displacements = zeros(size(positions))       # Zero initial displacements
initial_velocities = zeros(size(positions))       # Zero initial velocities
initial_accelerations = zeros(size(positions))       # Zero initial accelerations
initial_rotations = zeros(size(positions))        # Zero initial rotations
initial_angular_velocities = zeros(size(positions))   # Zero initial angular velocities
initial_angular_accelerations = zeros(size(positions))  # Zero initial angular accelerations
plane = "xy"                         # Plane of the problem

# Material and geometric properties for beams
# Material and geometric properties for beams
E = 32000                         # Young's modulus (MPa)
ν = 0.33                              # Poisson's ratio
ρ = 6.5e-9                                # Density (tonnes/mm³)
radius = 0.5                         # Beam radius (mm)
damping = 0.0                          # Damping coefficient
εL = 0.07/sqrt(2/3)                          # Max transformation strain
sigma_S_AS = 350                    # Start stress for A→S (MPa)
sigma_F_AS = 375                    # Finish stress for A→S (MPa)
sigma_S_SA = 200                    # Start stress for S→A (MPa)
sigma_F_SA = 10                    # Finish stress for S→A (MPa)
sigma_c_SAS = 350                     # Critical stress for SAS transformation (MPa)
Eᴹ = 12000 

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

# Displacement-controlled loading at the tip node (node 11), Y-direction DOF = (11-1)*6 + 2 = 62
# Ramps from 0 to 0.015 mm at t=2, then back to 0 at t=4
tip_dof = (nnodes - 1) * 6 + 2                                    # Y-DOF of the last node
u_max = 2.5                                                      # Maximum imposed displacement (mm)
displacement_function(t) = t <= 2.0 ? (u_max / 2.0) * t :
                            (t <= 4.0 ? u_max * (1.0 - (t - 2.0) / 2.0) : 0.0)
displacementimposed = ImposedDisplacement([tip_dof], displacement_function)

# No concentrated force (displacement-controlled)
force_function(t, i) = 0.0
concentrated_force = ConcentratedForce(force_function, Int[])

# Degrees of freedom (DOFs) definition
ndofs = nnodes * 6                     # Total number of DOFs (6 per node for displacement and rotation)
blocked_dofs = [1:6;]                  # First 6 DOFs are fixed (encastre at node 1)
encastre = Encastre(blocked_dofs)

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
initial_timestep = 1e-2    # Initial time step size
min_timestep = 1e-6   # Minimum allowed time step
max_timestep = 1e-2    # Maximum allowed time step (could be adjusted based on system behavior)
output_timestep = 1e-1   # Time step for output plotting or visualization
simulation_end_time = 3.99 # End time for the simulation (duration of the analysis)

# Convergence criteria for the solver
tolerance_residual = 1e-5   # Residual tolerance for convergence checks
tolerance_displacement = 1e-5    # Tolerance for changes in displacement (ΔD)
max_iterations = 10      # Maximum number of iterations for the solver

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
# PLOT FORCE VS DISPLACEMENT
#----------------------------------------------------

# Read internal force/moment data from CSV file
# Format: Time, FintX1, FintY1, FintZ1, MintX1, MintY1, MintZ1, FintX2, ...
force_moment_file = joinpath(params.output_dir, "internal_force_moment_data.csv")
force_moment_data, header = readdlm(force_moment_file, ',', header=true)

# Extract time and internal force in Y at the last node
tip_node = nnodes
fint_y_col = 1 + 6 * (tip_node - 1) + 2   # FintY{tip_node}
times = Float64.(force_moment_data[:, 1])
forces = Float64.(force_moment_data[:, fint_y_col])

# Use the imposed displacement law as the X-axis
imposed_displacements = [displacement_function(t) for t in times]

# Save the curve data to CSV
force_displacement_data_file = joinpath(params.output_dir, "force_vs_displacement_data.csv")
open(force_displacement_data_file, "w") do io
    println(io, "Displacement,FintY$tip_node")
    writedlm(io, hcat(imposed_displacements, forces), ',')
end

# Split into loading and unloading phases
loading_idx = times .<= 2.0
unloading_idx = times .> 2.0

# Create force vs displacement plot
p = plot(imposed_displacements[loading_idx], forces[loading_idx],
    xlabel="Imposed displacement (mm)",
    ylabel="Internal force Y at node $tip_node",
    title="Internal Force vs Imposed Displacement",
    linewidth=2,
    grid=true,
    label="Loading",
    color=:blue,
    dpi=300,
)
plot!(p, imposed_displacements[unloading_idx], forces[unloading_idx],
    linewidth=2,
    label="Unloading",
    color=:red,
)

# Save the plot
plot_file = joinpath(params.output_dir, "force_vs_displacement.png")
savefig(p, plot_file)
println("Force-displacement data saved to: $force_displacement_data_file")
println("Plot saved to: $plot_file")


