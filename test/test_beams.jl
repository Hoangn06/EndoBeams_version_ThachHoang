
#----------------------------------------------------
# PACKAGE INITIALIZATION
#----------------------------------------------------

using Pkg
using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, LinearSolve, Test, Revise, EndoBeams, BenchmarkTools

#----------------------------------------------------# BEAM DEFINITIONS
#----------------------------------------------------

# Define node positions and connectivity for beams
positions = [0.0  0.0  0.0;
0.0  2.5  0.0;
0.0  5.0  0.0;
0.0  7.5  0.0;
0.0 10.0  0.0;
-2.5 10.0  0.0;
-5.0 10.0  0.0;
-7.5 10.0  0.0;
-10.0 10.0  0.0]             
nnodes = size(positions, 1)    # Total number of nodes

nbeams = nnodes - 1                    # Number of beam elements
connectivity = [1 2;
2 3;
3 4;
4 5;
5 6;
6 7;
7 8;
8 9]  # Connectivity for beam elements # Connectivity for beam elements

# Initial conditions for displacements, velocities, accelerations, and rotations
initial_displacements = zeros(size(positions))       # Zero initial displacements
initial_velocities = zeros(size(positions))       # Zero initial velocities
initial_accelerations = zeros(size(positions))       # Zero initial accelerations
initial_rotations = zeros(size(positions))        # Zero initial rotations
initial_angular_velocities = zeros(size(positions))   # Zero initial angular velocities
initial_angular_accelerations = zeros(size(positions))  # Zero initial angular accelerations
plane = "xy"                         # Plane of the problem

# Material and geometric properties for beams
E = 1e6                              # Young's modulus (Pa)
ν = 0.3                              # Poisson's ratio
ρ = 1                                # Density (kg/m³)
radius = 0.3                         # Beam radius (m)
damping = 0                          # Damping coefficient

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
beams = ElasticBeams(nodes, connectivity, E, ν, ρ, radius, damping)

#----------------------------------
# BEAMS CONFIGURATION DEFINITIONS
#----------------------------------

# Displacement-controlled loading: impose Z-displacement at the tip node (node 9)
# DOF = (node - 1) * 6 + component  →  (9-1)*6 + 3 = 51  (global Z-translation)
tip_dof = (nnodes - 1) * 6 + 3
u_max   = 5.0                              # Maximum imposed displacement (m)
T_end   = 10.0                             # Duration of the linear ramp (s)
displacement_function(t) = u_max * t / T_end   # Linear ramp: 0 → u_max over T_end seconds
displacementimposed = ImposedDisplacement([tip_dof], displacement_function)

# No concentrated force (displacement-controlled)
force_function(t, i) = 0.0
concentrated_force = ConcentratedForce(force_function, Int[])

# Degrees of freedom (DOFs) definition
ndofs = nnodes * 6                     # Total number of DOFs (6 per node for displacement and rotation)
blocked_dofs = 1:6                     # First 6 DOFs are fixed (encastre at node 1)
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
# Δt = 1 s keeps the problem in the quasi-static regime (Δt >> sqrt(M/K) ≈ 0.016 s)
# so that elastic reaction forces dominate the inertial forces in each step.
initial_timestep = 1.0     # Initial time step size (s)
min_timestep     = 1e-10   # Minimum allowed time step
max_timestep     = 1.0     # Maximum allowed time step
output_timestep  = 1.0     # Time step for output/visualization
simulation_end_time = T_end   # End time matches the ramp duration

# Convergence criteria for the solver
tolerance_residual    = 1e-2   # Residual tolerance (slightly loose: res_norm normalized by elastic reactions)
tolerance_displacement = 1e-3  # Tolerance for displacement increment ΔD
max_iterations = 10            # Maximum Newton iterations per step

# Store solver parameters in a structured Params object
params = SimulationParams(;
α, β, γ, initial_timestep, min_timestep, max_timestep, output_timestep, simulation_end_time, tolerance_residual, tolerance_displacement, max_iterations, output_dir = "test/output3D"
, verbose = true)

#----------------------------------------------------
# START SIMULATION
#----------------------------------------------------

using Plots

# Run the solver with the defined configurations and parameters
run_simulation!(conf, params)

println("Simulation finished")

#----------------------------------------------------
# PLOT TIP Z-DISPLACEMENT vs TIME
#----------------------------------------------------

# Read displacement data (format: Time, DispX1, DispY1, DispZ1, DispX2, ...)
csv_file = joinpath("test", "output3D", "displacement_data.csv")
displacement_data, _ = readdlm(csv_file, ',', header=true)

times = Float64.(displacement_data[:, 1])

# Column for Z-displacement of tip node (node 9): 1 + (node-1)*3 + 3
tip_disp_col = 1 + (nnodes - 1) * 3 + 3
tip_disp_Z   = Float64.(displacement_data[:, tip_disp_col])

# Imposed displacement target (for comparison)
imposed_vals = [displacement_function(t) for t in times]

p = plot(times, [tip_disp_Z imposed_vals],
    xlabel  = "Time (s)",
    ylabel  = "Z-displacement of tip node (m)",
    title   = "Imposed Z-displacement — tip node 9",
    label   = ["FEM result" "Target (imposed)"],
    linewidth = 2,
    linestyle = [:solid :dash],
    grid    = true,
    dpi     = 300,
)
plot_file = joinpath("test", "output3D", "imposed_displacement_tip.png")
savefig(p, plot_file)
println("Plot saved to: $plot_file")
display(p)
