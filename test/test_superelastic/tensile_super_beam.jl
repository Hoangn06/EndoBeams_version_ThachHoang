
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
positions = [0.0  0.0  0.0;
0.0  20.0  0.0]           
nnodes = size(positions, 1)    # Total number of nodes

nbeams = nnodes - 1                    # Number of beam elements
connectivity = [1 2]  # Connectivity for beam elements

# Initial conditions for displacements, velocities, accelerations, and rotations
initial_displacements = zeros(size(positions))       # Zero initial displacements
initial_velocities = zeros(size(positions))       # Zero initial velocities
initial_accelerations = zeros(size(positions))       # Zero initial accelerations
initial_rotations = zeros(size(positions))        # Zero initial rotations
initial_angular_velocities = zeros(size(positions))   # Zero initial angular velocities
initial_angular_accelerations = zeros(size(positions))  # Zero initial angular accelerations
plane = "xy"                         # Plane of the problem

# Material and geometric properties for beams
E = 60000                         # Young's modulus (MPa)
ν = 0.3                              # Poisson's ratio
ρ = 6.5e-9                                # Density (tonnes/mm³)
radius = 0.745                         # Beam radius (mm)
damping = 0.0                          # Damping coefficient
εL = 0.075/sqrt(2/3)                          # Max transformation strain
sigma_S_AS = 520                    # Start stress for A→S (MPa)
sigma_F_AS = 600                    # Finish stress for A→S (MPa)
sigma_S_SA = 300                    # Start stress for S→A (MPa)
sigma_F_SA = 200                    # Finish stress for S→A (MPa)
sigma_c_SAS = 520                     # Critical stress for SAS transformation (MPa)
Eᴹ = 60000                         # Martensite Young's modulus (MPa)

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

# Force-controlled loading on DOF 8 (uy of node 2 — axial direction)
# Ramps from 0 to 3 N at t=2, then back to 0 at t=4
force_function(t, i) = t <= 2.0 ? 1000 * t : (t <= 4.0 ? 1000 * (4.0 - t) : 0.0)
concentrated_force = ConcentratedForce(force_function, [8])

# Degrees of freedom (DOFs) definition
ndofs = nnodes * 6                     # Total number of DOFs (6 per node for displacement and rotation)
blocked_dofs = [1:6;]                       # First 6 DOFs are fixed (e.g., at the boundary)
encastre = Encastre(blocked_dofs)

# Beam configuration struct initialization
conf = BeamsConfiguration(nodes, beams, Loads(concentrated_force), BoundaryConditions(encastre, ndofs))

#----------------------------------------------------
# SOLVER DEFINITIONS
#----------------------------------------------------

# HHT (Houbolt-Hughes-Taylor) time stepping parameters
α = -0.05   # Typically between 0 and -1, used for numerical stability
β = 0.25 * (1 - α)^2  # Damping parameter based on α
γ = 0.5 * (1 - 2 * α)  # Time-stepping parameter

# General time stepping parameters
initial_timestep = 1e-4    # Initial time step size
min_timestep = 1e-10   # Minimum allowed time step
max_timestep = 1e-1    # Maximum allowed time step (could be adjusted based on system behavior)
output_timestep = 1e-1   # Time step for output plotting or visualization
simulation_end_time = 4 # End time for the simulation (duration of the analysis)

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

# Read displacement data from CSV file
# Format: Time, DispX1, DispY1, DispZ1, DispX2, DispY2, DispZ2, ...
csv_file = joinpath("test", "output3D", "displacement_data.csv")
displacement_data, header = readdlm(csv_file, ',', header=true)

# Extract time and displacement columns
# For node 5 (tip): DispX5=col14, DispY5=col15, DispZ5=col16
times = Float64.(displacement_data[:, 1])  # Time column
displacements = Float64.(displacement_data[:, 6])  # DispY5 (Y displacement of node 5 / tip)

# Calculate force: piecewise function
# force = 1.5 * time if time <= 2, else force = 1.5 * (4.0 - time)
forces = [t <= 2 ? 1.5 * t : 1.5 * (4.0 - t) for t in times]

# Create force vs time plot
p_force_time = plot(times, forces,
    xlabel="Time (s)",
    ylabel="Force (N)",
    title="Force vs Time",
    linewidth=2,
    grid=true,
    legend=false,
    dpi=300,
)

# Save the force vs time plot
force_time_plot_file = joinpath("test", "output3D", "force_vs_time.png")
savefig(p_force_time, force_time_plot_file)
println("Plot saved to: $force_time_plot_file")

# Display the plot
display(p_force_time)

# Create force vs displacement plot
p = plot(displacements, forces,
    xlabel="Displacement (mm)",
    ylabel="Force (N)",
    title="Force vs Displacement at tip of the beam",
    linewidth=2,
    grid=true,
    legend=false,
    dpi=300,
)

# Save the plot
plot_file = joinpath("test", "output3D", "force_vs_displacement.png")
savefig(p, plot_file)
println("Plot saved to: $plot_file")

# Display the plot
display(p)


#----------------------------------------------------
# PLOT STRESS VS STRAIN
#----------------------------------------------------

# Read stress/strain data from CSV file
stress_strain_csv_file = joinpath("test", "output3D", "stress_strain_data.csv")
stress_strain_data, stress_strain_header = readdlm(stress_strain_csv_file, ',', header=true)

# CSV format: Time, stress_GP1, strain_GP1, stress_GP2, strain_GP2, ..., stress_GP27, strain_GP27
# Column indices: 1=Time, 2=stress_GP1, 3=strain_GP1, 4=stress_GP2, 5=strain_GP2, etc.

# Create stress vs strain plot with all 27 Gauss points
n_gauss_points = 27

# Generate distinct colors for each Gauss point using a color gradient
using ColorSchemes
colors = get(ColorSchemes.turbo, range(0, 1, length=n_gauss_points))

# Initialize plot with first Gauss point
strain_gp1 = Float64.(stress_strain_data[:, 3])  # strain_GP1 (column 3)
stress_gp1 = Float64.(stress_strain_data[:, 2])  # stress_GP1 (column 2)

p_stress_strain = plot(strain_gp1, stress_gp1,
    xlabel="Von-Mises Strain",
    ylabel="Von-Mises Stress (MPa)",
    title="              Von-Mises Stress vs. Von-Mises Strain (All 27 Gauss Points)",
    titlefontsize=30,
    xguidefontsize=28,
    yguidefontsize=28,
    xtickfontsize=24,
    ytickfontsize=24,
    linewidth=4,
    grid=true,
    label="GP 1",
    color=colors[1],
    size=(1772, 1181),  # 15cm x 10cm at 300dpi (no dpi parameter needed)
    left_margin=15Plots.mm,   # Space for Y-axis label
    bottom_margin=20Plots.mm, # Space for caption
    legend=:outerright,
    legendfontsize=19
)

# Add remaining Gauss points to the plot with distinct colors
for gp in 2:n_gauss_points
    stress_col = 2 + 2*(gp-1)  # stress column index
    strain_col = 3 + 2*(gp-1)  # strain column index
    
    strain_gp = Float64.(stress_strain_data[:, strain_col])
    stress_gp = Float64.(stress_strain_data[:, stress_col])
    
    plot!(p_stress_strain, strain_gp, stress_gp,
        linewidth=4,
        label="GP $gp",
        color=colors[gp]
    )
end

# Add figure caption below xlabel (compact version)
plot!(p_stress_strain, 
    xlabel="Von-Mises Strain\n──────────────────────────────────────────────────────────── \n Fig.1. Von-Mises Stress vs. Von-Mises Strain during tensile simulation")

# Save the plot
stress_strain_plot_file = joinpath("test", "output3D", "stress_vs_strain.png")
savefig(p_stress_strain, stress_strain_plot_file)
println("Stress vs Strain plot saved to: $stress_strain_plot_file")

# Display the plot
display(p_stress_strain)






