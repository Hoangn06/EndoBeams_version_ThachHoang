
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
0.0  10.0  0.0;
0.0  20.0  0.0]   

nnodes = size(positions, 1)    # Total number of nodes

nbeams = nnodes - 1                    # Number of beam elements
connectivity = [1 2;
                2 3]  # Connectivity for beam elements

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
ρ = 6.5e-9*1e6                                # Density (tonnes/mm³)
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
tip_dof = (nnodes - 1) * 6 + 2                                      # Y-DOF of the last node
u_max = 2.5                                                      # Maximum imposed displacement (mm)
displacement_function(t) = t <= 2.0 ? (u_max / 2.0) * t : (t <= 4.0 ? u_max * (1.0 - (t - 2.0) / 2.0) : u_max)
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
min_timestep = 1e-10   # Minimum allowed time step
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

# Read displacement data from CSV file
# Format: Time, DispX1, DispY1, DispZ1, DispX2, DispY2, DispZ2, ...
csv_file = joinpath("test", "output3D", "displacement_data.csv")
displacement_data, header = readdlm(csv_file, ',', header=true)

# Extract time and displacement columns (CSV: Time, then DispXk,DispYk,DispZk per node)
# Tip matches tip_dof / displacement_function: last node, Y component
tip_node = nnodes
dispy_col = 3 + 3 * (tip_node - 1)  # DispY{tip_node} column index (after Time)
times = Float64.(displacement_data[:, 1])
displacements = Float64.(displacement_data[:, dispy_col])

# Displacement vs time (tip Y)
p_disp_time = plot(times, displacements,
    xlabel="Time (s)",
    ylabel="Displacement (mm)",
    title="Tip Y displacement vs time",
    linewidth=2,
    grid=true,
    legend=false,
    dpi=300,
)
disp_time_plot_file = joinpath("test", "output3D", "displacement_vs_time.png")
savefig(p_disp_time, disp_time_plot_file)
println("Plot saved to: $disp_time_plot_file")
display(p_disp_time)


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
    xticks=0:0.05:1.0,
    yticks=0:100:1000,
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






