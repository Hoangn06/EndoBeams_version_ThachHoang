#----------------------------------------------------
# PACKAGE INITIALIZATION
#----------------------------------------------------

using Pkg
using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, LinearSolve, Test, Revise, EndoBeams, BenchmarkTools

#----------------------------------------------------# BEAM DEFINITIONS
#----------------------------------------------------

dx = 5; L = 10
aux = 0:dx:L
n = length(aux)
x = aux .- 4.31

#Generate positions of the beam
X0 = hcat(zeros(Float64, n), collect(x), fill(3.23, n))
X1 = hcat(fill(0.06,n), zeros(Float64, n), collect(aux))
positions  = vcat(X0, X1)
positions  = positions./100

# Define node positions and connectivity for beams
nnodes = size(positions, 1)    # Total number of nodes
nbeams = nnodes - 1     # Number of beam elements

half = Int(nnodes * 0.5)
conn1 = hcat(1:half-1,       2:half)  # build first half connectivity
conn2 = hcat(half+1:nnodes-1, half+2:nnodes) # build second half connectivity
connectivity = vcat(conn1, conn2) # combine into a single 2-column matrix


# Initial conditions for displacements, velocities, accelerations, and rotations
initial_displacements = zeros(size(positions))       # Zero initial displacements
initial_velocities = zeros(size(positions))       # Zero initial velocities
initial_accelerations = zeros(size(positions))       # Zero initial accelerations
initial_rotations = zeros(size(positions))        # Zero initial rotations
initial_angular_velocities = zeros(size(positions))   # Zero initial angular velocities
initial_angular_accelerations = zeros(size(positions))  # Zero initial angular accelerations
plane = "xy"                         # Plane of the problem

# Material and geometric properties for beams
E = 5*1e9                            # Young's modulus (Pa)
ν = 0.33                             # Poisson's ratio
ρ = 7850                             # Density (kg/m³)
radius = 0.3/1000                    # Beam radius (m)
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

# Degrees of freedom (DOFs) definition
ndofs = nnodes * 6                     # Total number of DOFs (6 per node for displacement and rotation)
blocked_dofs = Int[2:6; 6*(nnodes/2+1-1).+(1:6)]   # Fixed DOFs (e.g., at the boundary)
encastre = Encastre(blocked_dofs)       # Encastre structured

# Displacement imposed
imposed_dofs = Int[1]
dispA = 10/100           # Peak displacement 
tA   = 2.5             # Duration of each linear segment
k    = dispA / tA      # Constant slope for rising/falling legs
displacement_applied(t) = (k*t).*(t<=tA) + (-k*t + 2*dispA).*(t>tA && t<=2*tA) + (k*(t-2*tA)).*(t>2*tA && t<=3*tA) +(-k*t + 4*dispA).*(t>3*tA && t<=4*tA)

displacementimposed = ImposedDisplacement(imposed_dofs,displacement_applied)

# External force
loaded_dofs = Int[]   # Degree of freedom where the external force is applied
force_function(t,i) = 0 # Time-dependent external force function
concentrated_force = ConcentratedForce(force_function, loaded_dofs)  

# Beam configuration struct initialization
conf = BeamsConfiguration(nodes, beams, Loads(concentrated_force), BoundaryConditions(encastre,displacementimposed, ndofs))  # Store the configuration for beams

# Declare beam-to-beam contact
beam2beam = true
#----------------------------------------------------
# SOLVER DEFINITIONS
#----------------------------------------------------

# HHT (Houbolt-Hughes-Taylor) time stepping parameters
α = -0.05   # Typically between 0 and -1, used for numerical stability
β = 0.25 * (1 - α)^2  # Damping parameter based on α
γ = 0.5 * (1 - 2 * α)  # Time-stepping parameter

# General time stepping parameters
initial_timestep = 1E-4  # Initial time step size
min_timestep = 1E-10    # Minimum allowed time step
max_timestep = 1E-4    # Maximum allowed time step (could be adjusted based on system behavior)
output_timestep = 1E-2    # Time step for output plotting or visualization
simulation_end_time = 1 # End time for the simulation (duration of the analysis)

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
run_simulation!(conf, params,nothing, beam2beam)

