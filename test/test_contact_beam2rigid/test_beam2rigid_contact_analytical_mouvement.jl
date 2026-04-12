
#----------------------------------------------------
# PACKAGE INITIALIZATION
#----------------------------------------------------

using Pkg
using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, LinearSolve, Test, Revise, EndoBeams, BenchmarkTools, StaticArrays
using Plots

#----------------------------------------------------# BEAM DEFINITIONS
#----------------------------------------------------

# Define node positions and connectivity for beams
y_coords = collect(-5.0:1.0:5.0)          # 11 nodes from Y = -5 to Y = 5
positions = [zeros(11) y_coords zeros(11)]
nnodes = size(positions, 1)    # Total number of nodes

nbeams = nnodes - 1                    # Number of beam elements
connectivity = hcat(1:10, 2:11)        # Connectivity for 10 beam elements


# Initial conditions for displacements, velocities, accelerations, and rotations
initial_displacements = zeros(size(positions))       # Zero initial displacements
initial_velocities = zeros(size(positions))       # Zero initial velocities
initial_accelerations = zeros(size(positions))       # Zero initial accelerations
initial_rotations = zeros(size(positions))        # Zero initial rotations
initial_angular_velocities = zeros(size(positions))   # Zero initial angular velocities
initial_angular_accelerations = zeros(size(positions))  # Zero initial angular accelerations
plane = "xy"                         # Plane of the problem

# Material and geometric properties for beams
E = 50000                              # Young's modulus (MPa)
ν = 0.3                              # Poisson's ratio
ρ = 6.5e-9                                # Density (tonne/mm³)
radius = 0.25                         # Beam radius (mm)
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

# Define boundary conditions/ degrees of freedom
ndofs = nnodes * 6                     # Total number of DOFs
blocked_dofs = [2, 3, 62, 63]  # Node 1 and node 11: Y and Z displacements fixed
encastre = Encastre(blocked_dofs)

# Beam configuration struct initialization
conf = BeamsConfiguration(nodes, beams, nothing, BoundaryConditions(encastre, ndofs))


#----------------------------------------------------
# SOLVER DEFINITIONS
#----------------------------------------------------

# HHT (Houbolt-Hughes-Taylor) time stepping parameters
α = -0.2   # Typically between 0 and -1, used for numerical stability
β = 0.25 * (1 - α)^2  # Damping parameter based on α
γ = 0.5 * (1 - 2 * α)  # Time-stepping parameter

# General time stepping parameters
initial_timestep = 1e-3    # Initial time step size
min_timestep = 1e-10       # Minimum allowed time step
max_timestep = 1e-3        # Maximum allowed time step
output_timestep = 1e-2     # Time step for output plotting
simulation_end_time = 1.0     # End time for the simulation

# Convergence criteria for the solver
tolerance_residual = 1e-4       # Residual tolerance for convergence checks
tolerance_displacement = 1e-4   # Tolerance for changes in displacement (ΔD)
max_iterations = 10             # Maximum number of iterations for the solver

#----------------------------------
# INTERACTIONS
#----------------------------------

# Create interaction properties (using defaults for some values)
kₙ = 1e+3 # Penalty parameter
μ = 0.2 # Friction coefficient
εᵗ = 1e+1 # Regularized parameter for friction contact
ηₙ = 1e-3 
kₜ = kₙ
ηₜ = ηₙ
u̇ₛ = 0  
inter_properties = InteractionProperties(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ)

# Create master and slave surfaces

surface_slave = BeamElementSurface(connectivity)

# Moving plane: plane travels from x=1 to x=-1 over the simulation, pushing the beam in -X.
# Normal points toward the beam (-X direction). Plane at x=p → offset = -p.
plane_normal = SVector(-1.0, 0.0, 0.0)
offset_function(t) = -1.0 + 2.0 * t / simulation_end_time   # plane: x=1 at t=0, x=-1 at t=end
surface_master = MovingPlaneSurface(plane_normal, offset_function)

# Create the interaction instance
inter = MultiRigidInteraction([
    RigidInteraction(surface_master, surface_slave, inter_properties),
])

# Store solver parameters in a structured Params object
params = SimulationParams(;
α, β, γ, initial_timestep, min_timestep, max_timestep, output_timestep, simulation_end_time, tolerance_residual, tolerance_displacement, max_iterations, output_dir = "test/output3D"
, verbose = true)

#----------------------------------------------------
# START SIMULATION
#----------------------------------------------------

# Run the solver with the defined configurations and parameters
run_simulation!(conf, params, inter)


