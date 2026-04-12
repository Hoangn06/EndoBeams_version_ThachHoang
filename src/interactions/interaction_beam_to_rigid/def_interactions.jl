#----------------------------------
# INTERACTION PROPERTIES DEFINITION
#----------------------------------

# Defines the properties for contact interactions (e.g., friction, penalty method).
struct InteractionProperties
    kₙ::Float64   # Normal penalty parameter (stiffness)
    μ::Float64   # Friction coefficient
    εᵗ::Float64   # Regularized parameter for friction contact
    ηₙ::Float64   # Damping parameter in the normal direction
    kₜ::Float64   # Tangential penalty parameter
    ηₜ::Float64   # Damping parameter in the tangential direction
    u̇ₛ::Float64   # Slip velocity threshold for friction regularization
end

abstract type RigidBodySurface end
abstract type DeformableSurface end

#----------------------------------
# BEAMS SURFACE DEFINITIONS
#----------------------------------
# Defines a beam surface composed of beam segments (for contact with rigid bodies)
struct BeamElementSurface <: DeformableSurface
    contact_beams::Vector{Int}  
end

# Creates a BeamElementSurface from a connectivity matrix.
function BeamElementSurface(beams_connectivity::Array{Int, 2})
    contact_beams = collect(1:size(beams_connectivity, 1)) # Creates a StructArray of beam indices
    return BeamElementSurface(contact_beams)
end

#----------------------------------
# SURFACE FIXED DEFINITIONS
#----------------------------------

# Defines an analytical plane used as a master surface in contact problems.
struct PlaneSurface <: RigidBodySurface
    normal::Vec3{Float64}  # Outward unit normal of the plane
    offset::Float64        # Signed distance from origin: dot(normal, x) = offset on the plane
end


# Defines an analytical sphere used as a master surface in contact problems
struct SphereSurface <: RigidBodySurface
    center::Vec3{Float64}  # Sphere center (x, y, z)
    radius::Float64  # Sphere radius
end

# Defines an analytical cylinder used as a master surface in contact problems
struct CylinderSurface <: RigidBodySurface
    center::Vec3{Float64}  # Cylinder center (x, y, z)
    axis::Vec3{Float64}  # Cylinder axis
    radius::Float64  # Cylinder radius
end

#----------------------------------
# SURFACE MOVING DEFINITIONS
#----------------------------------
# Defines a sphere whose center moves over time according to a user-supplied function
mutable struct MovingSphereSurface{TF} <: RigidBodySurface
    center::MVector{3, Float64}  # Current sphere center (updated each time step)
    radius::Float64              # Sphere radius
    center_function::TF          # Function center_function(t) -> [x, y, z]
end

# Defines a moving sphere surface from a center function and a radius
function MovingSphereSurface(center_function::TF, radius::Float64) where TF
    c0 = center_function(0.0)
    return MovingSphereSurface{TF}(MVector{3, Float64}(c0), radius, center_function)
end

# Update the sphere center for the current time step
function update_surface!(s::MovingSphereSurface, t)
    c = s.center_function(t)
    s.center[1] = c[1]
    s.center[2] = c[2]
    s.center[3] = c[3]
end


# ---------------------------------
# Defines a cylinder whose center moves over time according to a user-supplied function
mutable struct MovingCylinderSurface{TF} <: RigidBodySurface
    center::MVector{3, Float64}  # Current cylinder center (updated each time step)
    axis::Vec3{Float64}  # Cylinder axis
    radius::Float64              # Cylinder radius
    center_function::TF          # Function center_function(t) -> [x, y, z]
end

# Defines a moving cylinder surface from a center function and a radius
function MovingCylinderSurface(center_function::TF, axis::Vec3{Float64}, radius::Float64) where TF
    c0 = center_function(0.0)
    return MovingCylinderSurface{TF}(MVector{3, Float64}(c0), axis, radius, center_function)
end

# Update the cylinder center for the current time step
function update_surface!(s::MovingCylinderSurface, t)
    c = s.center_function(t)
    s.center[1] = c[1]
    s.center[2] = c[2]
    s.center[3] = c[3]
end
#----------------------------------

# Defines a plane whose offset moves over time according to a user-supplied function.
mutable struct MovingPlaneSurface{TF} <: RigidBodySurface
    normal::Vec3{Float64}   # Fixed outward unit normal of the plane
    offset::Float64          # Current offset (updated each time step)
    offset_function::TF      # Function offset_function(t) -> Float64
end

# Defines a moving plane surface from a fixed normal and an offset function
function MovingPlaneSurface(normal::Vec3{Float64}, offset_function::TF) where TF
    return MovingPlaneSurface{TF}(normal, offset_function(0.0), offset_function)
end

# Update the plane offset for the current time step
function update_surface!(s::MovingPlaneSurface, t)
    s.offset = s.offset_function(t)
end

update_surface!(::RigidBodySurface, t) = nothing # No action for fixed surfaces

#----------------------------------
# SDF SURFACE DEFINITIONS
#----------------------------------

# Defines a discrete signed distance field used as a master surface in contact problems
struct DiscreteSignedDistanceField{Tsitp} <: RigidBodySurface
    flag_load_from_file_iterative::Bool
    sitp::Tsitp                  # Scaled interpolation of the SDF field
    dom::NTuple{6, Float64}  # Domain boundary coordinates
    dx::Float64              # Grid spacing in x direction
    dy::Float64              # Grid spacing in y direction
    dz::Float64              # Grid spacing in z direction
end

# Defines a triangulated surface used as a master surface in contact problems
struct TriangulatedSurface <: RigidBodySurface
    positions::Vector{Vec3{Float64}}  
    triangles::Vector{Vec3{Int}} 
end


#----------------------------------
# SDF SURFACE CREATION FUNCTIONS
#----------------------------------

# Creates a DiscreteSignedDistanceField reading the SDF data from a VTK file and setting up interpolation.
function DiscreteSignedDistanceField(filename, inside, flag_load_from_file_iterative=false)

    # Read SDF data from VTK file
    npx, npy, npz, dx, dy, dz, dom, sdf = read_vtk_sdf(filename)
    field = inside ? reshape(sdf, (npx, npy, npz)) : reshape(-sdf, (npx, npy, npz))

    # Define coordinate ranges for interpolation
    x, y, z = range(dom[1]; step=dx, stop=dom[2]), range(dom[3]; step=dy, stop=dom[4]), range(dom[5]; step=dz, stop=dom[6])

    # Create a quadratic interpolation for smooth gradients
    itp = interpolate(field, BSpline(Quadratic(Reflect(OnCell()))))
    sitp = scale(itp, x, y, z)  # Scaled interpolation over the coordinate grid

    return DiscreteSignedDistanceField{typeof(sitp)}(flag_load_from_file_iterative, sitp, dom, dx, dy, dz)  
    
end 

#------------------------------------------------
# BEAM-TRIANGLE CONNECTION DEFINITION 
#------------------------------------------------

# Represents a beam-triangle interaction connection.
struct BeamTriangleConnection
    beam_node::Int  # Beam node index
    triangle_nodes::Tuple{Int, Int, Int}  # Triangle node indices
    global_sparsity_map::SVector{144, Int}
end

#----------------------------------
# INTERACTION STRUCTURE DEFINITIONS
#----------------------------------

abstract type Interaction end

# Defines a rigid interaction between a master surface and a slave surface.
struct RigidInteraction <: Interaction
    master::RigidBodySurface
    slave::BeamElementSurface
    properties::InteractionProperties
end

# Defines a collection of rigid interactions (beams vs multiple surfaces).
struct MultiRigidInteraction <: Interaction
    interactions::Vector{RigidInteraction}
end
