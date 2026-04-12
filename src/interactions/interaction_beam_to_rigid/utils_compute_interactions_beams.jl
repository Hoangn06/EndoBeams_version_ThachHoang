#----------------------------------
# UTILS FUNCTIONS FOR COMPUTING THE CONTACT WITH PLANE SURFACE
#----------------------------------

# Compute contact gap, gradient, and Hessian for a plane surface.
@inline function contact_gap_beams(point, s::PlaneSurface, radius_beam)
    gₙ     = dot(s.normal, point) - s.offset - radius_beam
    ∂gₙ∂x  = s.normal
    ∂²gₙ∂x² = @SMatrix zeros(3,3)
    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
end

# Determine if a point is in contact with a plane surface.
@inline function incontact_beams(point, s::PlaneSurface, radius_beam::Float64)
    ḡₙ = radius_beam / 4
    gₙ = dot(s.normal, point) - s.offset - radius_beam
    return gₙ ≤ ḡₙ
end

#----------------------------------
# UTILS FUNCTIONS FOR COMPUTING THE CONTACT WITH SPHERE SURFACE
#----------------------------------

# Compute contact gap, gradient, and Hessian for a spherical
@inline function contact_gap_beams(point, master_surface::SphereSurface, radius_beam)
    # Vector from the point to the center of the sphere
    aux = point - master_surface.center
    
    # Compute the norm (distance) from the point to the surface center
    norm_aux = norm(aux)
    
    # Gap is the distance to the surface minus the radius of the beam
    gₙ = norm_aux - master_surface.radius - radius_beam
    
    # Gradient of the gap with respect to the point's coordinates
    invnorm = 1 / norm_aux  # Inverse of the distance
    ∂gₙ∂x = invnorm * aux  # Gradient (direction of normal)
    
    # Hessian of the gap (second derivative)
    ∂²gₙ∂x² = invnorm * ID3 - (invnorm^3) * (aux * aux')
    
    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
end

# Determine if a point is in contact with a spherical SDF based on a threshold distance.
@inline function incontact_beams(point, master_surface::SphereSurface, radius_beam::Float64)

    # Define the threshold distance
    ḡₙ = radius_beam/4

    # Vector from the point to the center of the sphere
    aux = point - master_surface.center
    
    # Compute the distance to the surface and compare to the threshold
    gₙ = norm(aux) - master_surface.radius - radius_beam
    
    # If gₙ is less than or equal to the threshold gap, then the point is in contact
    return gₙ ≤ ḡₙ
end

#----------------------------------
# UTILS FUNCTIONS FOR COMPUTING THE CONTACT WITH CYLINDER SURFACE
#----------------------------------

# Compute contact gap, gradient, and Hessian for a cylindrical surface
@inline function contact_gap_beams(point, s::CylinderSurface, radius_beam)
    v      = point - s.center              # vector from axis point to beam point
    v_ax   = dot(v, s.axis) * s.axis       # axial component (along cylinder axis)
    v_perp = v - v_ax                      # radial component (perpendicular to axis)
    d      = norm(v_perp)                  # radial distance from axis
    gₙ     = d - s.radius - radius_beam
    invd   = 1 / d
    n̂      = invd * v_perp                 # unit radial normal
    ∂gₙ∂x  = n̂
    # Hessian: curvature in circumferential direction only
    ∂²gₙ∂x² = invd * (ID3 - s.axis * s.axis' - n̂ * n̂')
    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
end

# Determine if a point is in contact with a cylindrical surface
@inline function incontact_beams(point, s::CylinderSurface, radius_beam::Float64)
    ḡₙ     = radius_beam / 4
    v      = point - s.center
    v_perp = v - dot(v, s.axis) * s.axis
    gₙ     = norm(v_perp) - s.radius - radius_beam
    return gₙ ≤ ḡₙ
end

#----------------------------------
# UTILS FUNCTIONS FOR COMPUTING THE CONTACT WITH MOVING CYLINDER SURFACE
#----------------------------------

# Compute contact gap, gradient, and Hessian for a moving cylindrical surface
@inline function contact_gap_beams(point, s::MovingCylinderSurface, radius_beam)
    v      = point - s.center
    v_ax   = dot(v, s.axis) * s.axis
    v_perp = v - v_ax
    d      = norm(v_perp)
    gₙ     = d - s.radius - radius_beam
    invd   = 1 / d
    n̂      = invd * v_perp
    ∂gₙ∂x  = n̂
    ∂²gₙ∂x² = invd * (ID3 - s.axis * s.axis' - n̂ * n̂')
    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
end

# Determine if a point is in contact with a moving cylindrical surface
@inline function incontact_beams(point, s::MovingCylinderSurface, radius_beam::Float64)
    ḡₙ     = radius_beam / 4
    v      = point - s.center
    v_perp = v - dot(v, s.axis) * s.axis
    gₙ     = norm(v_perp) - s.radius - radius_beam
    return gₙ ≤ ḡₙ
end

#----------------------------------
# UTILS FUNCTIONS FOR COMPUTING THE CONTACT WITH MOVING SPHERE SURFACE
#----------------------------------

# Contact gap, gradient, and Hessian for a moving sphere 
@inline function contact_gap_beams(point, master_surface::MovingSphereSurface, radius_beam)
    aux = point - master_surface.center
    norm_aux = norm(aux)
    gₙ = norm_aux - master_surface.radius - radius_beam
    invnorm = 1 / norm_aux
    ∂gₙ∂x = invnorm * aux
    ∂²gₙ∂x² = invnorm * ID3 - (invnorm^3) * (aux * aux')
    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
end

# Determine if a point is in contact with a moving sphere surface
@inline function incontact_beams(point, master_surface::MovingSphereSurface, radius_beam::Float64)
    ḡₙ = radius_beam / 4
    aux = point - master_surface.center
    gₙ = norm(aux) - master_surface.radius - radius_beam
    return gₙ ≤ ḡₙ
end

#----------------------------------
# UTILS FUNCTIONS FOR COMPUTING THE CONTACT WITH MOVING PLANE SURFACE
#----------------------------------

# Contact gap, gradient, and Hessian for a moving plane
@inline function contact_gap_beams(point, s::MovingPlaneSurface, radius_beam)
    gₙ       = dot(s.normal, point) - s.offset - radius_beam
    ∂gₙ∂x    = s.normal
    ∂²gₙ∂x²  = @SMatrix zeros(3,3)
    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
end

# Determine if a point is in contact with a moving plane surface
@inline function incontact_beams(point, s::MovingPlaneSurface, radius_beam::Float64)
    ḡₙ = radius_beam / 4
    gₙ = dot(s.normal, point) - s.offset - radius_beam
    return gₙ ≤ ḡₙ
end


#----------------------------------
# UTILS FUNCTIONS FOR COMPUTING THE CONTACT WITH DISCRETE SDF SURFACE
#----------------------------------

# Check if a point is within the specified domain.
@inline function isinside(point, dom)
    l_x = point[1] - dom[1]  
    flag_lx = l_x >= 0 && l_x <= (dom[2] - dom[1]) 
    l_y = point[2] - dom[3]
    flag_ly = l_y >= 0 && l_y <= (dom[4] - dom[3])
    l_z = point[3] - dom[5] 
    flag_lz = l_z >= 0 && l_z <= (dom[6] - dom[5])
    return flag_lx && flag_ly && flag_lz
end

# Compute contact gap, gradient, and Hessian for a discrete SDF.
@inline function contact_gap_beams(point, sdf::DiscreteSignedDistanceField, radius_beam::Float64)
    
    sitp = sdf.sitp
    gₙ = sitp(point...)
    ∂gₙ∂x = Interpolations.gradient(sitp, point...)
    ∂²gₙ∂x² = Interpolations.hessian(sitp, point...)
    
    # Normalize gradient for continuous contact calculations.
    nn = dot(∂gₙ∂x, ∂gₙ∂x)
    nmaginv = 1 / sqrt(nn)
    ∂gₙ∂x = ∂gₙ∂x * nmaginv 
    ∂²gₙ∂x² = ∂²gₙ∂x² * nmaginv * (ID3 - (∂gₙ∂x * ∂gₙ∂x') / nn)
    
    return gₙ - radius_beam, ∂gₙ∂x, ∂²gₙ∂x²
end 

# Check if a point is in contact with a discrete SDF.
@inline function incontact_beams(point, sdf::DiscreteSignedDistanceField, radius_beam::Float64)

    # Define the threshold distance
    ḡₙ = radius_beam/4

    # Check if the point is inside the SDF domain and compute the gap
    return isinside(point, sdf.dom) && sdf.sitp(point...) - radius_beam ≤ ḡₙ
    
end 
