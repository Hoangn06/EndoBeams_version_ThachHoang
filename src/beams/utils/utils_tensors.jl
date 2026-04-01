#---------------------------------------------------
# FUNCTIONS RELATED TO THE VOIGT TENSOR
#---------------------------------------------------
# Function to compute the volumetric part of a tensor (tr(t)=t11 + t22 + t33)
function volumetric_tensor(tensor::AbstractVector)
    return tensor[1]
end

# Function to compute the deviatoric part of a tensor (t - tr(t)/3 * I)
function deviatoric_tensor(tensor::AbstractVector)
    tr_over_3 = volumetric_tensor(tensor) / 3
    return @SVector [
        tensor[1] - tr_over_3,
        tensor[2],
        tensor[3]
    ]
end


#Function to compute the norm of a tensor in voigt notation which related to the strain
# The tensor's components in 12, 13 and 23 are scaled by 2
function norm_strain(tensor::AbstractVector)
    return sqrt(
        tensor[1]^2 + 
        2.0 * (tensor[2] / 2.0)^2 +
        2.0 * (tensor[3] / 2.0)^2
    )
end

# Function to compute the norm of a tensor in voigt notation which related to the stress
function norm_stress(tensor::AbstractVector)
    return sqrt(
        tensor[1]^2 + 
        2.0 * tensor[2]^2 +
        2.0 * tensor[3]^2
    )
end

# Function to compute the norm of the deviatoric stress tensor in voigt notation
function norm_stress_deviatoric(tensor::AbstractVector)
    return sqrt(
        (3/2)*tensor[1]^2 + 
        2.0 * tensor[2]^2 +
        2.0 * tensor[3]^2
    )
end

# Function to compute the norm of the deviatoric strain tensor in voigt notation
function norm_strain_deviatoric(tensor::AbstractVector)
    return sqrt(
        (3/2)*tensor[1]^2 + 
        2.0 * (tensor[2] / 2.0)^2 +
        2.0 * (tensor[3] / 2.0)^2
    )
end


# Compute the matrix constitutive for the constitutive law of the material
function matrix_constitutive(E::Float64, G::Float64)
    return @SMatrix [
        E  0.0  0.0;
        0.0  G  0.0;
        0.0  0.0  G
    ]
end

# Compute the matrix constitutive for the deviatoric part of the stress
function matrix_constitutive_deviatoric(E::Float64, G::Float64)
    return @SMatrix [
        (2.0/3.0)*E  0.0  0.0;
        0.0          G  0.0;
        0.0          0.0  G
    ]
end


#---------------------------------------------------
# Construction of the Matrix of shape function derivatives
#---------------------------------------------------
function matrix_Be(x1::Float64, x2::Float64, x3::Float64, l0::Float64)
    Be = zeros(3, 7)

    Be[1,1] = 1.0/l0
    Be[1,3] = x3 * ((6.0*x1)/(l0^2) - (4.0/l0) )
    Be[1,4] = x2 * (-(6.0*x1)/(l0^2) + (4.0/l0) )
    Be[1,6] = x3 * ((6.0*x1)/(l0^2) - (2.0/l0) )
    Be[1,7] = x2 * (-(6.0*x1)/(l0^2)+ (2.0/l0) )

    Be[2,2] = (x3/l0)
    Be[2,5] = -(x3/l0)   
    
    Be[3,2] = -(x2/l0)
    Be[3,5] = (x2/l0)

    
    return Mat37{Float64}(Be)
end
