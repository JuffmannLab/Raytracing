
"""
    pondinteraction!(lf::LightField1d, ray::Electron1d)

Return the scattering angle.

Calculate the new propagation direction of all the electron rays, when it interacts
with a given light field that is defined in `LightField1d`. The electron are contained
in the `ray` parameter.
"""
function pondinteraction!(lf::LightField1d, ray::Electron1d)
    
    # calculate the intensity gradient
    int_grad = _gradient(lf.intensity, lf.x)

    # calculate the ω parameter
    ω = 2*π*c / lf.λ

    # calculate the Δt
    Δt = abs(lf.t[1]-lf.t[2])

    # calculate the time integral
    t_int = sum(lf.envelope)
    
    # calculate the constant that is needed for units sake
    constant = - q^2 / (4 * m_e * ω^2)

    # iterate over all the electron beams
    for i = 1:size(ray.ψ, 1)
        
        # calculate the Force at the position of the electron
        Fx = constant * _interpolation(int_grad, lf.x, ray[i].ψ[1])

        # calculate Δpx
        Δpx = Fx * Δt * t_int

        # calculate the tangens alpha value
        tan_α = tan(ray[i].ψ[2])
        γ = sqrt(1 + tan_α^2)

        # update the angle α
        ray[i].ψ[2] = atan(tan_α + Δpx * γ / ray[i].p_ges)

        # update the new momentum for this ray
        ray[i].p_ges = sqrt( (tan_α * ray[i].p_ges / γ)^2 + (ray[i].p_ges / γ)^2 )
    end
end

"""
    pondinteraction!(lf::LightField2d, ray::Electron2d)

Calculate the scattering angle.

Calculate the angle that the electron has after the interaction with the light field
and save it in the ray object. The input parameters are `lf` which is a light field
object, aswell as `ray` which is the electron beams object.
"""
function pondinteraction!(lf::LightField2d, ray::Electron2d)
    
    # calculate the intensity gradients
    (int_grad_x, int_grad_y) = _gradient(lf.intensity, lf.x, lf.y)

    # calculate the ω parameter
    ω = 2*π*c / lf.λ

    # calculate the Δt
    Δt = abs(lf.t[1]-lf.t[2])

    # calculate the time integral
    t_int = sum(lf.envelope)

    # calculate the constant that is needed for units sake
    constant = - q^2 / (4 * m_e * ω^2)

    # iterate over the all the electron beams
    for i = 1:size(ray.ψ, 1)
        # interpolate the x and y force-direction
        Fx = constant * _interpolation(int_grad_x, lf.x, lf.y, ray[i].ψ[1], ray[i].ψ[2])
        Fy = constant * _interpolation(int_grad_y, lf.x, lf.y, ray[i].ψ[1], ray[i].ψ[2])

        # calculate Δpx and Δpy
        Δpx = Fx * Δt * t_int
        Δpy = Fy * Δt * t_int

        # calculate the tangens and a constant
        tan_α = tan(ray[i].ψ[3])
        tan_β = tan(ray[i].ψ[4])
        γ = sqrt(1 + tan_α^2 + tan_β^2)
        
        # calculate the new angle
        ray[i].ψ[3] = atan(tan_α + Δpx * γ / ray[i].p_ges)
        ray[i].ψ[4] = atan(tan_β + Δpy * γ / ray[i].p_ges)

        # calculate the new momentum
        ray[i].p_ges = sqrt( (tan_α * ray[i].p_ges / γ + Δpx)^2 + 
                             (tan_β * ray[i].p_ges / γ + Δpy)^2 +
                             (ray[i].p_ges / γ)^2 )
    end
end

"""
    _gradient(V::Vector{<:AbstractFloat}, x::Vector{<:Real})::Vector{<:AbstractFloat}

Calculate the gradient.

The gradient of a Vector `V` is calculated and returned. This method is
implemented with the dirichlet boundary conditions, that at the outer boundary
all the values are 0. The `x` value corresponds to the coordinate system we are in.
"""
function _gradient(V::Vector{<:AbstractFloat}, x::Vector{<:Real})::Vector{<:AbstractFloat}
    
    # define the output gradient field Vector
    F = similar(V)
    F[1] = 0
    F[end] = 0

    # calculate the Δx value
    Δx = abs(x[1]-x[2])

    # fill the output Vector with the calculation
    for i = 2:size(V, 1)-1
        F[i] = (V[i+1]-V[i-1]) / 2 / Δx
    end

    # return the gradientfield
    return F
end

"""
    _gradient(V::Matrix{<:AbstractFloat}, x::Vector{<:Real},
              y::Vector{<:Real})::Tuple{Matrix{<:AbstractFloat},Matrix{<:AbstractFloat}}

Calculate the gradient.

This function calculates and returns the gradient of a 2D input potential `V`, and will 
return the Gradient in the form of the tuple corresponding to `(Fx, Fy)`, where `Fx` and `Fy`
are Matrices with the same size and datatype as the `V` value. The `x` and `y` values
correspond to the coordinate system of the viewed Potentials/Forces.
"""
function _gradient(V::Matrix{<:AbstractFloat}, x::Vector{<:Real},
                   y::Vector{<:Real})::Tuple{Matrix{<:AbstractFloat},Matrix{<:AbstractFloat}}
    
    # define the output gradient vector field
    Fx = similar(V)
    Fy = similar(V)
    
    # set the boundaries
    Fx[1, :] .= 0
    Fy[1, :] .= 0
    Fx[end, :] .= 0
    Fy[end, :] .= 0
    Fx[:, 1] .= 0
    Fy[:, 1] .= 0
    Fx[:, end] .= 0
    Fy[:, end] .= 0

    # calculate the Δx and Δy value
    Δx = abs(x[1]-x[2])
    Δy = abs(y[1]-y[2])

    # calculate the gradient
    for j = 2:size(V, 2)-1
        for i = 2:size(V, 1)-1
            Fx[i, j] = (V[i+1,j] - V[i-1,j]) / 2 / Δx 
            Fy[i, j] = (V[i,j+1] - V[i,j-1]) / 2 / Δy
        end
    end

    return (Fx, Fy)
end

"""
    _interpolation(A::Vector{<:Real}, coords::Vector{<:Real}, x::Real)::Real

Return the linear interpolation.

Here the linear interpolation on a given Vector is calculated and returned.
The `A` value denotes the vector that will be interpolated, the `coords`
parameter is the coordinate system of the interpolation data `A`. `x` is
the point that should be interpolated.
"""
function _interpolation(A::Vector{<:Real}, coords::Vector{<:Real}, x::Real)::Real

    # if the x value is not in the area of coords return 0
    if x < coords[1] || x > coords[end]
        return zero(eltype(coords))
    end

    # calculate the Δx value
    Δx = abs(coords[1]-coords[2])

    # calculate the bin after which the x value lies
    n = round(Int, x / Δx, RoundDown) + 1

    # define the k value of the linear function
    k = (A[n] - A[n+1]) / (coords[n] - coords[n+1])

    # define the d value of the linear function
    d = A[n] - k * coords[n]

    # return the linear interpolation
    return k * x + d
end

"""
    _interpolation(A::Matrix{<:Real}, coords_x::Vector{<:Real},
                   coords_y::Vector{<:Real}, x::Real, y::Real)::Real

Return the bilinear interpolation.

Calculate and return the bilinear interpolation on the Matrix `A`
with the coordinates `coords_x` and `coords_y` at the point with
the position `x`, `y`.
"""
function _interpolation(A::Matrix{<:Real}, coords_x::Vector{<:Real},
                        coords_y::Vector{<:Real}, x::Real, y::Real)::Real

    # if the value is not in the area of coords return 0
    if x < coords_x[1] || x > coords_x[end]
        return zero(eltype(coords_x))
    elseif y < coords_y[1] || y > coords_y[end]
        return zero(eltype(coords_y))
    end

    # calculate the Δx and Δy values
    Δx = abs(coords_x[1]-coords_x[2])
    Δy = abs(coords_y[1]-coords_y[2])

    # calculate the bins in which the x and y values are
    n = round(Int, x / Δx, RoundDown) + 1
    m = round(Int, y / Δy, RoundDown) + 1

    # interpolation in x direction on both y-points
    k1 = (A[n, m] - A[n+1, m]) / (coords_x[n] - coords_y[n+1])
    d1 = A[n, m] - k1 * coords_x[n]

    k2 = (A[n, m+1] - A[n+1, m+1]) / (coords_x[n] - coords_y[n+1])
    d2 = A[n, m+1] - k2 * coords_x[n]

    A1 = k1 * x + d1
    A2 = k2 * x + d2

    # use the interpolated values on the x-axis for the interpolation in the y axis
    k = (A1 - A2) / (coords_y[m] - coords_y[m+1])
    d = A1 - k * coords_y[m]

    # return the interpolated value
    return k * y + d
end
