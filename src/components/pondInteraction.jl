
struct PondInteraction <: Component
    lf::AbstractLight
end


"""
    calculate!(ray::Electron, pond::PondInteraction)

Calculate the scattering angle.

Calculate the angle that the electron has after the interaction with the light field
and save it in the ray object. The input parameters are `pond` which is a light field
object, aswell as `ray` which is the electron beams object.
"""
function calculate!(ray::Electron, pond::PondInteraction)
 
    # get the light field out of the PondInteraction struct
    lf = pond.lf

    # calculate the intensity gradients
    (int_grad_x, int_grad_y) = _gradient(lf.intensity, lf.x, lf.y)

    # calculate the ω parameter
    ω = 2*π*c / lf.λ

    # calculate the Δt, Δx and Δy
    Δt = abs(lf.t[1]-lf.t[2])
    Δx = abs(lf.x[1]-lf.y[2])
    Δy = abs(lf.y[1]-lf.y[2])

    # calculate the time integral
    t_int = sum(lf.envelope) * Δt

    # normalization factor for the intensity
    I0 = lf.E / ( sum(lf.intensity) * Δx * Δy * t_int )

    # define α * ħ
    αħ = 1 / (4 * π * ε_0) * q^2 / c

    # define the beta (should it be v_z? i assume that it is not relevant, bc v ≈ v_z)
    β = ray.v / c

    # define the electron Energy and momentum
    γ = 1 / sqrt(1 - (ray.v / c)^2)
    Ee = γ * m_e * c^2
    p = γ * m_e * ray.v

    # calculate the constant that is needed for units sake
    # the plus before the beta denotes that the electron and the 
    # photon are counter propagating
    constant = αħ * lf.λ^2 / (2 * π * (1 + β)) * I0 / Ee

    # iterate over the all the electron beams
    for i = 1:size(ray.ψ, 1)

        # interpolate the x and y force-direction
        Fx = constant * _interpolation(int_grad_x, lf.x, lf.y, ray.ψ[i][1], ray.ψ[i][2])
        Fy = constant * _interpolation(int_grad_y, lf.x, lf.y, ray.ψ[i][1], ray.ψ[i][2])

        # calculate Δpx and Δpy
        Δpx = Fx * t_int
        Δpy = Fy * t_int

        # calculate the tangens and a constant
        tan_α = tan(ray.ψ[i][3])
        tan_β = tan(ray.ψ[i][4])
        A = sqrt(1 + tan_α^2 + tan_β^2)

        # calculate the momentum components
        pz = p / A
        px = pz * tan_α
        py = pz * tan_β
        
        # calculate the new angle
        # the complex - real transformation is to prevent floating point errors
        # due to floating point differences not being 0 when they should be
        ray.ψ[i][3] = atan((px + Δpx) / 
                           real(sqrt(complex(pz^2 * A^2 - (px + Δpx)^2 - (py + Δpy)^2))))
        ray.ψ[i][4] = atan((py + Δpy) / 
                           real(sqrt(complex(pz^2 * A^2 - (px + Δpx)^2 - (py + Δpy)^2))))

    end
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
    n = round(Int, (x-coords_x[1]) / Δx, RoundDown) + 1
    m = round(Int, (y-coords_x[1]) / Δy, RoundDown) + 1

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
