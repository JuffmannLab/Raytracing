
abstract type AbstractElectron end

mutable struct Electron1d <: AbstractElectron
    ψ::Vector{Array}
    intensity::Vector{<:Real}
    p_ges::Vector{<:Real}
    x::Vector{<:Real}
end

"""
    Electron(x::Vector{<:Real}, intensity::Vector{<:Real}, energy::Real, n::Int64)::Electron1d

Return a Electron1d type.

Create and return the 1D wave function. The `x` value corresponds to the
x-axis, and `intensity` corresponds to the intensity distribution along
the given axis. `n` correspond to the number of rays that should be simulated.
"""
function Electron(x::Vector{<:Real}, intensity::Vector{<:Real}, energy::Real, n::Int64)::Electron1d
   

    # TODO:
    # implement the random ray positions
    # create the new position angle vectors
    ψ = Vector{Array}(undef, size(x, 1))
    p_ges = similar(intensity)
    for i = 1:size(ψ, 1)
        ψ[i] = [x[i], 0.]
        momentum = sqrt( 2 * energy * q / m_e ) * m_e
        p_ges[i] = momentum
    end

    return Electron1d(ψ, intensity, p_ges, x)
end

mutable struct Electron2d <: AbstractElectron
    ψ::Vector{Array}
    intensity::Vector{<:Real}
    p_ges::Vector{<:Real}
    x::Vector{<:Real}
    y::Vector{<:Real}
end

"""
    Electron(x::Vector{<:Real}, y::Vector{<:Real}, intensity::Matrix{<:Real},
             energy::Real, n::Int64)::Electron2d

Return the Electron2d type.

Construct and return a 2D wave. The `x` and `y` vector correspond to the
x and y axis of the problem, and the `intensity` is the intensity corresponding
to the x and y axis given prior. The `energy` parameter denotes the electron
energy in eV. The number of rays that should be simulated are cenoted by the `n`
value.
"""
function Electron(x::Vector{<:Real}, y::Vector{<:Real}, intensity::Matrix{<:Real},
                  energy::Real, n::Int64)::Electron2d
    
    # create the array with the different angles and positions
    ψ = Vector{Array}(undef, n)

    # create a flattened intensity array
    int_flat = similar(intensity, n)
    p_ges = similar(int_flat)

    # calculate the momentum
    momentum = sqrt( 2 * energy * q / m_e ) * m_e

    # calculate the span of x and y, aswell as Δx and Δy
    span_x = abs(x[1]-x[end])
    span_y = abs(y[1]-y[end])
    Δx = abs(x[2]-x[1])
    Δy = abs(y[2]-y[1])

    # iterate over the ψ array, fill it with random coordinates
    for i = 1:n

        # generate random coordinates that are evenly distributed in x and y direction
        coords_rand = rand(Float64, 2)-0.5
        x_temp = coords_rand[1] * span_x
        y_temp = coords_rand[2] * span_y

        # get the position of the random coordinates in the intensity array
        m_x = round(Int, (x_temp-x[1]) / Δx, RoundDown)+1
        m_y = round(Int, (y_temp-y[1]) / Δy, RoundDown)+1

        # fill the ψ, p_ges and flat_int vectors
        ψ[n] = [x_temp, y_temp, 0. 0.]
        int_flat[n] = intensity[m_x, m_y]
        p_ges[n] = momentum

    end

    # return the wanted struct
    return Electron2d(ψ, int_flat, p_ges, x, y)
end

"""
    getintensity(wave::Electron)::Vector{<:Real}

Return the intensity.

Calculate the current intensity and return it. The `wave` value is the wave struct
from which the intensity should be returned.
"""
function getintensity(wave::Electron1d)::Vector{<:Real}
    
    # intensity is returned in the same coordinate system as 
    # the wave object has at the starting point!
    x = wave.x
    intensity = zeros(Float64, size(wave.intensity))
    dx = abs(x[1]-x[2])

    # loop over the rays
    for i = 1:size(wave.ψ, 1)

        # loop over the coordinate system
        for j = 1:size(x, 1)

            # check if the ray is near the current coorinate
            if abs(wave.ψ[i][1]-x[j]) <= dx/2

                # add the intensity from the incident plane
                # to the output plane
                intensity[j] += wave.intensity[i]

                # if a incident point is assigned to a output intensity
                # break the loop (for the sake of rounding and energy conservation)
                break
            else end
        end
    end

    # return the intensity at the different points
    return intensity
end

"""
    getintensity(ray::Electron2d), x::Vector{<:Real}, y::Vector{<:Real}::Matrix{<:Real}

Return the intensity.

Calculate the corrent intensity and return it. The `ray` value is the electron beam struct
from which the intensity should be returned, the `x` and `y` values are the coordinates
in which the intensity should be returned.
"""
function getintensity(ray::Electron2d, x::Vector{<:Real}, y::Vector{<:Real})::Matrix{<:Real}
    
    # create the empty intensity array
    intensity = zeros(Float64, size(x, 1), size(y, 1))
    Δx = abs(x[1]-x[2])
    Δy = abs(y[1]-y[2])
 
    # loop over the rays
    for i = 1:size(ray.ψ, 1)

        # calculate the position of the intensity in the new grid
        # the minus the first element "normalizes" the grid
        n = round(Int, (ray.ψ[i][1]-x[1]) / Δx, RoundDown) + 1
        m = round(Int, (ray.ψ[i][2]-y[1]) / Δy, RoundDown) + 1

        # add the intensity to the grid point
        try
            intensity[n, m] += ray.intensity[i]
        catch e
            # catch a BoundsError, and do nothing, because there is nothing to do
            if isa(e, BoundsError)
            end
        end
    end

    # return the intensity
    return intensity
end
