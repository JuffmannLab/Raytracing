
using FFTW

abstract type AbstractElectron end


mutable struct Electron <: AbstractElectron
    ψ::Vector{Array}
    intensity::Vector{<:Real}
    v::Real
    x::Vector{<:Real}
    y::Vector{<:Real}
end

"""
    Electron(x::Vector{<:Real}, y::Vector{<:Real}, intensity::Matrix{<:Real},
             U::Real, n::Int64)::Electron2d

Return the Electron2d type.

Construct and return a 2D wave. The `x` and `y` vector correspond to the
x and y axis of the problem, and the `intensity` is the intensity corresponding
to the x and y axis given prior. The `U` parameter denotes the electron
energy in eV. The number of rays that should be simulated are denoted by the `n`
value.
"""
function Electron(x::Vector{<:Real}, y::Vector{<:Real}, intensity::Matrix{<:Real},
                  U::Real, n::Int)::Electron
    
    # create the array with the different angles and positions
    ψ = Vector{Array}(undef, n)

    # create a flattened intensity array
    int_flat = similar(intensity, n)

    # calculate the relativistic electron velocity
    v = c * sqrt(1 - 1 / (1 + U * q / m_e / c^2)^2)

    # calculate the span of x and y, aswell as Δx and Δy
    span_x = abs(x[1]-x[end])
    span_y = abs(y[1]-y[end])
    Δx = abs(x[2]-x[1])
    Δy = abs(y[2]-y[1])

    # iterate over the ψ array, fill it with random coordinates
    for i = 1:n

        # generate random coordinates that are evenly distributed in x and y direction
        coords_rand = rand(Float64, 2) .- 0.5
        x_temp = coords_rand[1] * span_x
        y_temp = coords_rand[2] * span_y

        # get the position of the random coordinates in the intensity array
        m_x = round(Int, (x_temp-x[1]) / Δx, RoundDown) + 1
        m_y = round(Int, (y_temp-y[1]) / Δy, RoundDown) + 1

        # fill the ψ, p_ges and flat_int vectors
        ψ[i] = [x_temp, y_temp, 0., 0.]
        int_flat[i] = intensity[m_x, m_y]

    end

    # the wanted struct
    return Electron(ψ, int_flat, v, x, y)
end

"""
    getintensity(ray::Electron), x::Vector{<:Real}, y::Vector{<:Real}::Matrix{<:Real}

Return the intensity.

Calculate the current intensity and return it. The `ray` value is the electron beam struct
from which the intensity should be returned, the `x` and `y` values are the coordinates
in which the intensity should be returned.
"""
function getintensity(ray::Electron, x::Vector{<:Real}, y::Vector{<:Real})::Matrix{<:Real}
    
    # create the empty intensity array
    intensity = zeros(Float64, size(x, 1), size(y, 1))
    Δx = abs(x[1]-x[2])
    Δy = abs(y[1]-y[2])

    # loop over the rays
    for i in eachindex(ray.ψ)

        # calculate the position of the intensity in the new grid
        # the minus the first element "normalizes" the grid
        n = round(Int, (ray.ψ[i][1]-x[1]) / Δx, RoundDown) + 1
        m = round(Int, (ray.ψ[i][2]-y[1]) / Δy, RoundDown) + 1

        # add the intensity to the grid point
        if 0 < n < size(intensity, 1) + 1 && 0 < m < size(intensity, 2) + 1
            intensity[n, m] += ray.intensity[i]
        else
        end
    end

    # return the intensity
    return intensity
end

"""
mcp(intensity::Matrix{<:Real}, psf::Matrix{<:Real})::Matrix{<:Real}

Return the mcp electron distribution.

The light that leaves the mcp is calculated here. It is done by taking the input electron
intensity distribution in the `intensity` matrix and the point spread funciton `psf` of the
measuring device.
"""
function mcp(intensity::Matrix{<:Real}, psf::Matrix{<:Real})::Matrix{<:Real}
    INTENSITY = fft(intensity)
    PSF = fft(psf)
    PSF .*= INTENSITY
    return real.(ifftshift(ifft(PSF)))
end
