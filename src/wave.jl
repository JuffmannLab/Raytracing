
abstract type AbstractWave end

mutable struct Wave1d
    ψ::Vector{Array}
    intensity::Vector{<:Real}
    x::Vector{<:Real}
end

"""
    Wave(x::Vector{<:Real}, intensity::Vector{<:Real})::Wave1d

Return a Wave1d type.

Create and return the 1D wave function. The `x` value corresponds to the
x-axis, and `intensity` corresponds to the intensity distribution along
the given axis.
"""
function Wave(x::Vector{<:Real}, intensity::Vector{<:Real})::Wave1d
   
    # create the new position angle vectors
    ψ = Vector{Array}(undef, size(x, 1))
    for i = 1:size(ψ, 1)
        ψ[i] = [x[i] 0]
    end

    return Wave(ψ, intensity, x)
end

mutable struct Wave2d
    ψ::Vector{Array}
    intensity::Vector{<:Real}
    x::Vector{<:Real}
    y::Vector{<:Real}
end

"""
    Wave(x::Vector{<:Real}, y::Vector{<:Real}, intensity::Matrix{<:Real})::Wave2d

Return the Wave2d type.

Construct and return a 2D wave. The `x` and `y` vector correspond to the
x and y axis of the problem, and the `intensity` is the intensity corresponding
to the x and y axis given prior.
"""
function Wave(x::Vector{<:Real}, y::Vector{<:Real}, intensity::Matrix{<:Real})::Wave2d
    
    # create the array with the different angles and positions
    ψ = Vector{Array}(undef, size(x, 1) * size(y, 1))

    # create a flattened intensity array
    int_flat = similar(intensity, size(intensity, 1) * size(intensity, 2))

    # fill the ψ array
    for i = 1:size(x, 1)
        for j = 1:size(y, 1)
            ψ[i+(j-1)*size(y, 1)] = [x[i] y[j] 0 0]
            int_flat[i+(j-1)*size(y, 1)] = intensity[i, j]
        end
    end

    # return the wanted struct
    return Wave2d(ψ, int_flat, x, y)
end

"""
    getintensity(wave::Wave)::Vector{<:Real}

Return the intensity.

Calculate the current intensity and return it. The `wave` value is the wave struct
from which the intensity should be returned.
"""
function getintensity(wave::Wave1d)::Vector{<:Real}
    intensity = zeros(Float64, size(wave.intensity))
    dx = abs(wave.x[1]-wave.x[2])

    # loop over the incident points
    for i = 1:size(wave.ψ, 1)

        # loop over the endpoints
        for j = 1:size(wave.ψ, 1)

            # check if the
            if abs(wave.ψ[i][1]-wave.x[j]) <= dx/2

                # add the intensity from the incident plane
                # to the output plane
                intensity[j] += wave.intensity[i]

                # if a incident point is assigned to a output intensity
                # break the loop (for the sake of rounding and energy consercation)
                break
            else end
        end
    end

    # return the intensity at the different points
    return intensity
end
