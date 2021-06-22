
abstract type AbstractElectron end

mutable struct Electron1d <: AbstractElectron
    ψ::Vector{Array}
    intensity::Vector{<:Real}
    p_ges::Vector{<:Real}
    x::Vector{<:Real}
end

"""
    Electron(x::Vector{<:Real}, intensity::Vector{<:Real})::Electron1d

Return a Electron1d type.

Create and return the 1D wave function. The `x` value corresponds to the
x-axis, and `intensity` corresponds to the intensity distribution along
the given axis.
"""
function Electron(x::Vector{<:Real}, intensity::Vector{<:Real}, momentum::Real)::Electron1d
   
    # create the new position angle vectors
    ψ = Vector{Array}(undef, size(x, 1))
    p_ges = similar(intensity)
    for i = 1:size(ψ, 1)
        ψ[i] = [x[i] 0]
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
             momentum::Real)::Electron2d

Return the Electron2d type.

Construct and return a 2D wave. The `x` and `y` vector correspond to the
x and y axis of the problem, and the `intensity` is the intensity corresponding
to the x and y axis given prior.
"""
function Electron(x::Vector{<:Real}, y::Vector{<:Real}, intensity::Matrix{<:Real},
                  momentum::Real)::Electron2d
    
    # create the array with the different angles and positions
    ψ = Vector{Array}(undef, size(x, 1) * size(y, 1))

    # create a flattened intensity array
    int_flat = similar(intensity, size(intensity, 1) * size(intensity, 2))
    p_ges = similar(int_flat)

    # fill the ψ array
    for i = 1:size(x, 1)
        for j = 1:size(y, 1)
            ψ[i+(j-1)*size(y, 1)] = [x[i] y[j] 0 0]
            int_flat[i+(j-1)*size(y, 1)] = intensity[i, j]
            p_ges[i+(j-1)*size(y, 1)] = momentum
        end
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
    getintensity(wave::Electron2d)::Matrix{<:Real}

Return the intensity.

Calculate the corrent intensity and return it. The `wave` value is the wave struct
from which the intensity should be returned.
"""
function getintensity(wave::Electron2d)::Matrix{<:Real}
    x = wave.x
    y = wave.y
    intensity = zeros(Float64, size(x, 1), size(y, 1))
    dx = abs(x[1]-x[2])
    dy = abs(y[1]-y[2])

    # create a breakit variable: break the second loop if breakit is true
    breakit = false
 
    # loop over the rays
    for i = 1:size(wave.ψ, 1)

        # loop over the y axis
        for k = 1:size(y, 1)

            # loop over the x axis
            for j = 1:size(x, 1)

                # check if the ray is near the current coordinate position
                if abs(wave.ψ[i][1]-x[j]) <= dx/2 && abs(wave.ψ[i][2]-y[k]) <= dy/2

                    # add the rays to the intensity
                    intensity[j, k] = wave.intensity[i]

                    # break the loops if a ray is assigned to a new destination
                    breakit = true
                    break
                else end
            end

            # if the ray is assigned, i.e. breakit is true,
            # then break the second loop
            if breakit

                # reset the breakit variable
                breakit = false

                # break the second loop
                break
            else end
        end
    end

    # return the intensity
    return intensity
end
