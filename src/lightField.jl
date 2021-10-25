
# import the image toolbox
using Images

abstract type AbstractLight end

# create the 2D lightfield struct (see constructor)
struct LightField2d <: AbstractLight
    intensity::Matrix{<:Real}
    envelope::Vector{<:Real}
    x::Vector{<:Real}
    y::Vector{<:Real}
    t::Vector{<:Real}
    λ::Real
    E::Real
end

"""
    LightField(intensity::Matrix{<:Real}, envelope::Vector{<:Real}, x::Vector{<:Real},
               y::Vector{<:Real}, t::Vector{<:Real}, λ::Real, E::Real)::LightField2D

Return the 2D Light Field struct.

Create and return a 2D Light Field struct, where the intensity is given by the
`intensity` parameter, the temporal structure with the `envelope` parameter,
the coordinates in the `x` and `y` coordinates and the coordinates of the
temporal structure with the `t` parameter. The `λ` parameter denotes
the wavelength of the light field and the parameter `E` is the pulse energy.
"""
function LightField(intensity::Matrix{<:Real}, envelope::Vector{<:Real}, x::Vector{<:Real},
                    y::Vector{<:Real}, t::Vector{<:Real}, λ::Real, E::Real)::LightField2d
    return LightField2d(intensity, envelope, x, y, t, λ, E)
end

"""
    loadintensity(s::String)::Matrix{<:Real}

Return the intensity matrix.

Load the intensity image from the file at `s`. The image will then
be normalized, such that the values range from 0 to 1, and then returned as
a Matrix with the datatype being Float64.
"""
function loadintensity(s::String)::Matrix{<:Real}
    # load the image, save it in a matrix
    input = Float64.(Gray.(load(s)))

    # normalize the image between 0 and 1
    input .-= minimum(input)
    input ./= maximum(input)

    # return the wanted image
    return input
end

