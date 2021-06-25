
# import the image toolbox
using Images

abstract type AbstractLight end

# create the 1D lightfield strcut (see constructor)
struct LightField1d <: AbstractLight
    intensity::Vector{<:Real}
    envelope::Vector{<:Real}
    x::Vector{<:Real}
    t::Vector{<:Real}
    λ::Real
    E::Real
end

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
    LightField(intensity::Vector{<:Real}, envelope::Vector{<:Real}, x::Vector{<:Real},
               t::Vector{<:Real}, λ::Real, E::Real)::LightField1D

Return the 1D Light Field struct.

Create and return a Light Field struct (1D). The intensity is denoted with the
`intensity` parameter, the temporal envelope function with the `envelope`
parameter, the coordinates for the intensity pattern with the `x` parameter,
and the coordinates of the temporal envelope with the `t` parameter.
The `λ` parameter denotes the wavelength and `E` the pulse energy.
"""
function LightField(intensity::Vector{<:Real}, envelope::Vector{<:Real}, x::Vector{<:Real},
                    t::Vector{<:Real}, λ::Real, E::Real)::LightField1d
    return LightField1d(intensity, envelope, x, t, λ, E)
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

