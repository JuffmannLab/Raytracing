
# import the image toolbox
using Images

"""
    LightField(intensity::Matrix{<:Real}, norm::Real, x::Vector{<:Real},
               y::Vector{<:Real}, λ::Real, E::Real)::LightField

Return a Light Field struct.

Create and return a Light Field struct, where the intensity is given by the
`intensity` parameter, and the normalization factor (integral over the whole space)
`norm` of the intensity, the coordinates in the `x` and `y` direction The `λ` parameter 
denotes the wavelength of the light field and the parameter `E` is the pulse energy.
"""
struct LightField
    intensity::Matrix{<:Real}
    norm::Real
    x::Vector{<:Real}
    y::Vector{<:Real}
    λ::Real
    E::Real
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

