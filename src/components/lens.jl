
struct Lens <: Component
    mat1D::Matrix{<:Real}
    mat2D::Matrix{<:Real}
end

"""
    Lens(f::Real)::Lens

Return the Lens type.

Construct and return the lens type. The `f` denotes
the focal length of the lens.
"""
Lens(f::Real)::Lens = Lens([1. 0.; -1/f 1.], 
                           [1. 0. 0. 0.; 0. 1. 0. 0.; -1/f 0. 1. 0.; 0. -1/f 0. 1.])

"""
    calculate!(rays::Electron1d, lens::Lens)

Apply the effect of the lens on the electrons.

This function takes the `Electron1d` object and the `lens` object,
and uses the information provided to apply the effect of the Lens on
the electron beam.
"""
function calculate!(rays::Electron1d, lens::Lens)
    for ray in rays.ψ
        ray[:] = lens.mat1D * ray
    end
end

"""
    calculate!(rays::Electron1d, lens::Lens)

Apply the effect of the lens on the electrons.

This function takes the `Electron1d` object and the `lens` object,
and uses the information provided to apply the effect of the Lens on
the electron beam.
"""
function calculate!(rays::Electron2d, lens::Lens)
    for ray in rays.ψ
        ray[:] = lens.mat2D * ray
    end
end
