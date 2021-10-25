
struct Lens <: Component
    mat::Matrix{<:Real}
end

"""
    Lens(f::Real)::Lens

Return the Lens type.

Construct and return the lens type. The `f` denotes
the focal length of the lens.
"""
Lens(f::Real)::Lens = Lens([1. 0. 0. 0.; 0. 1. 0. 0.; -1/f 0. 1. 0.; 0. -1/f 0. 1.])


"""
    calculate!(rays::Electron, lens::Lens)

Apply the effect of the lens on the electrons.

This function takes the `Electron1d` object and the `lens` object,
and uses the information provided to apply the effect of the Lens on
the electron beam.
"""
function calculate!(rays::Electron, lens::Lens)
    for ray in rays.ψ
        ray[:] = lens.mat * ray
    end
end
