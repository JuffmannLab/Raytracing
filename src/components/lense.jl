
struct Lense <: Component
    mat1D::Matrix{<:Real}
    mat2D::Matrix{<:Real}
end

"""
    Lense(f::Real)::Lense

Return the Lense type.

Construct and return the lense type. The `f` denotes
the focal length of the lense.
"""
Lense(f::Real)::Lense = Lense([1 0; -1/f 1], [1 0 0 0; 0 1 0 0; -1/f 0 1 0; 0 -1/f 0 1])

"""
    calculate!(rays::Electron1d, lense::Lense)

Apply the effect of the lense on the electrons.

This function takes the `Electron1d` object and the `lense` object,
and uses the information provided to apply the effect of the Lense on
the electron beam.
"""
function calculate!(rays::Electron1d, lense::Lense)
    for ray in rays.ψ
        ray[:] = lense.mat1D * ray
    end
end

"""
    calculate!(rays::Electron1d, lense::Lense)

Apply the effect of the lense on the electrons.

This function takes the `Electron1d` object and the `lense` object,
and uses the information provided to apply the effect of the Lense on
the electron beam.
"""
function calculate!(rays::Electron2d, lense::Lense)
    for ray in rays.ψ
        ray[:] = lense.mat2D * ray
    end
end
