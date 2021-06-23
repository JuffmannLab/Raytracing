
# the free propagation struct and constructor
struct Free <: Component
    mat1D::Matrix{<:Real}
    mat2D::Matrix{<:Real}
end

"""
    Free(d::Real)::Free

Return the Free type.

Construct and return the Free type. The `d` value represents
the distance that the free propagation should be calculated.
"""
Free(d::Real)::Free = Free([1 d; 0 1], [1 0 d 0; 0 1 0 d; 0 0 1 0; 0 0 0 1])

"""
    calculate!(rays::Electron1d, free::Free)

Apply the free propagation on the electrons.

The free propagation is applied to all the differen rays that are contained
in the `rays` object, and applies the effect that is saved in the `free` object.
"""
function calculate!(rays::Electron1d, free::Free)
    for ray in rays.ψ
        ray[:] = free.mat1D * ray
    end
end

"""
    calculate!(rays::Electron1d, free::Free)

Apply the free propagation on the electrons.

The free propagation is applied to all the differen rays that are contained
in the `rays` object, and applies the effect that is saved in the `free` object.
"""
function calculate!(ray::Electron2d, free::Free)
    for ray in rays.ψ
        ray[:] = free.mat2D * ray
    end
end
