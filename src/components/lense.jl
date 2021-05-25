
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
