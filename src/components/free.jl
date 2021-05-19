
# the free propagation struct and constructor
struct Free <: Component
    mat1D::Matrix{<:Real}
    # mat2D::Matrix{<:Real} TODO
end

"""
    Free(d::Real)::Free

Return the Free type.

Construct and return the Free type. The `d` value represents
the distance that the free propagation should be calculated.
"""
Free(d::Real)::Free = Free([1 d; 0 1])

