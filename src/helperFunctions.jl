"""
    _interpolation(A::Matrix{<:Real}, coords_x::Vector{<:Real},
                   coords_y::Vector{<:Real}, x::Real, y::Real)::Real
Return the bilinear interpolation.
Calculate and return the bilinear interpolation on the Matrix `A`
with the coordinates `coords_x` and `coords_y` at the point with
the position `x`, `y`.
"""
function _interpolation(A::Matrix{<:Real}, coords_x::Vector{<:Real},
                        coords_y::Vector{<:Real}, x::Real, y::Real)::Real

    # if the value is not in the area of coords return 0
    if x < coords_x[1] || x > coords_x[end]
        return zero(eltype(A))
    elseif y < coords_y[1] || y > coords_y[end]
        return zero(eltype(A))
    end

    # calculate the Δx and Δy values
    Δx = abs(coords_x[1]-coords_x[2])
    Δy = abs(coords_y[1]-coords_y[2])

    # calculate the bins in which the x and y values are
    n = round(Int, (x-coords_x[1]) / Δx, RoundDown) + 1
    m = round(Int, (y-coords_y[1]) / Δy, RoundDown) + 1

    # interpolation in x direction on both y-points
    k1 = (A[n, m] - A[n+1, m]) / (coords_x[n] - coords_x[n+1])
    d1 = A[n, m] - k1 * coords_x[n]

    k2 = (A[n, m+1] - A[n+1, m+1]) / (coords_x[n] - coords_x[n+1])
    d2 = A[n, m+1] - k2 * coords_x[n]

    A1 = k1 * x + d1
    A2 = k2 * x + d2

    # use the interpolated values on the x-axis for the interpolation in the y axis
    k = (A1 - A2) / (coords_y[m] - coords_y[m+1])
    d = A1 - k * coords_y[m]

    # return the interpolated value
    return k * y + d

end


"""
    _gradient(V::Matrix{<:AbstractFloat}, x::Vector{<:Real},
              y::Vector{<:Real})::Tuple{Matrix{<:AbstractFloat},Matrix{<:AbstractFloat}}
Calculate the gradient.
This function calculates and returns the gradient of a 2D input potential `V`, and will
return the Gradient in the form of the tuple corresponding to `(Fx, Fy)`, where `Fx` and `Fy`
are Matrices with the same size and datatype as the `V` value. The `x` and `y` values
correspond to the coordinate system of the viewed Potentials/Forces.
"""
function _gradient(V::Matrix{<:AbstractFloat}, x::Vector{<:Real},
                   y::Vector{<:Real})::Tuple{Matrix{<:AbstractFloat},Matrix{<:AbstractFloat}}

    # define the output gradient vector field
    Fx = similar(V)
    Fy = similar(V)

    # set the boundaries
    Fx[1, :] .= 0
    Fy[1, :] .= 0
    Fx[end, :] .= 0
    Fy[end, :] .= 0
    Fx[:, 1] .= 0
    Fy[:, 1] .= 0
    Fx[:, end] .= 0
    Fy[:, end] .= 0

    # calculate the Δx and Δy value
    Δx = abs(x[1]-x[2])
    Δy = abs(y[1]-y[2])

    # calculate the gradient
    for j = 2:size(V, 2)-1
        for i = 2:size(V, 1)-1
            Fx[i, j] = (V[i+1,j] - V[i-1,j]) / 2 / Δx
            Fy[i, j] = (V[i,j+1] - V[i,j-1]) / 2 / Δy
        end
    end

    return (Fx, Fy)
end
