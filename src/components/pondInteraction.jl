
struct PondInteraction <: Component
    lf::LightField
end

"""
    calculate!(ray::Electron, pond::PondInteraction)

Calculate the scattering angle.

Calculate the angle that the electron has after the interaction with the light field
and save it in the ray object. The input parameters are `pond` which is a light field
object, aswell as `ray` which is the electron beams object.
"""
function calculate!(ray::Electron, pond::PondInteraction)
 
    # get the light field out of the PondInteraction struct
    lf = pond.lf

    # define the electron Energy and momentum
    γ = 1 / sqrt(1 - (ray.v / c)^2)
    β = ray.v / c
    p = γ * m_e * ray.v
    αħ = q^2 / (4 * π * ε_0 * c)
    E_e = γ * m_e * c^2

    # calculate the constant
    constant = - 1 / (1 + β) * αħ / (2 * π) * lf.E / E_e * lf.λ^2 / lf.norm

    gradx, grady = _gradient(lf.intensity, lf.x, lf.y)

    # iterate over the all the electron beams
    for i in eachindex(ray.ψ)

        # calculate the momentum change in x and y direction
        Δpx = constant * _interpolation(gradx, lf.x, lf.y, ray.ψ[i][1], ray.ψ[i][2])
        Δpy = constant * _interpolation(grady, lf.x, lf.y, ray.ψ[i][1], ray.ψ[i][2])

        # calculate the tangens and a constant
        tan_α = tan(ray.ψ[i][3])
        tan_β = tan(ray.ψ[i][4])
        A = sqrt(1 + tan_α^2 + tan_β^2)

        # calculate the momentum
        pz = p / A
        px = pz * tan_α
        py = pz * tan_β
        
        # calculate the new angle
        # the complex - real transformation is to prevent floating point errors
        # due to floating point differences not being 0 when they should be
        ray.ψ[i][3] = atan((px + Δpx) / 
                           real(sqrt(complex(pz^2 * A^2 - (px + Δpx)^2 - (py + Δpy)^2))))
        ray.ψ[i][4] = atan((py + Δpy) / 
                           real(sqrt(complex(pz^2 * A^2 - (px + Δpx)^2 - (py + Δpy)^2))))
    end
end


