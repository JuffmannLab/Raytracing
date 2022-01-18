
struct PondInteraction <: Component
    lf::AbstractLight
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

    # calculate the ω parameter
    ω = 2*π*c / lf.λ

    # define the electron Energy and momentum
    γ = 1 / sqrt(1 - (ray.v / c)^2)
    p = γ * m_e * ray.v

    # calculate the constant
    constant = - q^2 * lf.E / (1+ray.v / c) / (2 * m_e * γ * ω_L^2 * ε_0 * c) / lf.norm

    # create the interpolation object
    itp = interpolate(light_intensity, BSpline(Quadratic(Reflect(OnCell()))))
    nx = size(x_imprint, 1)
    ny = size(y_imprint, 1)

    # iterate over the all the electron beams
    for electron in ray.ψ

        # calculate the position that should be interpolated
        xi = (nx - 1) / (x_imprint[end] - x_imprint[1]) * (electron[1] - x_imprint[1]) + 1
        yi = (ny - 1) / (y_imprint[end] - y_imprint[1]) * (electron[2] - y_imprint[1]) + 1

        # calculate the momentum change
        Δpx, Δpy = constant * Interpolations.gradient(itp, xi, yi)

        # calculate the tangens and a constant
        tan_α = tan(electron[3])
        tan_β = tan(electron[4])
        A = sqrt(1 + tan_α^2 + tan_β^2)

        # calculate the momentum
        pz = p / A
        px = pz * tan_α
        py = pz * tan_β
        
        # calculate the new angle
        # the complex - real transformation is to prevent floating point errors
        # due to floating point differences not being 0 when they should be
        electron[3] = atan((px + Δpx) / 
                           real(sqrt(complex(pz^2 * A^2 - (px + Δpx)^2 - (py + Δpy)^2))))
        electron[4] = atan((py + Δpy) / 
                           real(sqrt(complex(pz^2 * A^2 - (px + Δpx)^2 - (py + Δpy)^2))))

    end
end
