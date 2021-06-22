abstract type Component end

include("./components/free.jl")
include("./components/lense.jl")

struct Setup
    setup::Vector{<:Component}
end

"""
    Setup(comps::Component...)::Setup

Return the Setup struct.

Create and return the setup struct. The varargs `comps` describe
the different components that are used in the setup.
"""
Setup(comps::Component...) = Setup(collect(comps))

"""
    propagation!(wave::Wave1d, setup::Setup)

Propagate the wave in the setup.

This function calculates the matrix that describes the setup `setup` in a
raytracing way. This matrix is then applied to all the rays contained
in `wave`. This function is responsible for the 1D propagation.
"""
function propagation!(wave::Electron1d, setup::Setup)
    
    # create the setup matrix, multiplication with all the
    # setup elements.
    setup_mat = setup.setup[end].mat1D

    # multiply all the elements
    for i = size(setup.setup, 1)-1:-1:1
        setup_mat *= setup.setup[i].mat1D
    end

    # apply the setup matrix onto the wavefunction
    for i = 1:size(wave.ψ, 1)
        wave.ψ[i] *= setup_mat
    end
end


"""
    propagation!(wave::Wave2d, setup::Setup)

Propagate the wave in the setup.

This function calculates the matrix that describes the setup `setup` in a
raytracing way. This matrix is then applied to all the rays contained
in `wave`. This function is responsible for the 2D propagation.
"""
function propagation!(wave::Electron2d, setup::Setup)
    
    setup_mat = setup.setup[end].mat2D

    # multiply all the elements
    for i = size(setup.setup, 1)-1:-1:1
        setup_mat *= setup.setup[i].mat2D
    end

    # apply the setup matrix onto the wavefunction
    for i = 1:size(wave.ψ, 1)
        wave.ψ[i] *= setup_mat
    end
end
